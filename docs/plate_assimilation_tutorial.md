# 4D 板块与俯冲板片数据同化实战教程

← 返回 [README](../README.md) | 完整参数表：[EarthByte 科学特性与参数总表](earthbyte_features.md)

本教程面向地球动力学模拟研究者，系统介绍如何在本引擎中使用吸收自 **EarthByte/CitcomS（Bower et al., 2015）** 的 **4D 板块重构、岩石圈热同化、俯冲板片热/运动学同化、EBA 耗散数缩放、下地幔粘度结构以及混合示踪粒子**功能，并充分发挥现代 C++17 版本独有的 **FGMRES + V-cycle + Block 3×3 求解器、精确地质年代步长对齐、无滞后边界同步与 ParaView 一键可视化（VTK/HDF5）**优势。

---

## 一、物理框架与 `citcoms_cpp` 架构增强

### 1.1 4D 地球动力学数据同化的四大核心模块

在“板块运动历史驱动地幔对流”（Data-Assimilating Mantle Convection）模型中，地表板块重构（如 GPlates / Müller et al. 模型）通过以下四个由浅入深的同化机制耦合到地幔对流求解器：

#### 1. 地表板块速度边界条件（`file_vbcs=1`）
在每一步根据当前地质年代 $t_{\text{age}}$（Ma），对相邻整百万年的地表速度场 $v_\theta, v_\phi$ 进行线性插值并作为顶面 Dirichlet 速度边界条件（`topvbc=1`）。

#### 2. 时变岩石圈热结构同化（`lith_age=1`, `lith_age_time=1`）
读取海底/大陆岩石圈热年龄场 $A(\theta, \phi, t)$（无量纲化后记为 $A_{\text{nd}} = A / \text{scalet}$ ），基于半空间冷却模型计算理论温度 $T_{\text{lith}} = T_m \operatorname{erf}\left(\frac{z}{2\sqrt{A_{\text{nd}}}}\right)$ ，并利用深度平滑混合权重将动力学温度与同化温度平滑过渡：

$$
T_{\text{new}} = \mathrm{daf} \cdot T_{\text{dyn}} + (1 - \mathrm{daf}) \cdot T_{\text{lith}}, \qquad \mathrm{daf} = \left(\frac{z}{z_{\text{lith}}}\right)^p
$$

其中指数 $p$ 由参数 `lith_age_exponent` 控制； $z_{\text{lith}}$ 可设为固定深度（`lith_age_depth_function=off`），或采用随年龄变化的半空间热边界层厚度 $z_{\text{lith}} = 2.32\sqrt{A_{\text{nd}}} \cdot d_0$ （`lith_age_depth_function=on`，其中前因子 $d_0$ 由 `lith_age_depth` 给出）。对于形变板块边界或洋脊等不希望强制同化的区域，只需在年龄文件中写入小于等于 `lith_age_stencil_value` 的数值（如 `-999`），求解器会自动跳过该区域使其自由动力学演化。

#### 3. 俯冲板片热结构同化（`slab_assim=1`）
读取三维俯冲板片温度 $T_{\text{slab}}$ 与连续掩膜权重 $\alpha \in [0, 1]$ ，在每个温度时间步结束时执行平滑热同化：

$$
T_{\text{new}} = (1 - \alpha)\, T_{\text{dyn}} + \alpha\, T_{\text{slab}}
$$

当 $\alpha = 1$ 时强制锁定为俯冲板片理论温度；当 $0 < \alpha < 1$ 时平滑过渡；当 $\alpha = 0$ 时完全自由对流。

#### 4. 俯冲板片内部运动学速度同化（`internal_vbcs_file=1`）
读取三维内部速度场 $(v_\theta, v_\phi, v_r)$ 与整型三态标志 `flag`（`0` = 自由求解，`1` = 强制施加内部速度 Dirichlet 边界，`2` = 保持原有边界标志不变），动态在所有多重网格层级（`levmin..levmax`）更新 `VBX/VBY/VBZ` 标志并重建残差约束与刚度矩阵。

### 1.2 与原生 EarthByte/CitcomS 的架构改进与增强

| 维度 | 原生 EarthByte/CitcomS | 本版本集成增强 |
|---|---|---|
| **内部速度同化时序** | 在时间步末尾更新 `internal_velocity_bc`，导致 Stokes 求解滞后 1 步，且重启首步易失步 | 在 `general_stokes_solver` **之前**（以及 `restart` 读入检查点后）同步更新 `internal_velocity_bc` 与刚度矩阵，**零步滞后** |
| **整百万年文件读取** | 在起始整年代（如 `start_age = 200 Ma`）会尝试打开 `201` 文件导致崩溃；非插值分支遗漏 `/scalet` | 对 `action = 1, 2, 5, 6, 7` 统一采用精确整年代单文件保护，不再越界查找 `max_age + 1`；修正 `age_t` 非插值分支量纲换算 |
| **时间输出与 0 Ma 终止** | 依赖 `storage_spacing_Myr` 过线触发，无法精确踩在整百万年 | `storage_spacing_Myr` 自动映射到 `citcoms_cpp` 的 `storage_interval_myr` + `exact_time_output = on`，**精确踩点整 Myr 与 `0.0 Ma`** |
| **高粘度对比 Stokes 求解** | 传统 `uzawa = cg` + `fmg`，强横向粘度对比（ $10^4\text{–}10^6$ ）及内部速度约束下易震荡满迭代 | 原生兼容 `uzawa = fgmres` + `mg_mode = vcycle` + `mg_smoother = block3x3_gs`，速度提升 **2.4×–2.7×** |
| **EBA 参考态向后兼容** | 强制读取 `refstate_file` 第 5 列作 `Di(z)`，导致传统 7 列（第 5 列为 `0`）文件加热项归零 | 新增显式开关 `eba_dis_scaling`（默认 `off`），100% 兼容旧 4/7 列文件；开启后同时支持 5 列与 7 列 EBA 文件 |
| **可视化与并行 I/O** | 仅传统分进程 ASCII / 旧版 VTK | 完整支持 `ascii,vtk`（`cap`/`single` 模式 + `.pvd` + `comp_nd`）与 `ascii,hdf5`（`.xdmf`）及全部 EarthByte 诊断输出 |

---

## 二、输入文件格式与命名规范

### 2.1 文件命名规则

时间演化的输入文件按**地质年代整数（Ma）**命名（注意：地质年代从过去向现今**递减**，如 `200` → `199` → ... → `0`）：

| 模型几何 | 文件命名模板 | 示例（`start_age = 200`） |
|---|---|---|
| **全球 12-cap (`CitcomSFull`)** | `<prefix><age>.<cap>`，其中 `<cap> = 0..11` | `bvel200.0` … `bvel200.11`，`slab200.0` … `slab0.11` |
| **区域 1-cap (`CitcomSRegional`)** | `<prefix><age>`（无 `.cap` 后缀） | `bvel200` … `bvel0`，`slab200` … `slab0` |

> **注意**：`citcoms_cpp` 会在 `<prefix>` 后直接拼接 `%0.0f`（不自动加点号），因此如果您的文件名为 `bvel.200.0`，请在配置中写 `vel_bound_file = "/path/to/bvel."`（末尾带点）；如果文件名为 `bvel200.0`，则写 `vel_bound_file = "/path/to/bvel"`。

### 2.2 网格节点排列顺序（极其重要！）

CitcomS 中所有 2D（地表）和 3D（全局单 cap）ASCII 输入文件必须严格遵循以下循环顺序（1-based 索引）：

- **2D 地表文件（`bvel`、`age`、`tbc`），总行数 = `nodex × nodey`**：
  ```text
  for j in 1..nodey:        # phi 方向（外层）
      for i in 1..nodex:    # theta 方向（内层，变化最快）
          输出节点 (i, j)
  ```
- **3D 体积文件（`slab`、`ivel`），总行数 = `nodex × nodey × nodez`**：
  ```text
  for j in 1..nodey:            # phi 方向（最外层）
      for i in 1..nodex:        # theta 方向（中层）
          for k in 1..nodez:    # 径向 r：k=1 为 CMB (radius_inner)，k=nodez 为地表 (radius_outer)（变化最快！）
              输出节点 (i, j, k)
  ```

### 2.3 五类输入文件的精确列格式

#### 1. 地表速度边界文件（`vel_bound_file`）
- **维度**：2D（`nodex × nodey` 行）
- **列格式**：`v_theta  v_phi`（2 列浮点数）
- **单位**：**cm/yr**（南向 $v_\theta$ 为正，东向 $v_\phi$ 为正；程序内部自动乘以 `scalev` 和 `timedir` 转为无量纲速度）。

#### 2. 岩石圈年龄文件（`lith_age_file`）
- **维度**：2D（`nodex × nodey` 行）
- **列格式**：`age_Ma`（1 列浮点数）
- **单位**：**Ma**（程序内部自动除以 `scalet` 转为无量纲扩散时间）。
- **掩膜约定**：形变带/俯冲带上方/洋中脊等不参与岩石圈热同化的节点，写入 ≤ `lith_age_stencil_value` 的值（如 `-999.0`）；当两端点插值任一端 ≤ `lith_age_stencil_value` 时，自动保持掩膜不插值。

#### 3. 俯冲板片热同化文件（`slab_assim_file`）
- **维度**：3D（`nodex × nodey × nodez` 行）
- **列格式**：`T_slab  alpha`（2 列浮点数）
  - `T_slab`：无量纲板片理论温度（范围通常在 `[0.0, 1.0]`）。
  - `alpha`：连续热同化掩膜权重，范围 `[0.0, 1.0]`（`0.0` = 无同化，`1.0` = 完全同化为 `T_slab`，`0~1` = 平滑过渡）。

#### 4. 内部速度同化文件（`velocity_internal_file`）
- **维度**：3D（`nodex × nodey × nodez` 行）
- **列格式**：`v_theta  v_phi  v_r  flag`（3 列浮点数 + 1 列整型标志）
  - `v_theta, v_phi, v_r`：单位 **cm/yr**（径向 $v_r$ 向上为正，俯冲下沉为负）。
  - `flag`（三态标志，**至关重要**）：
    - `0`：关闭该节点的内部速度/应力边界标志，恢复自由 Stokes 求解。
    - `1`：开启该节点的 `VBX | VBY | VBZ` 内部速度约束，速度锁定为该行指定的 `(v_theta, v_phi, v_r)`。
    - `2`：**跳过该节点**，保持已有边界标志不变。

> **重要安全规范**：**顶面（`k = nodez`）、底面 CMB（`k = 1`）以及 Regional 模型的四个侧面（`i = 1, nodex` 和 `j = 1, nodey`）必须统一标记为 `flag = 2`！** 否则会覆盖地表板块速度边界或 CMB 自由滑移边界导致求解发散。

#### 5. EBA 参考态文件（`refstate_file`，当 `reference_state = 0`）
- **行数**：`nodez` 行（从 CMB `k = 1` 到地表 `k = nodez`）
- **列格式**：
  - **5 列 EarthByte EBA 格式**（需设 `eba_dis_scaling = on`）：
    `rho(z)  g(z)  alpha(z)  Cp(z)  Di_scale(z)`
  - **7 列兼容格式**：`rho  g  alpha  Cp  col5  col6  col7`
    - 当 `eba_dis_scaling = off`（默认）时忽略第 5 列（等效 `Di_scale = 1.0`），与旧版 CitcomS 100% 数值一致。
    - 当 `eba_dis_scaling = on` 时读取第 5 列作为深度相关耗散数缩放因子 $Di_{\text{scale}}(z)$，自动作用于粘性耗散加热、绝热加热与相变潜热项。
- **归一化要求**：地表处（最后一行 `k = nodez`）各参考态物理量（`rho`, `g`, `alpha`, `Cp`）应归一化为 `1.0`。

---

## 三、数据准备与 GPlates / pyGPlates 工作流

### 3.1 典型 GPlates 前处理管线

对于全球板块重构同化（如 Müller et al., 2016/2019、Young et al., 2019、Bower et al., 2015）：

1. **导出网格坐标**：先运行 0 步（`maxstep = 0` 或 `stokes_flow_only = 1` 且 `output_optional = coord`）导出 CitcomS 各 cap 的节点坐标 $(\theta, \phi, r)$。
2. **使用 pyGPlates / GPlates 采样**：
   - 将每个重建年代（如 `200 Ma` 到 `0 Ma`，步长 `1 Ma`）的拓扑板块速度插值到 12 个 cap 的 2D 地表网格，输出 `bvel<age>.<cap>`。
   - 将对应的洋壳年龄格网（加上大陆岩石圈等效热年龄，如 150–300 Ma，以及形变带 `-999` 掩膜）采样到 2D 网格，输出 `age<age>.<cap>`。
   - 依据俯冲带位置与俯冲倾角/汇聚速率，构造三维半空间/平板热异常与平滑过渡掩膜，输出 `slab<age>.<cap>`（如需运动学导引，同时输出 `ivel<age>.<cap>`）。
3. **运行前完整性校验脚本（推荐）**：

```python
#!/usr/bin/env python3
"""快速校验 4D 同化输入文件的行数、范围与边界 flag=2 合规性"""
import os

def verify_assim_files(prefix_dir, ages, caps, nox, noy, noz):
    n2d, n3d = nox * noy, nox * noy * noz
    for age in ages:
        for cap in caps:
            suffix = f"{age}.{cap}" if cap is not None else f"{age}"
            # 1. 校验 2D 速度与年龄文件
            for name in ("bvel", "age"):
                p = os.path.join(prefix_dir, f"{name}{suffix}")
                assert os.path.exists(p), f"缺失文件: {p}"
                with open(p) as f:
                    n = sum(1 for _ in f)
                assert n == n2d, f"{p} 行数错误: 期望 {n2d}, 实际 {n}"
            # 2. 校验 3D 俯冲板片与内部速度文件
            pslab = os.path.join(prefix_dir, f"slab{suffix}")
            if os.path.exists(pslab):
                with open(pslab) as f:
                    lines = f.readlines()
                assert len(lines) == n3d, f"{pslab} 行数错误: {len(lines)} != {n3d}"
            pivel = os.path.join(prefix_dir, f"ivel{suffix}")
            if os.path.exists(pivel):
                with open(pivel) as f:
                    lines = f.readlines()
                assert len(lines) == n3d, f"{pivel} 行数错误: {len(lines)} != {n3d}"
                # 检查顶面 (k=noz) 与底面 (k=1) 必须全部为 flag=2
                for j in range(noy):
                    for i in range(nox):
                        for k in (0, noz - 1):
                            idx = k + i * noz + j * nox * noz
                            flag = int(float(lines[idx].split()[3]))
                            assert flag == 2, f"{pivel} 节点 (i={i+1},j={j+1},k={k+1}) 边界 flag 必须为 2，实际为 {flag}"
    print("所有 4D 同化输入文件格式与边界保护校验通过！")
```

---

## 四、开箱即用配置模板

### 4.1 模板 A：全球 12-cap 完整 4D 同化 + EBA + 现代高速求解器（推荐生产配置）

适合使用 `mpirun -np 96 ./build_linux/CitcomSFull global_4d_assim.cfg`（12 caps × 2×2×2 分块）运行：

```ini
# ============================================================
# citcoms_cpp 全球 4D 板块/俯冲同化 + EBA + FGMRES/V-cycle 模板
# ============================================================

# --- 1. 输出目录、格式与可视化（仅路径/前缀参数使用 ""）---
datadir="./output_global_4d"
datafile="earth4d"
output_format=ascii,vtk
output_optional_ascii=horiz_avg,divv,sten_temp,sten_velo
output_optional_vtk=comp_nd,surf
vtk_format=binary
vtk_output_mode=cap                  # 12 cap 合并 .vts + .pvd，ParaView 一键播放

# --- 2. 网格与 MPI 分解（12 caps × 2×2×2 = 96 ranks）---
nproc_surf = 12
nprocx = 2
nprocy = 2
nprocz = 2
nodex = 65
nodey = 65
nodez = 65
mgunitx = 4
mgunity = 4
mgunitz = 4
levels = 4                           # 4*2^(4-1)*2 + 1 = 65
radius_inner = 0.55
radius_outer = 1.0
coor = 1                             # 径向加密（可选）
coor_refine = 0.1,0.15,0.1,0.2

# --- 3. 地质年代控制与无损断点续算 ---
start_age = 200.0                    # 起始地质年代：200 Ma
storage_spacing_Myr = 5              # EarthByte 兼容参数：每 5 Myr 输出（自动启用 exact_time_output 精确踩点）
exact_time_output = on               # 自动裁剪 dt 精确落到 195, 190, ..., 0 Ma
exit_at_present = on                 # 到达 0.0 Ma 自动输出并安全退出
checkpoint_interval_myr = 10.0       # 每 10 Myr 存盘检查点
checkpoints_to_keep = 3              # 滑动窗口保留最近 3 个检查点
checkpoint_on_termination = on
restart = auto                       # 队列作业自动续算，无需手动改 cfg
minstep = 1
maxstep = 100000
maxtotstep = 100000
cpu_limits_in_seconds = 86000

# --- 4. 4D 板块、岩石圈、俯冲板片与内部速度同化 ---
topvbc = 1                           # 顶面指定板块速度
botvbc = 0                           # CMB 自由滑移
file_vbcs = 1
vel_bound_file = "./inputs/bvel/bvel"

lith_age = 1
lith_age_time = 1
lith_age_file = "./inputs/age/age"
lith_age_depth_function = on         # 2.32 * sqrt(age) * lith_age_depth
lith_age_depth = 1.0                 # 开启 depth_function 后作为无量纲前因子（1.0 = 标准半空间深度）
lith_age_exponent = 0.5              # EarthByte 经典平方根混合衰减（推荐 0.5 或 1.0）
lith_age_min = 0.01                  # 防止洋中脊 age=0 奇异（Ma）
lith_age_stencil_value = -100.0      # age <= -100 Ma（如 -999）的区域跳过岩石圈热同化

slab_assim = on
slab_assim_file = "./inputs/slab/slab"

internal_vbcs_file = on
velocity_internal_file = "./inputs/ivel/ivel"

# --- 5. EBA 参考态、粘性/绝热加热与相变 ---
rayleigh = 7.8e7
dissipation_number = 1.56            # 地表参考耗散数 Di_0
gruneisen = 1.0
surfaceT = 0.1
reference_state = 0
refstate_file = "./inputs/refstate_eba5col.dat"
eba_dis_scaling = on                 # 启用第 5 列深度相关 Di(z) 缩放
Q0 = 20.0                            # 内部放射性生热

# --- 6. EarthByte 下地幔粘度结构（Flament/Steinberger-Calderwood）---
Viscosity = system
VISC_UPDATE = on
num_mat = 4
z_layer = 0.025,0.06435,0.10516,0.45 # 约 160 / 410 / 670 / 2867 km
TDEPV = on
rheol = 22                           # 下地幔 N0 随深度线性 10 -> 100（也可试 20, 21, 25, 27）
visc0 = 1.0,1.0,1.0,10.0
viscE = 5.0,5.0,5.0,5.0
viscT = 0.5,0.5,0.5,0.5
VMIN = on
visc_min = 0.1
VMAX = on
visc_max = 1000.0

# --- 7. 混合示踪粒子（例如追踪 CMB LLSVP 化学致密层）---
tracer = on
tracer_ic_method = 0
tracers_per_element = 20
tracer_flavors = 2
ic_method_for_flavors = 0
z_interface = 0.425                  # 底部约 160 km 致密化学层
chemical_buoyancy = on
buoy_type = 1
hybrid_method = on                   # 空单元自动视为背景地幔 C=0，无需全地幔密铺粒子
buoyancy_ratio = 0.35
tracer_cushion = 100000
itracer_warnings = off               # 长时演化中若偶发个别丢失粒子不强制中止

# --- 8. citcoms_cpp 现代高速 Stokes + 多重网格求解器 ---
Solver = multigrid
aug_lagr = on
aug_number = 2.0e3
uzawa = fgmres                       # 强烈推荐：彻底消除高粘度对比 CG 停滞
fgmres_restart = 30
inner_rel_tol = 0.05                 # 相对内迭代残差
mg_mode = vcycle                     # 纯 V-cycle
mg_smoother = block3x3_gs            # 耦合 3x3 节点块平滑器
coarse_solver = cg                   # 最粗网格 PCG
visc_smooth_method = 4               # 几何平均粘度粗化（保留俯冲弱边界）
accuracy = 1.0e-3
piterations = 200
vhighstep = 3
vlowstep = 15

# --- 9. 温度边界与初始条件 ---
toptbc = 1
toptbcval = 0.0
bottbc = 1
bottbcval = 1.0
tic_method = 0
mantle_temp = 0.5
ADVECTION = on
finetunedt = 0.75
```

### 4.2 模板 B：纯地表速度 + 岩石圈热同化（无内部速度约束）

若您希望仅通过地表板块速度和岩石圈热年龄驱动俯冲，而不人为干预下沉板片的速度场，只需在模板 A 中设置：
```ini
slab_assim = off
internal_vbcs_file = off
lith_age = 1
lith_age_time = 1
lith_age_depth_function = on
lith_age_depth = 1.0
lith_age_exponent = 0.5
```

### 4.3 模板 C：扣除浅部岩石圈浮力的动态地形计算

在计算全球动态地形（Dynamic Topography）时，常需排除最浅部岩石圈（例如深度 0–300 km 对应的全局径向节点以上）的浮力贡献。在 `post_p = 1` 或 `stokes_flow_only = 1` 的后处理配置中加入：
```ini
# 假设 nodez = 65，地表为 global znode = 65，约 300 km 对应 znode = 56：
remove_buoyancy_above_znode = 56     # 将 global z-node > 56 的热/化学浮力全部置零
output_optional = surf,botm,geoid
```

---

## 五、关键参数选型参考与排错指南

### 5.1 `lith_age_exponent` 混合曲线对比

| `lith_age_exponent` 值 | 混合公式 $\mathrm{daf}(z/z_{\text{lith}})$ | 物理效果与适用场景 |
|---|---|---|
| 未设置（默认 `-1.0`） | $0.5 \cdot (z/z_{\text{lith}})$（传统 CIG 模式） | 100% 兼容老版本 CitcomS 算例 |
| `0.0` | $\mathrm{daf} = 1.0$ | 关闭岩石圈内部热混合（仅用 `lith_age` 初始化或只做 `slab_assim`） |
| `0.5`（**EarthByte 推荐**） | $\sqrt{z/z_{\text{lith}}}$ | 浅部快速贴近理论冷却曲线，深部温和过渡至自由地幔 |
| `1.0` | $z/z_{\text{lith}}$ | 线性过渡（浅部 100% 同化，底界 0% 同化） |
| `2.0` 或更大 | $(z/z_{\text{lith}})^2$ | 更强力地锁定整个岩石圈内部温度接近半空间冷却解 |

### 5.2 EarthByte 粘度结构（`rheol`）速查表

| `rheol` | 公式与下地幔（`l = 3`）前因子 $N_0(z)$ 特征 | 参考文献 |
|---|---|---|
| `20` | $\eta = N_0 \exp\big(\frac{E}{T+T_0} - \frac{E}{0.5+T_0}\big)$ ， $T\in[0,1]$ ，参考背景温度 $0.5$ | Flament et al. (2013, 2014) |
| `21` | 同 `20`，下地幔 $N_0(z)$ 从 `10` 线性增至 `50`（670–2867 km） | Flament et al. (2014) TC8 |
| `22` | 同 `20`，下地幔 $N_0(z)$ 从 `10` 线性增至 `100`（670–2867 km） | Flament et al. (2014) TC7 / Zhang et al. (2010) |
| `23` | 同 `20`，下地幔 $N_0(z)$ 从 `10` 线性增至 `30`（670–2867 km） | Flament et al. (2014) TC9 |
| `24` | Frank-Kamenetskii $\eta = N_0(z)\exp(E(T_0-T))$ ，下地幔自 `2` 线性递增 | Zhang et al. (2010) / Bower et al. (2013) |
| `25` | 同 `20`，下地幔 $N_0(z)$ 从 `5` 线性增至 `100`（EBA 专用） | Flament EBA 系列 |
| `26` | 同 `20`，下地幔 $N_0(z)$ 从 `5` 线性增至 `80`（EBA 专用） | Flament EBA 系列 |
| `27` | 同 `20`，下地幔 $N_0(z)$ 在 **1000–2867 km** 间从 `10` 线性增至 `100` | Rudolph et al. (2015) / Flament EBA |
| `28` | 同 `20`，下地幔 $N_0(z)$ 从 `2` 线性增至 `100` | Flament EBA 系列 |
| `112` | 下地幔指数增长 $N_0(z) = 0.5 \times 10^{4z/(r_o-r_i)-0.5}$ | Hassan et al. |
| `113` | 下地幔深度相关活化能项 + 底部热边界层粘度裁剪 `[20, 100]`（ $z > 0.43$ ） | Hassan et al. |
| `117` | Arrhenius（`rheol=7` 变体）， $T\in[0,1]$ 且浅部 $z < 0.06$ 限制 $\eta \le 1$ | Hassan et al. |
| `118` | 同 `117`，无浅部粘度裁剪 | Hassan et al. |

### 5.3 常见问题排查（FAQ）

1. **Q：开启 `internal_vbcs_file = on` 后 Stokes 求解不收敛或散度偏大？**
   - **检查 1**：确认所有顶面（`k = nodez`）、底面（`k = 1`）及 Regional 侧面的 `flag` 均为 `2`，不能为 `0` 或 `1`！
   - **检查 2**：确保内部指定速度场近似满足不可压缩连续性方程（ $\nabla \cdot \mathbf{v} \approx 0$ ），避免在相邻节点上强行给定强烈辐合/辐散的速度。
   - **检查 3**：在 `output_optional_ascii` 中开启 `divv,sten_velo` 检查散度高值区是否正好位于 `flag = 1` 边缘；同时开启 `uzawa = fgmres` + `mg_smoother = block3x3_gs` + `aug_lagr = on`。

2. **Q：开启 `eba_dis_scaling = on` 后报错提示列数不足，或者加热项为 0？**
   - 检查 `refstate_file` 是否至少有 5 列，且第 5 列为正的耗散数深度缩放系数（如 `1.0` 左右，而非旧文件留空的 `0.0`）。若使用传统 7 列（第 5 列填 0）参考态文件，请保持 `eba_dis_scaling = off`。

3. **Q：在 `output_format=ascii,vtk` 下找不到 `.sten_temp` 或 `.divv` 文件？**
   - 复合输出模式下，ASCII 诊断字段由 `output_optional_ascii` 单独控制，请在配置中显式加入：
     `output_optional_ascii=horiz_avg,sten_temp,sten_velo,divv`。
   - 在单独 `output_format=ascii` 模式下，直接写在 `output_optional` 中即可生效。

4. **Q：如何验证所有 EarthByte 功能在本地编译环境中工作正常？**
   - 运行内置自动化测试套件：
     ```bash
     python3 tests/test_earthbyte_features.py
     ```
     该脚本会在临时目录中依次验证 EBA `Di` 缩放与旧 7 列兼容性、全部 13 种 EarthByte 粘度结构、完整 4D 岩石圈/板片/内部速度同化与重启、动态地形浮力截断、混合示踪粒子以及 12 核全球球谐诊断输出。

---

## 六、来源声明与学术引用

本教程涉及的 4D 板块/俯冲板片数据同化、EBA 耗散数深度缩放、下地幔粘度结构（`rheol=20–28, 112–118`）、混合示踪粒子及诊断输出机制吸收改编自 **EarthByte/CitcomS**（`https://github.com/EarthByte/citcoms`，GPL-2.0 开源协议）。使用上述功能发表科研成果时，请在引用原始 CitcomS 文献（Moresi et al., 2014; Zhong et al., 2000, 2008; Tan et al., 2006）的同时，按需引用以下 EarthByte 核心文献：

- Bower, D. J., Gurnis, M., & Flament, N. (2015). *Assimilating lithosphere and slab history in 4-D Earth models*. Physics of the Earth and Planetary Interiors, 238, 8–22. [https://doi.org/10.1016/j.pepi.2014.10.013](https://doi.org/10.1016/j.pepi.2014.10.013)
- Flament, N., Gurnis, M., Williams, S., Seton, M., Skogseid, J., Heine, C., & Müller, R. D. (2014). *Topographic asymmetry of the South Atlantic from global models of mantle flow and lithospheric stretching*. Earth and Planetary Science Letters, 387, 107–119. [https://doi.org/10.1016/j.epsl.2013.11.017](https://doi.org/10.1016/j.epsl.2013.11.017)
- Zhang, N., Zhong, S., Leng, W., & Li, Z.-X. (2010). *A model for the evolution of the Earth's mantle structure since the Early Paleozoic*. JGR: Solid Earth, 115, B06401. [https://doi.org/10.1029/2009JB006896](https://doi.org/10.1029/2009JB006896)
- Rudolph, M. L., Lekić, V., & Lithgow-Bertelloni, C. (2015). *Viscosity jump in Earth's mid-mantle*. Science, 350(6266), 1349–1352. [https://doi.org/10.1126/science.aad1929](https://doi.org/10.1126/science.aad1929)
