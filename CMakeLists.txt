cmake_minimum_required(VERSION 3.16)
project(CitcomS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags for GCC and Clang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare)
endif()

# MPI is required for CitcomS
find_package(MPI REQUIRED)

# Optional packages
find_package(OpenMP)

option(ENABLE_HDF5 "Enable HDF5 parallel I/O support" OFF)
if(ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C)
    if(HDF5_FOUND)
        message(STATUS "HDF5 found. Enabling HDF5 output support.")
        add_definitions(-DUSE_HDF5)
        include_directories(${HDF5_INCLUDE_DIRS})
    else()
        message(WARNING "HDF5 requested but not found.")
    endif()
endif()

# Core CitcomS library sources
set(CITCOMS_LIB_SOURCES
    lib/Advection_diffusion.cpp
    lib/Anisotropic_viscosity.cpp
    lib/BC_util.cpp
    lib/Checkpoints.cpp
    lib/Citcom_init.cpp
    lib/Composition_related.cpp
    lib/Construct_arrays.cpp
    lib/Convection.cpp
    lib/Determine_net_rotation.cpp
    lib/Drive_solvers.cpp
    lib/Element_calculations.cpp
    lib/General_matrix_functions.cpp
    lib/Ggrd_handling.cpp
    lib/Global_operations.cpp
    lib/Initial_temperature.cpp
    lib/Instructions.cpp
    lib/Interuption.cpp
    lib/Lith_age.cpp
    lib/Material_properties.cpp
    lib/Mineral_physics_models.cpp
    lib/Nodal_mesh.cpp
    lib/Output.cpp
    lib/Output_gzdir.cpp
    lib/Output_h5.cpp
    lib/Output_vtk.cpp
    lib/Pan_problem_misc_functions.cpp
    lib/Parallel_util.cpp
    lib/Parsing.cpp
    lib/Phase_change.cpp
    lib/Problem_related.cpp
    lib/Process_buoyancy.cpp
    lib/Shape_functions.cpp
    lib/Size_does_matter.cpp
    lib/Solver_conj_grad.cpp
    lib/Solver_multigrid.cpp
    lib/Sphere_harmonics.cpp
    lib/Sphere_util.cpp
    lib/Stokes_flow_Incomp.cpp
    lib/Topo_gravity.cpp
    lib/Tracer_setup.cpp
    lib/Viscosity_structures.cpp
    lib/Full_boundary_conditions.cpp
    lib/Full_geometry_cartesian.cpp
    lib/Full_lith_age_read_files.cpp
    lib/Full_parallel_related.cpp
    lib/Full_read_input_from_files.cpp
    lib/Full_solver.cpp
    lib/Full_sphere_related.cpp
    lib/Full_tracer_advection.cpp
    lib/Full_version_dependent.cpp
    lib/Regional_boundary_conditions.cpp
    lib/Regional_geometry_cartesian.cpp
    lib/Regional_lith_age_read_files.cpp
    lib/Regional_parallel_related.cpp
    lib/Regional_read_input_from_files.cpp
    lib/Regional_solver.cpp
    lib/Regional_sphere_related.cpp
    lib/Regional_tracer_advection.cpp
    lib/Regional_version_dependent.cpp
)

add_library(CitcomS STATIC ${CITCOMS_LIB_SOURCES})
target_include_directories(CitcomS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(CitcomS PUBLIC MPI::MPI_CXX)
if(ENABLE_HDF5 AND HDF5_FOUND)
    target_link_libraries(CitcomS PUBLIC ${HDF5_LIBRARIES})
endif()

# Executables
# 1. CitcomSRegional
add_executable(CitcomSRegional bin/Citcom.cpp bin/CitcomSRegional.cpp)
target_link_libraries(CitcomSRegional PRIVATE CitcomS MPI::MPI_CXX)

# 2. CitcomSFull
add_executable(CitcomSFull bin/Citcom.cpp bin/CitcomSFull.cpp)
target_link_libraries(CitcomSFull PRIVATE CitcomS MPI::MPI_CXX)

# 3. Py2C (Standalone tool)
add_executable(Py2C Py2C/Py2CMain.cpp Py2C/Py2C.cpp)
target_include_directories(Py2C PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Py2C)

# 4. citcoms_combine (Modern high-performance C++17 combiner)
add_executable(citcoms_combine bin/citcoms_combine.cpp)
if(OpenMP_FOUND)
    target_link_libraries(citcoms_combine PRIVATE OpenMP::OpenMP_CXX)
endif()

