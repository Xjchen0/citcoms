/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <algorithm>
#include <vector>
#include "element_definitions.h"
#include "global_defs.h"

namespace {
inline double* ensure_mg_buf(std::vector<double>& vec, size_t n) {
    if (vec.size() < n) {
        vec.resize(n);
    }
    return vec.data();
}
} // namespace

#ifndef USE_CUDA
static void block3x3_gauss_seidel(struct All_variables *E, double **d0, double **F, double **Ad,
                                  double acc, int *cycles, int level, int guess,
                                  int symmetric, int start_dir);
static void chebyshev_smoother(struct All_variables *E, double **d0, double **F, double **Ad,
                               double acc, int *cycles, int level, int guess);
#endif

static void mg_smooth(struct All_variables *E, double **d0, double **F, double **Ad,
                      double acc, int *cycles, int level, int guess, int start_dir)
{
#ifndef USE_CUDA
    if (strcmp(E->control.mg_smoother, "block3x3_gs") == 0) {
        block3x3_gauss_seidel(E, d0, F, Ad, acc, cycles, level, guess, 0, 0);
        return;
    }
    if (strcmp(E->control.mg_smoother, "block3x3_sgs") == 0) {
        block3x3_gauss_seidel(E, d0, F, Ad, acc, cycles, level, guess, 1, start_dir);
        return;
    }
    if (strcmp(E->control.mg_smoother, "chebyshev") == 0) {
        chebyshev_smoother(E, d0, F, Ad, acc, cycles, level, guess);
        return;
    }
#endif
    (void)start_dir;
    gauss_seidel(E, d0, F, Ad, acc, cycles, level, guess);
}

static void mg_solve_coarse(struct All_variables *E, double **vel_low, double **rhs_low,
                            double **AU_low, double acc, int levmin)
{
    int cycles = E->control.v_steps_low;
    if (strcmp(E->control.coarse_solver, "cg") == 0) {
        const double r0 = sqrt(global_vdot(E, rhs_low, rhs_low, levmin));
        if (r0 <= 1.0e-30) {
            const int neq = E->lmesh.NEQ[levmin];
            for (int m = 1; m <= E->sphere.caps_per_proc; m++)
                for (int i = 0; i < neq; i++)
                    vel_low[m][i] = 0.0;
            return;
        }
        const double coarse_acc = (acc * 0.01 < 1.0e-3 * r0) ? (acc * 0.01) : (1.0e-3 * r0);
        conj_grad(E, vel_low, rhs_low, coarse_acc, &cycles, levmin);
    } else {
        mg_smooth(E, vel_low, rhs_low, AU_low, acc * 0.01, &cycles, levmin, 0, 0);
    }
}

#ifdef _UNICOS
#include <fortran.h>
#endif

int epsilon[4][4] = {   /* Levi-Cita epsilon */
  {0, 0, 0, 0},
  {0, 1,-1, 1},
  {0,-1, 1,-1},
  {0, 1,-1, 1} };


/*  ===========================================================
    Iterative solver also using multigrid  ........
    ===========================================================  */

int solve_del2_u(struct All_variables *E, double **d0, double **F, double acc, int high_lev)
{
  /* modern c++: void assemble_del2_u() prototype in header */
  /* modern c++: void e_assemble_del2_u() prototype in header */
  /* modern c++: void n_assemble_del2_u() prototype in header */
  /* modern c++: void strip_bcs_from_residual() prototype in header */
  /* modern c++: void gauss_seidel() prototype in header */

  /* modern c++: double conj_grad() prototype in header */
  /* modern c++: double multi_grid() prototype in header */
  /* modern c++: double global_vdot() prototype in header */
  /* modern c++: void record() prototype in header */
  /* modern c++: void report() prototype in header */

  int count,counts,cycles,convergent,valid;
  int i, neq, m;

  char message[200];

  double initial_time, time; /* funcs moved to header: CPU_time0() */
  double residual,prior_residual,r0;
  double *D1[NCS], *r[NCS], *Au[NCS];

  neq  = E->lmesh.NEQ[high_lev];

  for (m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=0;i<neq;i++)  {
	d0[m][i] = 0.0;
      }

  r0=residual=sqrt(global_vdot(E,F,F,high_lev));

  /* Optional relative inner tolerance: prevents early inaccurate inner exits when ||grad(s)|| is small */
  double target_acc = acc;
  if (E->control.inner_rel_tol > 0.0) {
      const double rel_acc = E->control.inner_rel_tol * r0;
      if (rel_acc < target_acc)
          target_acc = rel_acc;
  }

  prior_residual=2*residual;
  count = 0;
  initial_time=CPU_time0();

  if (!E->control.NMULTIGRID) {
    /* conjugate gradient solution */

    cycles = E->control.v_steps_low;
    residual = conj_grad(E,d0,F,target_acc,&cycles,high_lev);
    valid = (residual < target_acc)? 1:0;
  } else  {
    
    /* solve using multigrid  */
    const int max_cycles = (E->control.inner_mg_cycles > 0)
                           ? E->control.inner_mg_cycles
                           : E->control.max_mg_cycles;

    counts =0;
    if(E->parallel.me==0){	/* output */
      snprintf(message,200,"resi = %.6e for iter %d acc %.6e",residual,counts,target_acc);
      record(E,message);
      report(E,message);
    }

    do {
      residual=multi_grid(E,d0,F,target_acc,high_lev);
      valid = (residual < target_acc || E->control.inner_mg_cycles > 0)?1:0;
      counts ++;
      if(E->parallel.me==0){	/* output  */
	snprintf(message,200,"resi = %.6e for iter %d acc %.6e",residual,counts,target_acc);
	record(E,message);
	report(E,message);
      }
    }  while ((residual >= target_acc) && counts < max_cycles);

    cycles = counts;
  }


  /* Convergence check .....
     We should give it a chance to recover if it briefly diverges initially, and
     don't worry about slower convergence if it is close to the answer   */

  if(((count > 0) && (residual > r0*2.0))  ||
     ((fabs(residual-prior_residual) < acc*0.1) && (residual > acc * 10.0))   )
    convergent=0;
  else {
    convergent=1;
    prior_residual=residual;
  }

  if(E->control.print_convergence&&E->parallel.me==0)   {
    fprintf(E->fp,"%s residual (%03d) = %.3e from %.3e to %.3e in %5.2f secs \n",
	    (convergent ? " * ":"!!!"),cycles,residual,r0,acc,CPU_time0()-initial_time);
    fflush(E->fp);
  }

  count++;

  E->monitor.momentum_residual = residual;
  E->control.total_iteration_cycles += count;
  E->control.total_v_solver_calls += 1;

  return(valid);
}

/* =================================
   recursive multigrid function ....
   ================================= */

double multi_grid(struct All_variables *E, double **d1, double **F, double acc, int hl) /* higher level of two */
{
    double residual,AudotAu;
    /* modern c++: void interp_vector() prototype in header */
    /* modern c++: void project_vector() prototype in header */
    int m,i,j,Vn,Vnmax,cycles;
    double alpha,beta;
    /* modern c++: void gauss_seidel() prototype in header */
    /* modern c++: void element_gauss_seidel() prototype in header */
    /* modern c++: void e_assemble_del2_u() prototype in header */
    /* modern c++: void strip_bcs_from_residual() prototype in header */
    /* modern c++: void n_assemble_del2_u() prototype in header */

    /* modern c++: conj_grad(), global_vdot() prototypes in header */

    FILE *fp;
    char filename[1000];
    int lev,ic,ulev,dlev;

    const int levmin = E->mesh.levmin;
    const int levmax = E->mesh.levmax;

    double time1; /* funcs moved to header: CPU_time0() */
    double *res[MAX_LEVELS][NCS],*AU[MAX_LEVELS][NCS];
    double *vel[MAX_LEVELS][NCS],*del_vel[MAX_LEVELS][NCS];
    double *rhs[MAX_LEVELS][NCS],*fl[MAX_LEVELS][NCS];
				/* because it's recursive, need a copy at
				    each level */

    static std::vector<double> ws_del_vel[MAX_LEVELS][NCS], ws_AU[MAX_LEVELS][NCS];
    static std::vector<double> ws_vel[MAX_LEVELS][NCS], ws_res[MAX_LEVELS][NCS], ws_fl[MAX_LEVELS][NCS];

    for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)
      for(m=1;m<=E->sphere.caps_per_proc;m++)    {
	del_vel[i][m] = ensure_mg_buf(ws_del_vel[i][m], (size_t)(E->lmesh.NEQ[i]+1));
	AU[i][m]      = ensure_mg_buf(ws_AU[i][m],      (size_t)(E->lmesh.NEQ[i]+1));
	vel[i][m]     = ensure_mg_buf(ws_vel[i][m],     (size_t)(E->lmesh.NEQ[i]+1));
	res[i][m]     = ensure_mg_buf(ws_res[i][m],     (size_t)(E->lmesh.NEQ[i]+1));
	if (i<E->mesh.levmax)
	  fl[i][m]    = ensure_mg_buf(ws_fl[i][m],      (size_t)(E->lmesh.NEQ[i]+1));
      }

    Vnmax = E->control.mg_cycle;

    if (strcmp(E->control.mg_mode, "vcycle") == 0) {
        /* Standard V-cycle: single fine-to-coarse-to-fine traversal without FMG warm-up.
         * Cuts coarse-grid MPI latency by ~4x and saves fine-grid MatVecs when used as a
         * preconditioner inside FGMRES / FCG. */
        for (m = 1; m <= E->sphere.caps_per_proc; m++)
            for (j = 0; j < E->lmesh.NEQ[levmax]; j++)
                res[levmax][m][j] = F[m][j];

        for (Vn = 1; Vn <= Vnmax; Vn++) {
            /* Downward stroke of the V */
            for (dlev = levmax; dlev >= levmin + 1; dlev--) {
                cycles = ((dlev == levmax) ? E->control.v_steps_high : E->control.down_heavy);
                mg_smooth(E, vel[dlev], res[dlev], AU[dlev], 0.01, &cycles, dlev, 0, 0);

                for (m = 1; m <= E->sphere.caps_per_proc; m++)
                    for (i = 0; i < E->lmesh.NEQ[dlev]; i++)
                        res[dlev][m][i] -= AU[dlev][m][i];

                project_vector(E, dlev, res[dlev], res[dlev - 1], 1);
                strip_bcs_from_residual(E, res[dlev - 1], dlev - 1);
            }

            /* Bottom of the V */
            mg_solve_coarse(E, vel[levmin], res[levmin], AU[levmin], acc, levmin);

            /* Upward stroke of the V */
            for (ulev = levmin + 1; ulev <= levmax; ulev++) {
                cycles = ((ulev == levmax) ? E->control.v_steps_high : E->control.up_heavy);

                interp_vector(E, ulev - 1, vel[ulev - 1], del_vel[ulev]);
                strip_bcs_from_residual(E, del_vel[ulev], ulev);
                mg_smooth(E, del_vel[ulev], res[ulev], AU[ulev], 0.01, &cycles, ulev, 1, 1);

                double AudotRes = 0.0;
                global_vdot2(E, AU[ulev], AU[ulev], res[ulev], ulev, &AudotAu, &AudotRes);
                alpha = (AudotAu != 0.0) ? (AudotRes / AudotAu) : 0.0;

                for (m = 1; m <= E->sphere.caps_per_proc; m++)
                    for (i = 0; i < E->lmesh.NEQ[ulev]; i++)
                        vel[ulev][m][i] += alpha * del_vel[ulev][m][i];

                if (ulev == levmax) {
                    for (m = 1; m <= E->sphere.caps_per_proc; m++)
                        for (i = 0; i < E->lmesh.NEQ[ulev]; i++)
                            res[ulev][m][i] -= alpha * AU[ulev][m][i];
                }
            }

            for (m = 1; m <= E->sphere.caps_per_proc; m++)
                for (j = 0; j < E->lmesh.NEQ[levmax]; j++)
                    d1[m][j] += vel[levmax][m][j];
        }

        for (m = 1; m <= E->sphere.caps_per_proc; m++)
            for (j = 0; j < E->lmesh.NEQ[levmax]; j++)
                F[m][j] = res[levmax][m][j];

        residual = sqrt(global_vdot(E, F, F, hl));
        return residual;
    }

        /* Full Multigrid (FMG) - legacy default */
        /* Project residual onto all the lower levels */

    project_vector(E,levmax,F,fl[levmax-1],1);
    strip_bcs_from_residual(E,fl[levmax-1],levmax-1);
    for(lev=levmax-1;lev>levmin;lev--) {
        project_vector(E,lev,fl[lev],fl[lev-1],1);
        strip_bcs_from_residual(E,fl[lev-1],lev-1);
       }

        /* Solve for the lowest level */

/*    time=CPU_time0(); */
    mg_solve_coarse(E,vel[levmin],fl[levmin],AU[levmin],acc,levmin);

    for(lev=levmin+1;lev<=levmax;lev++) {
      //time=CPU_time0();

                         /* Utilize coarse solution and smooth at this level */
      interp_vector(E,lev-1,vel[lev-1],vel[lev]);
      strip_bcs_from_residual(E,vel[lev],lev);

      if (lev==levmax)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(j=0;j<E->lmesh.NEQ[lev];j++)
             res[lev][m][j]=F[m][j];
      else
        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(j=0;j<E->lmesh.NEQ[lev];j++)
             res[lev][m][j]=fl[lev][m][j];

      for(Vn=1;Vn<=Vnmax;Vn++)   {
                                        /*    Downward stoke of the V    */
        for (dlev=lev;dlev>=levmin+1;dlev--)   {

                                      /* Pre-smoothing  */
          cycles=((dlev==levmax)?E->control.v_steps_high:E->control.down_heavy);
          ic = ((dlev==lev)?1:0);
          mg_smooth(E,vel[dlev],res[dlev],AU[dlev],0.01,&cycles,dlev,ic,0);

          for(m=1;m<=E->sphere.caps_per_proc;m++)
             for(i=0;i<E->lmesh.NEQ[dlev];i++)  {
               res[dlev][m][i]  = res[dlev][m][i] - AU[dlev][m][i];
	       }

          project_vector(E,dlev,res[dlev],res[dlev-1],1);
          strip_bcs_from_residual(E,res[dlev-1],dlev-1);
          }

                                        /*    Bottom of the V    */
       mg_solve_coarse(E,vel[levmin],res[levmin],AU[levmin],acc,levmin);
                                        /*    Upward stoke of the V    */
        for (ulev=levmin+1;ulev<=lev;ulev++)   {
            cycles=((ulev==levmax)?E->control.v_steps_high:E->control.up_heavy);

            interp_vector(E,ulev-1,vel[ulev-1],del_vel[ulev]);
            strip_bcs_from_residual(E,del_vel[ulev],ulev);
            mg_smooth(E,del_vel[ulev],res[ulev],AU[ulev],0.01,&cycles,ulev,1,1);

            double AudotRes = 0.0;
            global_vdot2(E, AU[ulev], AU[ulev], res[ulev], ulev, &AudotAu, &AudotRes);
            alpha = (AudotAu != 0.0) ? (AudotRes / AudotAu) : 0.0;

            for(m=1;m<=E->sphere.caps_per_proc;m++)
              for(i=0;i<E->lmesh.NEQ[ulev];i++)   {
               vel[ulev][m][i] += alpha*del_vel[ulev][m][i];
               }

            if (ulev ==levmax)
              for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(i=0;i<E->lmesh.NEQ[ulev];i++)   {
                  res[ulev][m][i] -= alpha*AU[ulev][m][i];
                  }

            }
        }
      }

    for(m=1;m<=E->sphere.caps_per_proc;m++){
      for(j=0;j<E->lmesh.NEQ[levmax];j++)   {
        F[m][j]=res[levmax][m][j];
        d1[m][j]+=vel[levmax][m][j];
      }
    }
    residual = sqrt(global_vdot(E,F,F,hl));
    
    /* Persistent multigrid workspace retained across calls */
    return(residual);
}


/*  ===========================================================
    Conjugate gradient relaxation for the matrix equation Kd = f
    Returns the residual reduction after itn iterations ...
    ===========================================================  */


#ifndef USE_CUDA
double conj_grad(struct All_variables *E, double **d0, double **F, double acc, int *cycles, int level)
{
    double *r0[NCS],*r1[NCS],*r2[NCS];
    double *z0[NCS],*z1[NCS],*z2[NCS];
    double *p1[NCS],*p2[NCS];
    double *Ap[NCS];
    double *shuffle[NCS];

    int m,count,i,steps;
    double residual;
    double alpha,beta,dotprod,dotr1z1,dotr0z0;

    double time; /* funcs moved to header: CPU_time0() */

    /* modern c++: void parallel_process_termination() prototype in header */
    /* modern c++: void assemble_del2_u() prototype in header */
    /* modern c++: void strip_bcs_from_residual() prototype in header */
    /* modern c++: double global_vdot() prototype in header */

    const int mem_lev=E->mesh.levmax;
    const int high_neq = E->lmesh.NEQ[level];

    steps = *cycles;

    static std::vector<double> ws_r0[NCS], ws_r1[NCS], ws_r2[NCS], ws_z0[NCS], ws_z1[NCS];
    static std::vector<double> ws_p1[NCS], ws_p2[NCS], ws_Ap[NCS];
    const size_t alloc_neq = (size_t)(high_neq + 1);

    for(m=1;m<=E->sphere.caps_per_proc;m++)    {
      r0[m] = ensure_mg_buf(ws_r0[m], alloc_neq);
      r1[m] = ensure_mg_buf(ws_r1[m], alloc_neq);
      r2[m] = ensure_mg_buf(ws_r2[m], alloc_neq);
      z0[m] = ensure_mg_buf(ws_z0[m], alloc_neq);
      z1[m] = ensure_mg_buf(ws_z1[m], alloc_neq);
      p1[m] = ensure_mg_buf(ws_p1[m], alloc_neq);
      p2[m] = ensure_mg_buf(ws_p2[m], alloc_neq);
      Ap[m] = ensure_mg_buf(ws_Ap[m], alloc_neq);
    }

    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=0;i<high_neq;i++) {
	r1[m][i] = F[m][i];
        d0[m][i] = 0.0;
	}

    residual = sqrt(global_vdot(E,r1,r1,level));

    if (residual <= 1.0e-30) {
        *cycles = 0;
        return 0.0;
    }
    count = 0;

    while (((residual > acc) && (count < steps)) || count == 0)  {

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
      int blk_nno = 0;
      const double *blk_m = get_block3x3_BI_ptr(level, m, &blk_nno);
      if (blk_m && blk_nno == E->lmesh.NNO[level]) {
          const struct ID *ID_m = E->ID[level][m];
          for (i = 1; i <= blk_nno; i++) {
              const int e1 = ID_m[i].doff[1];
              const int e2 = ID_m[i].doff[2];
              const int e3 = ID_m[i].doff[3];
              const double v1 = r1[m][e1];
              const double v2 = r1[m][e2];
              const double v3 = r1[m][e3];
              const double *b = blk_m + (size_t)i * 6;
              z1[m][e1] = b[0] * v1 + b[1] * v2 + b[2] * v3;
              z1[m][e2] = b[1] * v1 + b[3] * v2 + b[4] * v3;
              z1[m][e3] = b[2] * v1 + b[4] * v2 + b[5] * v3;
          }
      } else {
          for(i=0;i<high_neq;i++)
             z1[m][i] = E->BI[level][m][i] * r1[m][i];
      }
    }

	dotr1z1 = global_vdot(E,r1,z1,level);

	if (count == 0 )
          for(m=1;m<=E->sphere.caps_per_proc;m++)
	    for(i=0;i<high_neq;i++)
	      p2[m][i] = z1[m][i];
	else {
	    if (dotr0z0 == 0.0) break;
	    beta = dotr1z1/dotr0z0;
            for(m=1;m<=E->sphere.caps_per_proc;m++)
	      for(i=0;i<high_neq;i++)
		p2[m][i] = z1[m][i] + beta * p1[m][i];
	}

    dotr0z0 = dotr1z1;

	assemble_del2_u(E,p2,Ap,level,1);

	dotprod=global_vdot(E,p2,Ap,level);

	if(0.0==dotprod)
	    alpha=1.0e-3;
	else
	    alpha = dotr1z1/dotprod;

        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(i=0;i<high_neq;i++) {
	    d0[m][i] += alpha * p2[m][i];
	    r2[m][i] = r1[m][i] - alpha * Ap[m][i];
	    }

	residual = sqrt(global_vdot(E,r2,r2,level));

        for(m=1;m<=E->sphere.caps_per_proc;m++)    {
	  shuffle[m] = r0[m]; r0[m] = r1[m]; r1[m] = r2[m]; r2[m] = shuffle[m];
	  shuffle[m] = z0[m]; z0[m] = z1[m]; z1[m] = shuffle[m];
	  shuffle[m] = p1[m]; p1[m] = p2[m]; p2[m] = shuffle[m];
          }

	count++;
	/* end of while-loop */

    }

    *cycles=count;

    strip_bcs_from_residual(E,d0,level);

    /* Persistent CG workspace retained across calls */
    return(residual);   }
#endif /* !USE_CUDA */


/* ========================================================================================
   An element by element version of the gauss-seidel routine. Initially this is a test
   platform, we want to know if it handles discontinuities any better than the node/equation
   versions
   =========================================================================================*/

void element_gauss_seidel(struct All_variables *E, double **d0, double **F, double **Ad, double acc, int *cycles, int level, int guess)
{
    int count,i,j,k,l,m,ns,nc,d,steps,loc;
    int p1,p2,p3,q1,q2,q3;
    int e,eq,node,node1;
    int element,eqn1,eqn2,eqn3,eqn11,eqn12,eqn13;

    /* modern c++: void e_assemble_del2_u() prototype in header */
    /* modern c++: void n_assemble_del2_u() prototype in header */
    /* modern c++: void strip_bcs_from_residual() prototype in header */

    double U1[24],AD1[24],F1[24];
    double w1,w2,w3;
    double w11,w12,w13;
    double w[24];

    double *dd[NCS],*elt_k;
    int *vis[NCS];

    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[E->mesh.nsd];
    const int neq=E->lmesh.NEQ[level];
    const int nel=E->lmesh.NEL[level];
    const int nno=E->lmesh.NNO[level];


    steps=*cycles;

    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      dd[m] = (double *)malloc(neq*sizeof(double));
      vis[m] = (int *)malloc((nno+1)*sizeof(int));
    }
    elt_k=(double *)malloc((24*24)*sizeof(double));

    if(guess){
	e_assemble_del2_u(E,d0,Ad,level,1);
    }
    else {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
	for(i=0;i<neq;i++)
	     Ad[m][i]=d0[m][i]=0.0;
    }

   count=0;
   while (count <= steps)      {
     for (m=1;m<=E->sphere.caps_per_proc;m++) {
	for(i=1;i<=nno;i++)
	    vis[m][i]=0;

	for(e=1;e<=nel;e++) {

	    elt_k = E->elt_k[level][m][e].k;

	    for(i=1;i<=ends;i++)  {
		node=E->IEN[level][m][e].node[i];
		p1=(i-1)*dims;
		w[p1] = w[p1+1] = w[p1+2] = 1.0;
		if(E->NODE[level][m][node] & VBX)
		    w[p1] = 0.0;
		if(E->NODE[level][m][node] & VBY)
		    w[p1+1] = 0.0;
		if(E->NODE[level][m][node] & VBZ)
		    w[p1+2] = 0.0;

	        }


	    for(i=1;i<=ends;i++)   {
		node=E->IEN[level][m][e].node[i];
		if(!vis[m][node])
		    continue;

		eqn1=E->ID[level][m][node].doff[1];
		eqn2=E->ID[level][m][node].doff[2];
                eqn3=E->ID[level][m][node].doff[3];
		p1=(i-1)*dims*n;
		p2=p1+n;
		p3=p2+n;


		/* update Au */
		for(j=1;j<=ends;j++) {
		    node1=E->IEN[level][m][e].node[j];

	            eqn11=E->ID[level][m][node1].doff[1];
		    eqn12=E->ID[level][m][node1].doff[2];
		    eqn13=E->ID[level][m][node1].doff[3];
		    q1=(j-1)*3;

                    Ad[m][eqn11] += w[q1]*(elt_k[p1+q1] * dd[m][eqn1] + elt_k[p2+q1] * dd[m][eqn2] + elt_k[p3+q1] * dd[m][eqn3]);
                    Ad[m][eqn12] += w[q1+1]*(elt_k[p1+q1+1] * dd[m][eqn1] + elt_k[p2+q1+1] * dd[m][eqn2] + elt_k[p3+q1+1] * dd[m][eqn3]);
	            Ad[m][eqn13] += w[q1+2]*(elt_k[p1+q1+2] * dd[m][eqn1] + elt_k[p2+q1+2] * dd[m][eqn2] + elt_k[p3+q1+2] * dd[m][eqn3]);

		   }
	        }


	    for(i=1;i<=ends;i++)   {
		node=E->IEN[level][m][e].node[i];
		if(vis[m][node])
		    continue;

		eqn1=E->ID[level][m][node].doff[1];
		eqn2=E->ID[level][m][node].doff[2];
		eqn3=E->ID[level][m][node].doff[3];
		p1=(i-1)*dims*n;
		p2=p1+n;
		p3=p2+n;

		/* update dd, d0 */
		d0[m][eqn1] += (dd[m][eqn1] = w[(i-1)*dims]*(F[m][eqn1]-Ad[m][eqn1])*E->BI[level][m][eqn1]);
		d0[m][eqn2] += (dd[m][eqn2] = w[(i-1)*dims+1]*(F[m][eqn2]-Ad[m][eqn2])*E->BI[level][m][eqn2]);
		d0[m][eqn3] += (dd[m][eqn3] = w[(i-1)*dims+2]*(F[m][eqn3]-Ad[m][eqn3])*E->BI[level][m][eqn3]);

		vis[m][node]=1;

		/* update Au */
		for(j=1;j<=ends;j++) {
		   node1=E->IEN[level][m][e].node[j];

		   eqn11=E->ID[level][m][node1].doff[1];
		   eqn12=E->ID[level][m][node1].doff[2];
		   eqn13=E->ID[level][m][node1].doff[3];
		   q1=(j-1)*3;
		   q2=q1+1;
		   q3=q1+2;

		   Ad[m][eqn11] += w[q1]*(elt_k[p1+q1] * dd[m][eqn1] + elt_k[p2+q1] * dd[m][eqn2] + elt_k[p3+q1] * dd[m][eqn3]);
		   Ad[m][eqn12] += w[q2]*(elt_k[p1+q2] * dd[m][eqn1] + elt_k[p2+q2] * dd[m][eqn2] + elt_k[p3+q2] * dd[m][eqn3]);
		   Ad[m][eqn13] += w[q3]*(elt_k[p1+q3] * dd[m][eqn1] + elt_k[p2+q3] * dd[m][eqn2] + elt_k[p3+q3] * dd[m][eqn3]);
		   }
	       }

	    }          /* end for el  */
	 }               /* end for m */

        (E->solver.exchange_id_d)(E, Ad, level);
        (E->solver.exchange_id_d)(E, d0, level);

	/* completed cycle */

	    count++;

    }

    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      free((double*) dd[m]);
      free((int*) vis[m]);
    }
    free((double*) elt_k);

    return;
}


#ifndef USE_CUDA

/* ============================================================================
   Multigrid Gauss-Seidel relaxation scheme which requires the storage of local
   information, otherwise some other method is required. NOTE this is a bit worse
   than real gauss-seidel because it relaxes all the equations for a node at one
   time (Jacobi at a node). It does the job though.
   ============================================================================ */

void gauss_seidel(struct All_variables *E, double **d0, double **F, double **Ad, double acc, int *cycles, int level, int guess)
{

    int count,i,j,k,l,m,ns,steps;
    int *C;
    int eqn1,eqn2,eqn3;

    /* modern c++: void parallel_process_termination() prototype in header */
    /* modern c++: void n_assemble_del2_u() prototype in header */

    double U1,U2,U3,UU;
    double residual; /* funcs moved to header: global_vdot() */

    higher_precision *B1,*B2,*B3;


    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[E->mesh.nsd];
    const int neq=E->lmesh.NEQ[level];
    const int num_nodes=E->lmesh.NNO[level];
    const int nox=E->lmesh.NOX[level];
    const int noz=E->lmesh.NOY[level];
    const int noy=E->lmesh.NOZ[level];
    const int max_eqn=14*dims;

    const double zeroo = 0.0;

    steps=*cycles;
    //sor = 1.3;

    if(guess) {
      n_assemble_del2_u(E,d0,Ad,level,1);
    }
    else
      for (m=1;m<=E->sphere.caps_per_proc;m++)
	for(i=0;i<neq;i++) {
	    d0[m][i]=Ad[m][i]=zeroo;
	}

    count = 0;


    while (count < steps) {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(j=0;j<=E->lmesh.NEQ[level];j++)
          E->temp[m][j] = zeroo;

      for (m=1;m<=E->sphere.caps_per_proc;m++)
          Ad[m][neq] = zeroo;

      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(i=1;i<=E->lmesh.NNO[level];i++)
          if(E->NODE[level][m][i] & OFFSIDE)   {

	    eqn1=E->ID[level][m][i].doff[1];
	    eqn2=E->ID[level][m][i].doff[2];
	    eqn3=E->ID[level][m][i].doff[3];
	    E->temp[m][eqn1] = (F[m][eqn1] - Ad[m][eqn1])*E->BI[level][m][eqn1];
	    E->temp[m][eqn2] = (F[m][eqn2] - Ad[m][eqn2])*E->BI[level][m][eqn2];
	    E->temp[m][eqn3] = (F[m][eqn3] - Ad[m][eqn3])*E->BI[level][m][eqn3];
	    E->temp1[m][eqn1] = Ad[m][eqn1];
	    E->temp1[m][eqn2] = Ad[m][eqn2];
	    E->temp1[m][eqn3] = Ad[m][eqn3];
            }

      for (m=1;m<=E->sphere.caps_per_proc;m++) {
        double * __restrict__ Ad_m = Ad[m];
        double * __restrict__ temp_m = E->temp[m];
        double * __restrict__ d0_m = d0[m];
        const double * __restrict__ F_m = F[m];
        const double * __restrict__ BI_m = E->BI[level][m];
        const unsigned int * __restrict__ NODE_m = E->NODE[level][m];
        const struct ID * __restrict__ ID_m = E->ID[level][m];
        const int * __restrict__ nmap_base = E->Node_map[level][m];
        const higher_precision * __restrict__ ek1_base = E->Eqn_k1[level][m];
        const higher_precision * __restrict__ ek2_base = E->Eqn_k2[level][m];
        const higher_precision * __restrict__ ek3_base = E->Eqn_k3[level][m];
        const int nno_lev = E->lmesh.NNO[level];

        for(i=1;i<=nno_lev;i++)     {

	    eqn1=ID_m[i].doff[1];
	    eqn2=ID_m[i].doff[2];
	    eqn3=ID_m[i].doff[3];
            const size_t row_off = (size_t)(i-1)*max_eqn;
            const int * __restrict__ C_r = nmap_base + row_off;
	    const higher_precision * __restrict__ B1_r = ek1_base + row_off;
	    const higher_precision * __restrict__ B2_r = ek2_base + row_off;
 	    const higher_precision * __restrict__ B3_r = ek3_base + row_off;

                 /* Ad on boundaries differs after the following operation, but
                  no communications are needed yet, because boundary Ad will
                  not be used for the G-S iterations for interior nodes */

            double a1 = Ad_m[eqn1];
            double a2 = Ad_m[eqn2];
            double a3 = Ad_m[eqn3];
            for(j=3;j<max_eqn;j++)  {
                 const double uu = temp_m[C_r[j]];
                 a1 += (double)B1_r[j]*uu;
                 a2 += (double)B2_r[j]*uu;
                 a3 += (double)B3_r[j]*uu;
                 }
            Ad_m[eqn1] = a1;
            Ad_m[eqn2] = a2;
            Ad_m[eqn3] = a3;

            double t1, t2, t3;
            if (!(NODE_m[i]&OFFSIDE))   {
               t1 = (F_m[eqn1] - a1)*BI_m[eqn1];
               t2 = (F_m[eqn2] - a2)*BI_m[eqn2];
               t3 = (F_m[eqn3] - a3)*BI_m[eqn3];
               temp_m[eqn1] = t1;
               temp_m[eqn2] = t2;
               temp_m[eqn3] = t3;
	       } else {
               t1 = temp_m[eqn1];
               t2 = temp_m[eqn2];
               t3 = temp_m[eqn3];
               }

                 /* Ad on boundaries differs after the following operation */
	    for(j=0;j<max_eqn;j++)
		    Ad_m[C_r[j]] += (double)B1_r[j]*t1
                                 +  (double)B2_r[j]*t2
                                 +  (double)B3_r[j]*t3;

	    d0_m[eqn1] += t1;
	    d0_m[eqn2] += t2;
	    d0_m[eqn3] += t3;
  	    }
      }

      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(i=1;i<=E->lmesh.NNO[level];i++)
          if(E->NODE[level][m][i] & OFFSIDE)   {
	    eqn1=E->ID[level][m][i].doff[1];
	    eqn2=E->ID[level][m][i].doff[2];
	    eqn3=E->ID[level][m][i].doff[3];
	    Ad[m][eqn1] -= E->temp1[m][eqn1];
	    Ad[m][eqn2] -= E->temp1[m][eqn2];
	    Ad[m][eqn3] -= E->temp1[m][eqn3];
	    }

      (E->solver.exchange_id_d)(E, Ad, level);

      for (m=1;m<=E->sphere.caps_per_proc;m++){
 	for(i=1;i<=E->lmesh.NNO[level];i++)
          if(E->NODE[level][m][i] & OFFSIDE)   {
	    eqn1=E->ID[level][m][i].doff[1];
	    eqn2=E->ID[level][m][i].doff[2];
	    eqn3=E->ID[level][m][i].doff[3];
	    Ad[m][eqn1] += E->temp1[m][eqn1];
	    Ad[m][eqn2] += E->temp1[m][eqn2];
	    Ad[m][eqn3] += E->temp1[m][eqn3];
	    }
      }

      count++;

/*     for (m=1;m<=E->sphere.caps_per_proc;m++)
	for(i=0;i<neq;i++)          {
	   F[m][i] -= Ad[m][i];
	   Ad[m][i] = 0.0;
	  }
*/
      }

    *cycles=count;
    return;

}

/* ============================================================================
 * 3x3 Coupled Nodal Block Gauss-Seidel / Symmetric Gauss-Seidel Smoother.
 * Solves the exact 3x3 nodal velocity block (u_theta, u_phi, u_r) simultaneously
 * on all interior nodes, eliminating intra-node component decoupling under strong
 * augmented-Lagrangian penalty and high viscosity contrast.
 * When symmetric != 0, alternates forward (1..nno) and backward (nno..1) sweeps
 * starting with start_dir (0 = forward first, 1 = backward first) so pre- and
 * post-smoothing are exact adjoints.
 * ============================================================================ */
static void block3x3_gauss_seidel(struct All_variables *E, double **d0, double **F, double **Ad,
                                  double acc, int *cycles, int level, int guess,
                                  int symmetric, int start_dir)
{
    (void)acc;
    const int dims = E->mesh.nsd;
    const int neq = E->lmesh.NEQ[level];
    const int nno_lev = E->lmesh.NNO[level];
    const int max_eqn = 14 * dims;
    const int steps = *cycles;
    const double zeroo = 0.0;

    if (guess) {
        n_assemble_del2_u(E, d0, Ad, level, 1);
    } else {
        for (int m = 1; m <= E->sphere.caps_per_proc; m++)
            for (int i = 0; i < neq; i++)
                d0[m][i] = Ad[m][i] = zeroo;
    }

    int count = 0;
    while (count < steps) {
        for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (int j = 0; j <= neq; j++)
                E->temp[m][j] = zeroo;
            Ad[m][neq] = zeroo;
        }

        for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (int i = 1; i <= nno_lev; i++) {
                if (E->NODE[level][m][i] & OFFSIDE) {
                    const int eqn1 = E->ID[level][m][i].doff[1];
                    const int eqn2 = E->ID[level][m][i].doff[2];
                    const int eqn3 = E->ID[level][m][i].doff[3];
                    E->temp[m][eqn1] = (F[m][eqn1] - Ad[m][eqn1]) * E->BI[level][m][eqn1];
                    E->temp[m][eqn2] = (F[m][eqn2] - Ad[m][eqn2]) * E->BI[level][m][eqn2];
                    E->temp[m][eqn3] = (F[m][eqn3] - Ad[m][eqn3]) * E->BI[level][m][eqn3];
                    E->temp1[m][eqn1] = Ad[m][eqn1];
                    E->temp1[m][eqn2] = Ad[m][eqn2];
                    E->temp1[m][eqn3] = Ad[m][eqn3];
                }
            }
        }

        const int is_backward = (symmetric != 0) ? ((count + start_dir) & 1) : 0;

        for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
            double * __restrict__ Ad_m = Ad[m];
            double * __restrict__ temp_m = E->temp[m];
            double * __restrict__ d0_m = d0[m];
            const double * __restrict__ F_m = F[m];
            const unsigned int * __restrict__ NODE_m = E->NODE[level][m];
            const struct ID * __restrict__ ID_m = E->ID[level][m];
            const int * __restrict__ nmap_base = E->Node_map[level][m];
            const higher_precision * __restrict__ ek1_base = E->Eqn_k1[level][m];
            const higher_precision * __restrict__ ek2_base = E->Eqn_k2[level][m];
            const higher_precision * __restrict__ ek3_base = E->Eqn_k3[level][m];
            int blk_nno = 0;
            const double * __restrict__ blk_m = get_block3x3_BI_ptr(level, m, &blk_nno);
            const double * __restrict__ BI_m = E->BI[level][m];

            if (!is_backward) {
                /* Forward sweep: i = 1 .. nno_lev */
                for (int i = 1; i <= nno_lev; i++) {
                    const int eqn1 = ID_m[i].doff[1];
                    const int eqn2 = ID_m[i].doff[2];
                    const int eqn3 = ID_m[i].doff[3];
                    const size_t row_off = (size_t)(i - 1) * max_eqn;
                    const int * __restrict__ C_r = nmap_base + row_off;
                    const higher_precision * __restrict__ B1_r = ek1_base + row_off;
                    const higher_precision * __restrict__ B2_r = ek2_base + row_off;
                    const higher_precision * __restrict__ B3_r = ek3_base + row_off;

                    double a1 = Ad_m[eqn1];
                    double a2 = Ad_m[eqn2];
                    double a3 = Ad_m[eqn3];
                    for (int j = 3; j < max_eqn; j++) {
                        const double uu = temp_m[C_r[j]];
                        a1 += (double)B1_r[j] * uu;
                        a2 += (double)B2_r[j] * uu;
                        a3 += (double)B3_r[j] * uu;
                    }
                    Ad_m[eqn1] = a1;
                    Ad_m[eqn2] = a2;
                    Ad_m[eqn3] = a3;

                    double t1, t2, t3;
                    if (!(NODE_m[i] & OFFSIDE)) {
                        const double r1 = F_m[eqn1] - a1;
                        const double r2 = F_m[eqn2] - a2;
                        const double r3 = F_m[eqn3] - a3;
                        if (blk_m && blk_nno == nno_lev) {
                            const double *b = blk_m + (size_t)i * 6;
                            t1 = b[0] * r1 + b[1] * r2 + b[2] * r3;
                            t2 = b[1] * r1 + b[3] * r2 + b[4] * r3;
                            t3 = b[2] * r1 + b[4] * r2 + b[5] * r3;
                        } else {
                            t1 = r1 * BI_m[eqn1];
                            t2 = r2 * BI_m[eqn2];
                            t3 = r3 * BI_m[eqn3];
                        }
                        temp_m[eqn1] = t1;
                        temp_m[eqn2] = t2;
                        temp_m[eqn3] = t3;
                    } else {
                        t1 = temp_m[eqn1];
                        t2 = temp_m[eqn2];
                        t3 = temp_m[eqn3];
                    }

                    for (int j = 0; j < max_eqn; j++) {
                        Ad_m[C_r[j]] += (double)B1_r[j] * t1
                                      + (double)B2_r[j] * t2
                                      + (double)B3_r[j] * t3;
                    }

                    d0_m[eqn1] += t1;
                    d0_m[eqn2] += t2;
                    d0_m[eqn3] += t3;
                }
            } else {
                /* Backward sweep: i = nno_lev .. 1.
                 * As i decreases, every i' > i has already scattered L^T * t_{i'}
                 * directly into Ad_m[eqn(i)], so Ad_m[i] holds (K*d0 + L^T*t)_i! */
                for (int i = nno_lev; i >= 1; i--) {
                    const int eqn1 = ID_m[i].doff[1];
                    const int eqn2 = ID_m[i].doff[2];
                    const int eqn3 = ID_m[i].doff[3];
                    const size_t row_off = (size_t)(i - 1) * max_eqn;
                    const int * __restrict__ C_r = nmap_base + row_off;
                    const higher_precision * __restrict__ B1_r = ek1_base + row_off;
                    const higher_precision * __restrict__ B2_r = ek2_base + row_off;
                    const higher_precision * __restrict__ B3_r = ek3_base + row_off;

                    double t1, t2, t3;
                    if (!(NODE_m[i] & OFFSIDE)) {
                        const double r1 = F_m[eqn1] - Ad_m[eqn1];
                        const double r2 = F_m[eqn2] - Ad_m[eqn2];
                        const double r3 = F_m[eqn3] - Ad_m[eqn3];
                        if (blk_m && blk_nno == nno_lev) {
                            const double *b = blk_m + (size_t)i * 6;
                            t1 = b[0] * r1 + b[1] * r2 + b[2] * r3;
                            t2 = b[1] * r1 + b[3] * r2 + b[4] * r3;
                            t3 = b[2] * r1 + b[4] * r2 + b[5] * r3;
                        } else {
                            t1 = r1 * BI_m[eqn1];
                            t2 = r2 * BI_m[eqn2];
                            t3 = r3 * BI_m[eqn3];
                        }
                        temp_m[eqn1] = t1;
                        temp_m[eqn2] = t2;
                        temp_m[eqn3] = t3;
                    } else {
                        t1 = temp_m[eqn1];
                        t2 = temp_m[eqn2];
                        t3 = temp_m[eqn3];
                    }

                    for (int j = 0; j < max_eqn; j++) {
                        Ad_m[C_r[j]] += (double)B1_r[j] * t1
                                      + (double)B2_r[j] * t2
                                      + (double)B3_r[j] * t3;
                    }

                    d0_m[eqn1] += t1;
                    d0_m[eqn2] += t2;
                    d0_m[eqn3] += t3;
                }

                /* Complete Ad += L * temp so Ad == K_loc * d0 before boundary exchange */
                for (int i = 1; i <= nno_lev; i++) {
                    const int eqn1 = ID_m[i].doff[1];
                    const int eqn2 = ID_m[i].doff[2];
                    const int eqn3 = ID_m[i].doff[3];
                    const size_t row_off = (size_t)(i - 1) * max_eqn;
                    const int * __restrict__ C_r = nmap_base + row_off;
                    const higher_precision * __restrict__ B1_r = ek1_base + row_off;
                    const higher_precision * __restrict__ B2_r = ek2_base + row_off;
                    const higher_precision * __restrict__ B3_r = ek3_base + row_off;

                    double a1 = Ad_m[eqn1];
                    double a2 = Ad_m[eqn2];
                    double a3 = Ad_m[eqn3];
                    for (int j = 3; j < max_eqn; j++) {
                        const double uu = temp_m[C_r[j]];
                        a1 += (double)B1_r[j] * uu;
                        a2 += (double)B2_r[j] * uu;
                        a3 += (double)B3_r[j] * uu;
                    }
                    Ad_m[eqn1] = a1;
                    Ad_m[eqn2] = a2;
                    Ad_m[eqn3] = a3;
                }
            }
        }

        for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (int i = 1; i <= nno_lev; i++) {
                if (E->NODE[level][m][i] & OFFSIDE) {
                    const int eqn1 = E->ID[level][m][i].doff[1];
                    const int eqn2 = E->ID[level][m][i].doff[2];
                    const int eqn3 = E->ID[level][m][i].doff[3];
                    Ad[m][eqn1] -= E->temp1[m][eqn1];
                    Ad[m][eqn2] -= E->temp1[m][eqn2];
                    Ad[m][eqn3] -= E->temp1[m][eqn3];
                }
            }
        }

        (E->solver.exchange_id_d)(E, Ad, level);

        for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (int i = 1; i <= nno_lev; i++) {
                if (E->NODE[level][m][i] & OFFSIDE) {
                    const int eqn1 = E->ID[level][m][i].doff[1];
                    const int eqn2 = E->ID[level][m][i].doff[2];
                    const int eqn3 = E->ID[level][m][i].doff[3];
                    Ad[m][eqn1] += E->temp1[m][eqn1];
                    Ad[m][eqn2] += E->temp1[m][eqn2];
                    Ad[m][eqn3] += E->temp1[m][eqn3];
                }
            }
        }

        count++;
    }

    *cycles = count;
}

/* ============================================================================
 * 1st-Kind Chebyshev Polynomial Smoother preconditioned by exact global 3x3
 * nodal block Jacobi (M = block_diag(K)).
 * - Inherently symmetric and MPI-rank-count invariant (identical numerics on
 *   any number of ranks, zero MPI_Allreduce inside smoothing loops).
 * - Spectral radius lambda_max(M^{-1} K) is estimated once per stiffness
 *   rebuild via 8 power iterations and cached per level.
 * ============================================================================ */
static void apply_block3x3_inv(struct All_variables *E, int level, double **in, double **out)
{
    const int nno_lev = E->lmesh.NNO[level];
    const int neq = E->lmesh.NEQ[level];
    for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
        int blk_nno = 0;
        const double * __restrict__ blk_m = get_block3x3_BI_ptr(level, m, &blk_nno);
        const double * __restrict__ in_m = in[m];
        double * __restrict__ out_m = out[m];
        if (blk_m && blk_nno == nno_lev) {
            const struct ID * __restrict__ ID_m = E->ID[level][m];
            for (int i = 1; i <= nno_lev; i++) {
                const int e1 = ID_m[i].doff[1];
                const int e2 = ID_m[i].doff[2];
                const int e3 = ID_m[i].doff[3];
                const double r1 = in_m[e1];
                const double r2 = in_m[e2];
                const double r3 = in_m[e3];
                const double *b = blk_m + (size_t)i * 6;
                out_m[e1] = b[0] * r1 + b[1] * r2 + b[2] * r3;
                out_m[e2] = b[1] * r1 + b[3] * r2 + b[4] * r3;
                out_m[e3] = b[2] * r1 + b[4] * r2 + b[5] * r3;
            }
        } else {
            const double * __restrict__ BI_m = E->BI[level][m];
            for (int i = 0; i < neq; i++)
                out_m[i] = BI_m[i] * in_m[i];
        }
    }
}

static void chebyshev_smoother(struct All_variables *E, double **d0, double **F, double **Ad,
                               double acc, int *cycles, int level, int guess)
{
    (void)acc;
    const int neq = E->lmesh.NEQ[level];
    const int steps = *cycles;

    if (steps <= 0) {
        if (guess) {
            n_assemble_del2_u(E, d0, Ad, level, 1);
        } else {
            for (int m = 1; m <= E->sphere.caps_per_proc; m++)
                for (int i = 0; i < neq; i++)
                    d0[m][i] = Ad[m][i] = 0.0;
        }
        *cycles = 0;
        return;
    }

    static std::vector<double> ws_r[NCS], ws_z[NCS], ws_dx[NCS], ws_v[NCS], ws_Av[NCS];
    double *r[NCS], *z[NCS], *dx[NCS], *v[NCS], *Av[NCS];
    const size_t alloc_neq = (size_t)(neq + 1);
    for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
        r[m]  = ensure_mg_buf(ws_r[m],  alloc_neq);
        z[m]  = ensure_mg_buf(ws_z[m],  alloc_neq);
        dx[m] = ensure_mg_buf(ws_dx[m], alloc_neq);
        v[m]  = ensure_mg_buf(ws_v[m],  alloc_neq);
        Av[m] = ensure_mg_buf(ws_Av[m], alloc_neq);
    }

    /* Estimate and cache lambda_max(M^{-1} K) once per stiffness rebuild.
     * Initialize v from physical coordinates (x,y,z) so shared MPI boundary
     * nodes are 100% consistent across all ranks with zero extra communication. */
    double lmax = get_chebyshev_lmax(level);
    if (lmax <= 0.0) {
        const int nno_lev = E->lmesh.NNO[level];
        const double kfreq = (double)(E->mesh.ELZ[level] + 1) * 3.141592653589793;
        for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
            const struct ID *ID_m = E->ID[level][m];
            const double *th_arr = E->SX[level][m][1];
            const double *ph_arr = E->SX[level][m][2];
            const double *r_arr  = E->SX[level][m][3];
            for (int i = 1; i <= nno_lev; i++) {
                const double th = th_arr[i];
                const double ph = ph_arr[i];
                const double rr = r_arr[i];
                const double st = sin(th);
                const double cx = rr * st * cos(ph);
                const double cy = rr * st * sin(ph);
                const double cz = rr * cos(th);
                v[m][ID_m[i].doff[1]] = sin(kfreq * (1.13 * cx + 1.71 * cy + 2.37 * cz));
                v[m][ID_m[i].doff[2]] = cos(kfreq * (1.91 * cx - 1.37 * cy + 1.53 * cz));
                v[m][ID_m[i].doff[3]] = sin(kfreq * (2.23 * cx + 1.19 * cy - 1.83 * cz));
            }
        }
        strip_bcs_from_residual(E, v, level);
        double vn = sqrt(global_vdot(E, v, v, level));
        if (vn > 0.0) {
            const double inv_vn = 1.0 / vn;
            for (int m = 1; m <= E->sphere.caps_per_proc; m++)
                for (int i = 0; i < neq; i++)
                    v[m][i] *= inv_vn;

            double est = 2.0;
            const int n_power = 8;
            for (int it = 0; it < n_power; it++) {
                n_assemble_del2_u(E, v, Av, level, 1);
                apply_block3x3_inv(E, level, Av, z);
                double zz = 0.0, vz = 0.0;
                global_vdot2(E, z, z, v, level, &zz, &vz);
                const double zn = sqrt(zz);
                if (zn <= 1.0e-30) break;
                est = zn;
                const double inv_zn = 1.0 / zn;
                for (int m = 1; m <= E->sphere.caps_per_proc; m++)
                    for (int i = 0; i < neq; i++)
                        v[m][i] = z[m][i] * inv_zn;
            }
            lmax = est;
        } else {
            lmax = 2.0;
        }
        if (lmax < 1.0e-6) lmax = 2.0;
        set_chebyshev_lmax(level, lmax);
    }

    const double eig_max = 1.1 * lmax;
    const double eig_min = 0.12 * lmax;
    const double theta = 0.5 * (eig_max + eig_min);
    const double delta = 0.5 * (eig_max - eig_min);
    const double sigma1 = theta / delta;
    double rho = 1.0 / sigma1;

    if (guess) {
        n_assemble_del2_u(E, d0, Ad, level, 1);
        for (int m = 1; m <= E->sphere.caps_per_proc; m++)
            for (int i = 0; i < neq; i++)
                r[m][i] = F[m][i] - Ad[m][i];
    } else {
        for (int m = 1; m <= E->sphere.caps_per_proc; m++)
            for (int i = 0; i < neq; i++) {
                d0[m][i] = 0.0;
                r[m][i] = F[m][i];
            }
    }

    /* k = 0 */
    apply_block3x3_inv(E, level, r, z);
    const double inv_theta = 1.0 / theta;
    for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
        for (int i = 0; i < neq; i++) {
            const double d = inv_theta * z[m][i];
            dx[m][i] = d;
            d0[m][i] += d;
        }
    }

    /* k = 1 .. steps - 1 */
    for (int k = 1; k < steps; k++) {
        n_assemble_del2_u(E, d0, Ad, level, 1);
        for (int m = 1; m <= E->sphere.caps_per_proc; m++)
            for (int i = 0; i < neq; i++)
                r[m][i] = F[m][i] - Ad[m][i];

        apply_block3x3_inv(E, level, r, z);
        const double rho_next = 1.0 / (2.0 * sigma1 - rho);
        const double c_dx = rho * rho_next;
        const double c_z = 2.0 * rho_next / delta;
        for (int m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (int i = 0; i < neq; i++) {
                const double d = c_dx * dx[m][i] + c_z * z[m][i];
                dx[m][i] = d;
                d0[m][i] += d;
            }
        }
        rho = rho_next;
    }

    /* Maintain exact postcondition Ad = K * d0 */
    n_assemble_del2_u(E, d0, Ad, level, 1);
    *cycles = steps;
}
#endif /* !USE_CUDA */

/* Fast (conditional) determinant for 3x3 or 2x2 ... otherwise calls general routine */

double determinant(double A[4][4], int n)
{switch (n)
    { case 1:
	return(A[1][1]);
      case 2:
	return(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
      case 3:
	return(A[1][1]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])-
	       A[1][2]*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+
	       A[1][3]*(A[2][1]*A[3][2]-A[2][2]*A[3][1]));
      default:
	return(1);
/*	return(gen_determinant(A,n)); */
      }
}



double cofactor(double A[4][4], int i, int j, int n)
{
    int k,l,p,q;
  /* modern c++: double determinant() prototype in header */

  double B[4][4]; /* because of recursive behaviour of det/cofac, need to use
			       new copy of B at each 'n' level of this routine */

  if (n>3) printf("Error, no cofactors for matrix more than 3x3\n");

  p=q=1;

  for(k=1;k<=n;k++)    {
     if(k==i) continue;
     for(l=1;l<=n;l++)      {
	   if (l==j) continue;
           B[p][q]=A[k][l];
	   q++ ;
	   }
     q=1;p++;
     }


  return(epsilon[i][j]*determinant(B,n-1));


}

long double lg_pow(long double a, int n)
{
    /* compute the value of "a" raised to the power of "n" */
    long double b = 1.0;
    int i;

    for(i=0; i<n; i++) {
        b = b*a;
    }

    return(b);
}


