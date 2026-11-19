/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <mpi.h>
#include "element_definitions.h"
#include "global_defs.h"

namespace {
inline double* ensure_buf(std::vector<double>& vec, size_t n) {
    if (vec.size() < n) {
        vec.resize(n);
    }
    return vec.data();
}
} // namespace

void myerror(struct All_variables *E, const char *message);

static void solve_Ahat_p_fhat(struct All_variables *E,
                              double **V, double **P, double **F,
                              double imp, int *steps_max);
static void solve_Ahat_p_fhat_CG(struct All_variables *E,
                                 double **V, double **P, double **F,
                                 double imp, int *steps_max);
static void solve_Ahat_p_fhat_FCG(struct All_variables *E,
                                  double **V, double **P, double **F,
                                  double imp, int *steps_max);
static void solve_Ahat_p_fhat_FGMRES(struct All_variables *E,
                                     double **V, double **P, double **F,
                                     double imp, int *steps_max);
static void solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                    double **V, double **P, double **F,
                                    double imp, int *steps_max);
static void solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                      double **V, double **P, double **F,
                                      double imp, int *steps_max);
static void initial_vel_residual(struct All_variables *E,
                                 double **V, double **P, double **F,
                                 double imp);


/* Master loop for pressure and (hence) velocity field */

void solve_constrained_flow_iterative(struct All_variables *E)
{
    /* modern c++: void v_from_vector() prototype in header */
    /* modern c++: void v_from_vector_pseudo_surf() prototype in header */
    /* modern c++: void p_to_nodes() prototype in header */

    int cycles;

    cycles=E->control.p_iterations;

    /* Solve for velocity and pressure, correct for bc's */

    solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    if(E->control.pseudo_free_surf)
        v_from_vector_pseudo_surf(E);
    else
        v_from_vector(E);

    p_to_nodes(E,E->P,E->NP,E->mesh.levmax);

    return;
}


/* ========================================================================= */

static double momentum_eqn_residual(struct All_variables *E,
                                    double **V, double **P, double **F)
{
    /* Compute the norm of (F - grad(P) - K*V)
     * This norm is ~= E->monitor.momentum_residual */
    /* modern c++: void assemble_del2_u() prototype in header */
    /* modern c++: void assemble_grad_p() prototype in header */
    /* modern c++: void strip_bcs_from_residual() prototype in header */
    /* modern c++: double global_v_norm2() prototype in header */

    int i, m;
    double *r1[NCS], *r2[NCS];
    double res;
    const int lev = E->mesh.levmax;
    const int neq = E->lmesh.neq;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        r1[m] = (double *)malloc((neq+1)*sizeof(double));
        r2[m] = (double *)malloc((neq+1)*sizeof(double));
    }

    /* r2 = F - grad(P) - K*V */
    assemble_grad_p(E, P, E->u1, lev);
    assemble_del2_u(E, V, r1, lev, 1);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            r2[m][i] = F[m][i] - E->u1[m][i] - r1[m][i];

    strip_bcs_from_residual(E, r2, lev);

    res = sqrt(global_v_norm2(E, r2));

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        free(r1[m]);
        free(r2[m]);
    }
    return(res);
}


static void print_convergence_progress(struct All_variables *E,
                                       int count, double time0,
                                       double v_norm, double p_norm,
                                       double dv, double dp,
                                       double div)
{
    double t; /* funcs moved to header: CPU_time0() */
    t = CPU_time0() - time0;

    fprintf(E->fp, "(%03d) %5.1f s v=%e p=%e "
            "div/v=%.2e dv/v=%.2e dp/p=%.2e step %d\n",
            count, t, v_norm, p_norm, div, dv, dp,
            E->monitor.solution_cycles);
    fprintf(stderr, "(%03d) %5.1f s v=%e p=%e "
            "div/v=%.2e dv/v=%.2e dp/p=%.2e step %d\n",
            count, t, v_norm, p_norm, div, dv, dp,
            E->monitor.solution_cycles);
    fflush(stderr);
    return;
}


static int keep_iterating(struct All_variables *E,
                          double acc, int converging)
{
    const int required_converging_loops = 2;

    if(E->control.check_continuity_convergence)
        return (E->monitor.incompressibility > acc) ||
	    (converging < required_converging_loops);
    else
        return (E->monitor.incompressibility > acc) &&
	    (converging < required_converging_loops);
}


static void solve_Ahat_p_fhat(struct All_variables *E,
                               double **V, double **P, double **F,
                               double imp, int *steps_max)
{
    if(strcmp(E->control.uzawa, "fgmres") == 0) {
        solve_Ahat_p_fhat_FGMRES(E, V, P, F, imp, steps_max);
    }
    else if(E->control.inv_gruneisen == 0) {
        if(strcmp(E->control.uzawa, "fcg") == 0)
            solve_Ahat_p_fhat_FCG(E, V, P, F, imp, steps_max);
        else
            solve_Ahat_p_fhat_CG(E, V, P, F, imp, steps_max);
    }
    else {
        if(strcmp(E->control.uzawa, "cg") == 0 || strcmp(E->control.uzawa, "fcg") == 0)
            solve_Ahat_p_fhat_iterCG(E, V, P, F, imp, steps_max);
        else if(strcmp(E->control.uzawa, "bicg") == 0)
            solve_Ahat_p_fhat_BiCG(E, V, P, F, imp, steps_max);
        else
            myerror(E, "Error: unknown Uzawa iteration\n");
    }

    return;
}


/* Solve incompressible Stokes flow using
 * conjugate gradient (CG) iterations
 */

static void solve_Ahat_p_fhat_CG(struct All_variables *E,
                                 double **V, double **P, double **FF,
                                 double imp, int *steps_max)
{
    int m, j, count, valid, lev, npno, neq;

    double *r1[NCS], *r2[NCS], *z1[NCS], *s1[NCS], *s2[NCS], *cu[NCS];
    double *F[NCS];
    double *shuffle[NCS];
    double alpha, delta, r0dotz0, r1dotz1;
    double v_res;
    double inner_imp;
    /* modern c++: double global_pdot() prototype in header */
    /* modern c++: global_v_norm2(), global_p_norm2(), global_div_norm2() prototypes in header */

    double time0; /* funcs moved to header: CPU_time0() */
    double v_norm, p_norm;
    double dvelocity, dpressure;
    int converging;
    /* modern c++: void assemble_c_u() prototype in header */
    /* modern c++: void assemble_div_u() prototype in header */
    /* modern c++: void assemble_del2_u() prototype in header */
    /* modern c++: void assemble_grad_p() prototype in header */
    /* modern c++: void strip_bcs_from_residual() prototype in header */
    /* modern c++: int solve_del2_u() prototype in header */
    /* modern c++: void parallel_process_termination() prototype in header */
    /* modern c++: void v_from_vector() prototype in header */
    /* modern c++: void v_from_vector_pseudo_surf() prototype in header */
    /* modern c++: void assign_v_to_vector() prototype in header */

    inner_imp = imp * E->control.inner_accuracy_scale; /* allow for different innner loop accuracy */

    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    static std::vector<double> ws_F[NCS], ws_r1[NCS], ws_r2[NCS], ws_z1[NCS], ws_s1[NCS], ws_s2[NCS], ws_cu[NCS];
    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[m]  = ensure_buf(ws_F[m],  (size_t)neq);
        r1[m] = ensure_buf(ws_r1[m], (size_t)(npno+1));
        r2[m] = ensure_buf(ws_r2[m], (size_t)(npno+1));
        z1[m] = ensure_buf(ws_z1[m], (size_t)(npno+1));
        s1[m] = ensure_buf(ws_s1[m], (size_t)(npno+1));
        s2[m] = ensure_buf(ws_s2[m], (size_t)(npno+1));
        cu[m] = ensure_buf(ws_cu[m], (size_t)(npno+1));
    }

    time0 = CPU_time0();
    count = 0;
    v_res = E->monitor.fdotf;

    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[m][j] = FF[m][j];


    /* calculate the contribution of compressibility in the continuity eqn */
    if(E->control.inv_gruneisen != 0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++)
                cu[m][j] = 0.0;

        assemble_c_u(E, V, cu, lev);
    }


    /* calculate the initial velocity residual */
    /* In the compressible case, the initial guess of P might be bad.
     * Do not correct V with it. */
    if(E->control.inv_gruneisen == 0)
        initial_vel_residual(E, V, P, F, inner_imp*v_res);


    /* initial residual r1 = div(V) */
    assemble_div_u(E, V, r1, lev);


    /* add the contribution of compressibility to the initial residual */
    if(E->control.inv_gruneisen != 0)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++) {
                r1[m][j] += cu[m][j];
            }

    E->monitor.vdotv = global_v_norm2(E, V);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, r1)
                                        / (1e-32 + E->monitor.vdotv));

    v_norm = sqrt(E->monitor.vdotv);
    p_norm = sqrt(E->monitor.pdotp);
    dvelocity = 1.0;
    dpressure = 1.0;
    converging = 0;

    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0,
                                   v_norm, p_norm,
                                   dvelocity, dpressure,
                                   E->monitor.incompressibility);
    }

  
    r0dotz0 = 0;

    while( (count < *steps_max) && keep_iterating(E, imp, converging) ) {
        /* require two consecutive converging iterations to quit the while-loop */

        /* preconditioner BPI ~= inv(K), z1 = BPI*r1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                z1[m][j] = E->BPI[lev][m][j] * r1[m][j];


        /* r1dotz1 = <r1, z1> */
        r1dotz1 = global_pdot(E, r1, z1, lev);
        assert(r1dotz1 != 0.0  /* Division by zero in head of incompressibility iteration */);

        /* update search direction */
        if(count == 0)
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    s2[m][j] = z1[m][j];
        else {
            /* s2 = z1 + s1 * <r1,z1>/<r0,z0> */
            delta = r1dotz1 / r0dotz0;
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    s2[m][j] = z1[m][j] + delta * s1[m][j];
        }

        /* solve K*u1 = grad(s2) for u1 */
        assemble_grad_p(E, s2, F, lev);
        valid = solve_del2_u(E, E->u1, F, inner_imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
        }
        strip_bcs_from_residual(E, E->u1, lev);


        /* F = div(u1) */
        assemble_div_u(E, E->u1, F, lev);


        /* alpha = <r1, z1> / <s2, F> */
        alpha = r1dotz1 / global_pdot(E, s2, F, lev);


        /* r2 = r1 - alpha * div(u1) */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                r2[m][j] = r1[m][j] - alpha * F[m][j];


        /* P = P + alpha * s2 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                P[m][j] += alpha * s2[m][j];


        /* V = V - alpha * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[m][j] -= alpha * E->u1[m][j];


        /* compute velocity, pressure, delta, and incompressibility residuals in 1 fused MPI_Allreduce */
        assemble_div_u(E, V, z1, lev);
        if(E->control.inv_gruneisen != 0)
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(j=1;j<=npno;j++) {
                    z1[m][j] += cu[m][j];
            }

        {
            double norms[5];
            global_stokes_norms(E, V, P, E->u1, s2, z1, norms);
            E->monitor.vdotv = norms[0];
            E->monitor.pdotp = norms[1];
            v_norm = sqrt(norms[0]);
            p_norm = sqrt(norms[1]);
            dvelocity = alpha * sqrt(norms[2] / (1e-32 + norms[0]));
            dpressure = alpha * sqrt(norms[3] / (1e-32 + norms[1]));
            E->monitor.incompressibility = sqrt(norms[4] / (1e-32 + norms[0]));
        }

        count++;


        if (E->control.print_convergence && E->parallel.me==0)  {
            print_convergence_progress(E, count, time0,
                                       v_norm, p_norm,
                                       dvelocity, dpressure,
                                       E->monitor.incompressibility);
        }

	if(!valid){
            /* reset consecutive converging iterations */
            converging = 0;
	}else{
            /* how many consecutive converging iterations? */
            if(E->control.check_pressure_convergence) {
                /* check dv and dp */
                if(dvelocity < imp && dpressure < imp)
                    converging++;
                else
                    converging = 0;
            }else{
                /* check dv only */
                if(dvelocity < imp)
                    converging++;
                else
                    converging = 0;
            }
	  
	}

        /* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[m] = s1[m];
            s1[m] = s2[m];
            s2[m] = shuffle[m];

            shuffle[m] = r1[m];
            r1[m] = r2[m];
            r2[m] = shuffle[m];
        }

        /* shift <r0, z0> = <r1, z1> */
        r0dotz0 = r1dotz1;
	if((E->sphere.caps == 12) && (E->control.inner_remove_rigid_rotation)){
	  /* allow for removal of net rotation at each iterative step
	     (expensive) */
	  if(E->control.pseudo_free_surf) /* move from U to V */
	    v_from_vector_pseudo_surf(E);
	  else
	    v_from_vector(E);
	  remove_rigid_rot(E);	/* correct V */
	  assign_v_to_vector(E); /* assign V to U */
	}
    } /* end loop for conjugate gradient */

    assemble_div_u(E, V, z1, lev);
    if(E->control.inv_gruneisen != 0)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++) {
                z1[m][j] += cu[m][j];
            }


    /* Persistent workspace retained across timesteps */
    *steps_max=count;

    return;
}


/* =========================================================================
 * Solve Stokes Schur complement using Truncated Flexible Conjugate Gradient
 * (Notay's FCG(m_max)) with fused multi-dot MPI_Allreduce and periodic
 * true divergence residual refresh.
 * Tolerates variable/inexact/asymmetric inner Multigrid solves without
 * losing S-conjugacy or stalling/rebounding div/v.
 * ========================================================================= */
static void solve_Ahat_p_fhat_FCG(struct All_variables *E,
                                  double **V, double **P, double **FF,
                                  double imp, int *steps_max)
{
    int m, j, count, valid, lev, npno, neq;
    double *r1[NCS], *z1[NCS], *s2[NCS], *w_cur[NCS], *cu[NCS], *F[NCS];
    double alpha, v_res, inner_imp, time0;
    double v_norm, p_norm, dvelocity, dpressure;
    int converging;

    inner_imp = imp * E->control.inner_accuracy_scale;
    npno = E->lmesh.npno;
    neq  = E->lmesh.neq;
    lev  = E->mesh.levmax;

    const int caps = E->sphere.caps_per_proc;
    const int mmax = (E->control.fcg_mmax > 0) ? E->control.fcg_mmax : 10;
    const int refresh_every = (E->control.residual_refresh_every > 0) ? E->control.residual_refresh_every : 20;

    static std::vector<double> ws_F[NCS], ws_r1[NCS], ws_z1[NCS], ws_s2[NCS], ws_w[NCS], ws_cu[NCS];
    static std::vector<double> ws_hist_s[NCS], ws_hist_w[NCS];
    static std::vector<double> hist_d, local_g, glob_g;

    for (m = 1; m <= caps; m++) {
        F[m]     = ensure_buf(ws_F[m],  (size_t)neq);
        r1[m]    = ensure_buf(ws_r1[m], (size_t)(npno + 1));
        z1[m]    = ensure_buf(ws_z1[m], (size_t)(npno + 1));
        s2[m]    = ensure_buf(ws_s2[m], (size_t)(npno + 1));
        w_cur[m] = ensure_buf(ws_w[m],  (size_t)(npno + 1));
        cu[m]    = ensure_buf(ws_cu[m], (size_t)(npno + 1));
        ensure_buf(ws_hist_s[m], (size_t)mmax * (size_t)(npno + 1));
        ensure_buf(ws_hist_w[m], (size_t)mmax * (size_t)(npno + 1));
    }
    hist_d.resize((size_t)mmax, 0.0);
    local_g.resize((size_t)mmax, 0.0);
    glob_g.resize((size_t)mmax, 0.0);

    time0 = CPU_time0();
    count = 0;
    v_res = E->monitor.fdotf;

    for (m = 1; m <= caps; m++)
        memcpy(F[m], FF[m], (size_t)neq * sizeof(double));

    if (E->control.inv_gruneisen != 0) {
        for (m = 1; m <= caps; m++)
            for (j = 1; j <= npno; j++)
                cu[m][j] = 0.0;
        assemble_c_u(E, V, cu, lev);
    }

    if (E->control.inv_gruneisen == 0)
        initial_vel_residual(E, V, P, F, inner_imp * v_res);

    assemble_div_u(E, V, r1, lev);
    if (E->control.inv_gruneisen != 0)
        for (m = 1; m <= caps; m++)
            for (j = 1; j <= npno; j++)
                r1[m][j] += cu[m][j];

    E->monitor.vdotv = global_v_norm2(E, V);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, r1) / (1e-32 + E->monitor.vdotv));

    v_norm = sqrt(E->monitor.vdotv);
    p_norm = sqrt(E->monitor.pdotp);
    dvelocity = 1.0;
    dpressure = 1.0;
    converging = 0;

    if (E->control.print_convergence && E->parallel.me == 0) {
        print_convergence_progress(E, count, time0, v_norm, p_norm,
                                   dvelocity, dpressure, E->monitor.incompressibility);
    }

    int hist_len = 0;
    int hist_head = 0;

    while ((count < *steps_max) && keep_iterating(E, imp, converging)) {
        /* 1. Precondition residual: z1 = BPI * r1 */
        for (m = 1; m <= caps; m++) {
            const double * __restrict__ bpi = E->BPI[lev][m];
            const double * __restrict__ rm  = r1[m];
            double * __restrict__ zm        = z1[m];
            for (j = 1; j <= npno; j++)
                zm[j] = bpi[j] * rm[j];
        }

        /* 2. Truncated S-orthogonalization against last hist_len search directions */
        if (hist_len > 0) {
            for (int h = 0; h < hist_len; h++)
                local_g[h] = 0.0;
            for (m = 1; m <= caps; m++) {
                const double * __restrict__ zm = z1[m];
                const double *hw_base = ws_hist_w[m].data();
                for (int h = 0; h < hist_len; h++) {
                    const double * __restrict__ hw = hw_base + (size_t)h * (size_t)(npno + 1);
                    double acc = 0.0;
                    for (j = 1; j <= npno; j++)
                        acc += zm[j] * hw[j];
                    local_g[h] += acc;
                }
            }
            MPI_Allreduce(local_g.data(), glob_g.data(), hist_len, MPI_DOUBLE, MPI_SUM, E->parallel.world);

            for (m = 1; m <= caps; m++) {
                double * __restrict__ sm = s2[m];
                const double * __restrict__ zm = z1[m];
                for (j = 1; j <= npno; j++)
                    sm[j] = zm[j];
                const double *hs_base = ws_hist_s[m].data();
                for (int h = 0; h < hist_len; h++) {
                    const double coeff = glob_g[h] / hist_d[h];
                    const double * __restrict__ hs = hs_base + (size_t)h * (size_t)(npno + 1);
                    for (j = 1; j <= npno; j++)
                        sm[j] -= coeff * hs[j];
                }
            }
        } else {
            for (m = 1; m <= caps; m++)
                memcpy(&s2[m][1], &z1[m][1], (size_t)npno * sizeof(double));
        }

        /* 3. Inner velocity solve: K * u1 = grad(s2) */
        assemble_grad_p(E, s2, F, lev);
        valid = solve_del2_u(E, E->u1, F, inner_imp * v_res, lev);
        if (!valid && (E->parallel.me == 0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
        }
        strip_bcs_from_residual(E, E->u1, lev);

        /* 4. Schur image: w_cur = div(u1) */
        assemble_div_u(E, E->u1, w_cur, lev);

        /* 5. Fused 2-double Allreduce for d_k = <s2, w_cur> and num_k = <s2, r1> */
        double loc2[2] = {0.0, 0.0}, glob2[2] = {0.0, 0.0};
        for (m = 1; m <= caps; m++) {
            const double * __restrict__ sm = s2[m];
            const double * __restrict__ wm = w_cur[m];
            const double * __restrict__ rm = r1[m];
            double sd = 0.0, sn = 0.0;
            for (j = 1; j <= npno; j++) {
                sd += sm[j] * wm[j];
                sn += sm[j] * rm[j];
            }
            loc2[0] += sd;
            loc2[1] += sn;
        }
        MPI_Allreduce(loc2, glob2, 2, MPI_DOUBLE, MPI_SUM, E->parallel.world);
        double d_k = glob2[0];
        double num_k = glob2[1];

        /* Guard against loss of positive curvature from extreme inexactness */
        if (d_k <= 1e-30) {
            hist_len = 0;
            hist_head = 0;
            /* Refresh true residual and retry on next iteration */
            assemble_div_u(E, V, r1, lev);
            if (E->control.inv_gruneisen != 0)
                for (m = 1; m <= caps; m++)
                    for (j = 1; j <= npno; j++)
                        r1[m][j] += cu[m][j];
            count++;
            continue;
        }

        alpha = num_k / d_k;

        /* 6. Update P, V, and recursive residual r1 */
        for (m = 1; m <= caps; m++) {
            double * __restrict__ Pm = P[m];
            double * __restrict__ rm = r1[m];
            const double * __restrict__ sm = s2[m];
            const double * __restrict__ wm = w_cur[m];
            for (j = 1; j <= npno; j++) {
                Pm[j] += alpha * sm[j];
                rm[j] -= alpha * wm[j];
            }
            double * __restrict__ Vm = V[m];
            const double * __restrict__ um = E->u1[m];
            for (j = 0; j < neq; j++)
                Vm[j] -= alpha * um[j];
        }

        /* 7. Store (s2, w_cur, d_k) in circular FCG history */
        for (m = 1; m <= caps; m++) {
            double *hs = ws_hist_s[m].data() + (size_t)hist_head * (size_t)(npno + 1);
            double *hw = ws_hist_w[m].data() + (size_t)hist_head * (size_t)(npno + 1);
            memcpy(&hs[1], &s2[m][1], (size_t)npno * sizeof(double));
            memcpy(&hw[1], &w_cur[m][1], (size_t)npno * sizeof(double));
        }
        hist_d[hist_head] = d_k;
        hist_head = (hist_head + 1) % mmax;
        if (hist_len < mmax)
            hist_len++;

        /* 8. Periodic true residual refresh (eliminates recursive residual drift) */
        const int do_refresh = (((count + 1) % refresh_every) == 0);
        if (do_refresh) {
            assemble_div_u(E, V, r1, lev);
            if (E->control.inv_gruneisen != 0)
                for (m = 1; m <= caps; m++)
                    for (j = 1; j <= npno; j++)
                        r1[m][j] += cu[m][j];
        }

        /* 9. Evaluate all 5 norms in 1 fused MPI_Allreduce */
        {
            double norms[5];
            global_stokes_norms(E, V, P, E->u1, s2, r1, norms);
            E->monitor.vdotv = norms[0];
            E->monitor.pdotp = norms[1];
            v_norm = sqrt(norms[0]);
            p_norm = sqrt(norms[1]);
            dvelocity = fabs(alpha) * sqrt(norms[2] / (1e-32 + norms[0]));
            dpressure = fabs(alpha) * sqrt(norms[3] / (1e-32 + norms[1]));
            E->monitor.incompressibility = sqrt(norms[4] / (1e-32 + norms[0]));
        }

        count++;

        if (!valid) {
            converging = 0;
        } else {
            if (E->control.check_pressure_convergence) {
                if (dvelocity < imp && dpressure < imp)
                    converging++;
                else
                    converging = 0;
            } else {
                if (dvelocity < imp)
                    converging++;
                else
                    converging = 0;
            }
        }

        /* Before exiting on convergence, verify true physical divergence residual once */
        if (!do_refresh && !keep_iterating(E, imp, converging)) {
            assemble_div_u(E, V, r1, lev);
            if (E->control.inv_gruneisen != 0)
                for (m = 1; m <= caps; m++)
                    for (j = 1; j <= npno; j++)
                        r1[m][j] += cu[m][j];
            E->monitor.incompressibility = sqrt(global_div_norm2(E, r1) / (1e-32 + E->monitor.vdotv));
        }

        if (E->control.print_convergence && E->parallel.me == 0) {
            print_convergence_progress(E, count, time0, v_norm, p_norm,
                                       dvelocity, dpressure, E->monitor.incompressibility);
        }

        if ((E->sphere.caps == 12) && (E->control.inner_remove_rigid_rotation)) {
            if (E->control.pseudo_free_surf)
                v_from_vector_pseudo_surf(E);
            else
                v_from_vector(E);
            remove_rigid_rot(E);
            assign_v_to_vector(E);
        }
    }

    *steps_max = count;
    return;
}


/* =========================================================================
 * Solve Stokes Schur complement using Right-Preconditioned Restarted
 * Flexible GMRES (FGMRES(m)) in CitcomS's element-volume-weighted W-norm
 * <a, b>_W = sum_e (a_e * b_e / V_e).
 * - Guarantees strictly monotonic decrease of E->monitor.incompressibility!
 * - Works identically for BOTH incompressible (BA/EBA) and compressible (ALA/TALA).
 * - Uses CGS2 (2-pass Classical Gram-Schmidt) so each Krylov step needs only
 *   2 fused multi-dot MPI_Allreduce calls regardless of Krylov dimension m.
 * ========================================================================= */
static void solve_Ahat_p_fhat_FGMRES(struct All_variables *E,
                                     double **V, double **P, double **FF,
                                     double imp, int *steps_max)
{
    int m, j, count, valid, lev, npno, neq;
    double *F[NCS], *r0[NCS], *w_vec[NCS], *dP_step[NCS], *dV_step[NCS];
    double v_res, inner_imp, time0;
    double v_norm, p_norm, dvelocity, dpressure;
    int converging;

    inner_imp = imp * E->control.inner_accuracy_scale;
    npno = E->lmesh.npno;
    neq  = E->lmesh.neq;
    lev  = E->mesh.levmax;

    const int caps = E->sphere.caps_per_proc;
    const int m_res = (E->control.fgmres_restart >= 2) ? E->control.fgmres_restart : 30;
    const int compressible = (E->control.inv_gruneisen != 0);
    const double inv_vol = 1.0 / E->mesh.volume;

    static std::vector<double> ws_F[NCS], ws_r0[NCS], ws_w[NCS], ws_dP[NCS], ws_dV[NCS];
    static std::vector<double> ws_Vb[NCS], ws_Zb[NCS], ws_Ub[NCS];
    static std::vector<double> H, cs, sn, g, y, y_prev, loc_h, glob_h;

    for (m = 1; m <= caps; m++) {
        F[m]       = ensure_buf(ws_F[m],  (size_t)neq);
        r0[m]      = ensure_buf(ws_r0[m], (size_t)(npno + 1));
        w_vec[m]   = ensure_buf(ws_w[m],  (size_t)(npno + 1));
        dP_step[m] = ensure_buf(ws_dP[m], (size_t)(npno + 1));
        dV_step[m] = ensure_buf(ws_dV[m], (size_t)neq);
        ensure_buf(ws_Vb[m], (size_t)(m_res + 1) * (size_t)(npno + 1));
        ensure_buf(ws_Zb[m], (size_t)m_res * (size_t)(npno + 1));
        ensure_buf(ws_Ub[m], (size_t)m_res * (size_t)neq);
    }

    const int ldH = m_res + 1;
    H.assign((size_t)ldH * (size_t)m_res, 0.0);
    cs.assign((size_t)m_res, 0.0);
    sn.assign((size_t)m_res, 0.0);
    g.assign((size_t)(m_res + 1), 0.0);
    y.assign((size_t)m_res, 0.0);
    y_prev.assign((size_t)m_res, 0.0);
    loc_h.resize((size_t)(m_res + 1), 0.0);
    glob_h.resize((size_t)(m_res + 1), 0.0);

    time0 = CPU_time0();
    count = 0;
    v_res = E->monitor.fdotf;

    for (m = 1; m <= caps; m++)
        memcpy(F[m], FF[m], (size_t)neq * sizeof(double));

    initial_vel_residual(E, V, P, F, inner_imp * v_res);

    if (compressible)
        assemble_div_rho_u(E, V, r0, lev);
    else
        assemble_div_u(E, V, r0, lev);

    E->monitor.vdotv = global_v_norm2(E, V);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, r0) / (1e-32 + E->monitor.vdotv));

    v_norm = sqrt(E->monitor.vdotv);
    p_norm = sqrt(E->monitor.pdotp);
    dvelocity = 1.0;
    dpressure = 1.0;
    converging = 0;

    if (E->control.print_convergence && E->parallel.me == 0) {
        print_convergence_progress(E, count, time0, v_norm, p_norm,
                                   dvelocity, dpressure, E->monitor.incompressibility);
    }

    /* Outer restart loop */
    while ((count < *steps_max) && keep_iterating(E, imp, converging)) {
        if (count > 0) {
            /* Recompute exact true physical residual at start of every restart cycle */
            if (compressible)
                assemble_div_rho_u(E, V, r0, lev);
            else
                assemble_div_u(E, V, r0, lev);
        }

        /* beta = ||r0||_W */
        double loc_beta2 = 0.0, glob_beta2 = 0.0;
        for (m = 1; m <= caps; m++) {
            const double * __restrict__ rm = r0[m];
            const struct COORD * __restrict__ ecom = E->eco[m];
            double s = 0.0;
            for (j = 1; j <= npno; j++)
                s += (rm[j] * rm[j]) / ecom[j].area;
            loc_beta2 += s;
        }
        MPI_Allreduce(&loc_beta2, &glob_beta2, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
        double beta = sqrt(glob_beta2);
        if (beta <= 1e-32)
            break;

        /* v_0 = r0 / beta */
        const double inv_beta = 1.0 / beta;
        for (m = 1; m <= caps; m++) {
            double * __restrict__ v0 = ws_Vb[m].data();
            const double * __restrict__ rm = r0[m];
            for (j = 1; j <= npno; j++)
                v0[j] = rm[j] * inv_beta;
        }

        std::fill(g.begin(), g.end(), 0.0);
        std::fill(y_prev.begin(), y_prev.end(), 0.0);
        g[0] = beta;

        /* Inner Arnoldi cycle k = 0 .. m_res - 1 */
        for (int k = 0; k < m_res && count < *steps_max && keep_iterating(E, imp, converging); k++) {
            double *Zk_ptrs[NCS];
            double *Uk_ptrs[NCS];
            for (m = 1; m <= caps; m++) {
                const double * __restrict__ vk = ws_Vb[m].data() + (size_t)k * (size_t)(npno + 1);
                double * __restrict__ zk = ws_Zb[m].data() + (size_t)k * (size_t)(npno + 1);
                const double * __restrict__ bpi = E->BPI[lev][m];
                for (j = 1; j <= npno; j++)
                    zk[j] = bpi[j] * vk[j];
                Zk_ptrs[m] = zk;
                Uk_ptrs[m] = ws_Ub[m].data() + (size_t)k * (size_t)neq;
            }

            /* Flexible right preconditioner action: K * U_k = grad(Z_k) */
            assemble_grad_p(E, Zk_ptrs, F, lev);
            valid = solve_del2_u(E, Uk_ptrs, F, inner_imp * v_res, lev);
            if (!valid && (E->parallel.me == 0)) {
                fputs("Warning: solver not converging! 1\n", stderr);
                fputs("Warning: solver not converging! 1\n", E->fp);
            }
            strip_bcs_from_residual(E, Uk_ptrs, lev);

            /* Schur action: w = div(U_k) or div(rho * U_k) */
            if (compressible)
                assemble_div_rho_u(E, Uk_ptrs, w_vec, lev);
            else
                assemble_div_u(E, Uk_ptrs, w_vec, lev);

            /* CGS2 orthogonalization in W-norm (2 fused Allreduces, O(eps) orthogonality) */
            double *Hk = &H[(size_t)k * (size_t)ldH];
            for (int i = 0; i <= k; i++) Hk[i] = 0.0;

            for (int pass = 0; pass < 2; pass++) {
                for (int i = 0; i <= k; i++) loc_h[i] = 0.0;
                for (m = 1; m <= caps; m++) {
                    const double * __restrict__ wm = w_vec[m];
                    const struct COORD * __restrict__ ecom = E->eco[m];
                    const double *vb_base = ws_Vb[m].data();
                    for (int i = 0; i <= k; i++) {
                        const double * __restrict__ vi = vb_base + (size_t)i * (size_t)(npno + 1);
                        double acc = 0.0;
                        for (j = 1; j <= npno; j++)
                            acc += (wm[j] * vi[j]) / ecom[j].area;
                        loc_h[i] += acc;
                    }
                }
                MPI_Allreduce(loc_h.data(), glob_h.data(), k + 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
                for (int i = 0; i <= k; i++) Hk[i] += glob_h[i];

                for (m = 1; m <= caps; m++) {
                    double * __restrict__ wm = w_vec[m];
                    const double *vb_base = ws_Vb[m].data();
                    for (int i = 0; i <= k; i++) {
                        const double hij = glob_h[i];
                        const double * __restrict__ vi = vb_base + (size_t)i * (size_t)(npno + 1);
                        for (j = 1; j <= npno; j++)
                            wm[j] -= hij * vi[j];
                    }
                }
            }

            /* h_{k+1, k} = ||w||_W */
            double loc_w2 = 0.0, glob_w2 = 0.0;
            for (m = 1; m <= caps; m++) {
                const double * __restrict__ wm = w_vec[m];
                const struct COORD * __restrict__ ecom = E->eco[m];
                double s = 0.0;
                for (j = 1; j <= npno; j++)
                    s += (wm[j] * wm[j]) / ecom[j].area;
                loc_w2 += s;
            }
            MPI_Allreduce(&loc_w2, &glob_w2, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
            double h_next = sqrt(glob_w2);
            Hk[k + 1] = h_next;

            if (h_next > 1e-32) {
                const double inv_h = 1.0 / h_next;
                for (m = 1; m <= caps; m++) {
                    double * __restrict__ vnext = ws_Vb[m].data() + (size_t)(k + 1) * (size_t)(npno + 1);
                    const double * __restrict__ wm = w_vec[m];
                    for (j = 1; j <= npno; j++)
                        vnext[j] = wm[j] * inv_h;
                }
            }

            /* Apply previous Givens rotations to Hk */
            for (int i = 0; i < k; i++) {
                const double tmp = cs[i] * Hk[i] + sn[i] * Hk[i + 1];
                Hk[i + 1] = -sn[i] * Hk[i] + cs[i] * Hk[i + 1];
                Hk[i] = tmp;
            }

            /* Compute new Givens rotation for (Hk[k], Hk[k+1]) */
            const double denom = hypot(Hk[k], Hk[k + 1]);
            if (denom > 1e-32) {
                cs[k] = Hk[k] / denom;
                sn[k] = Hk[k + 1] / denom;
            } else {
                cs[k] = 1.0;
                sn[k] = 0.0;
            }
            Hk[k] = cs[k] * Hk[k] + sn[k] * Hk[k + 1];
            Hk[k + 1] = 0.0;

            const double gk_new = cs[k] * g[k];
            g[k + 1] = -sn[k] * g[k];
            g[k] = gk_new;

            /* Back-substitute (k+1)x(k+1) upper triangular H * y = g */
            for (int i = k; i >= 0; i--) {
                double rhs = g[i];
                for (int col = i + 1; col <= k; col++)
                    rhs -= H[(size_t)col * (size_t)ldH + (size_t)i] * y[col];
                const double diag = H[(size_t)i * (size_t)ldH + (size_t)i];
                y[i] = (fabs(diag) > 1e-32) ? (rhs / diag) : 0.0;
            }

            /* Incremental update: dP = sum_{i=0..k} (y[i] - y_prev[i]) * Z_i, dV = sum ... * U_i */
            for (m = 1; m <= caps; m++) {
                double * __restrict__ dPm = dP_step[m];
                double * __restrict__ dVm = dV_step[m];
                for (j = 1; j <= npno; j++) dPm[j] = 0.0;
                for (j = 0; j < neq; j++)   dVm[j] = 0.0;

                const double *zb_base = ws_Zb[m].data();
                const double *ub_base = ws_Ub[m].data();
                for (int i = 0; i <= k; i++) {
                    const double dy = y[i] - y_prev[i];
                    if (dy == 0.0) continue;
                    const double * __restrict__ zi = zb_base + (size_t)i * (size_t)(npno + 1);
                    const double * __restrict__ ui = ub_base + (size_t)i * (size_t)neq;
                    for (j = 1; j <= npno; j++) dPm[j] += dy * zi[j];
                    for (j = 0; j < neq; j++)   dVm[j] += dy * ui[j];
                }

                double * __restrict__ Pm = P[m];
                double * __restrict__ Vm = V[m];
                for (j = 1; j <= npno; j++) Pm[j] += dPm[j];
                for (j = 0; j < neq; j++)   Vm[j] -= dVm[j];
            }
            for (int i = 0; i <= k; i++) y_prev[i] = y[i];

            /* Fused 4-norm reduction (div norm comes directly from exact |g[k+1]|!) */
            {
                double loc4[4] = {0.0, 0.0, 0.0, 0.0};
                double glob4[4] = {0.0, 0.0, 0.0, 0.0};
                const int nno = E->lmesh.nno;
                for (m = 1; m <= caps; m++) {
                    const double * __restrict__ Vm = V[m];
                    const double * __restrict__ dVm = dV_step[m];
                    const double * __restrict__ nmass = E->NMass[m];
                    const struct ID * __restrict__ idm = E->id[m];
                    double sv = 0.0, sdv = 0.0;
                    for (j = 1; j <= nno; j++) {
                        const int e1 = idm[j].doff[1], e2 = idm[j].doff[2], e3 = idm[j].doff[3];
                        const double w = nmass[j];
                        sv  += (Vm[e1]*Vm[e1] + Vm[e2]*Vm[e2] + Vm[e3]*Vm[e3]) * w;
                        sdv += (dVm[e1]*dVm[e1] + dVm[e2]*dVm[e2] + dVm[e3]*dVm[e3]) * w;
                    }
                    loc4[0] += sv;
                    loc4[2] += sdv;

                    const double * __restrict__ Pm = P[m];
                    const double * __restrict__ dPm = dP_step[m];
                    const struct COORD * __restrict__ ecom = E->eco[m];
                    double sp = 0.0, sdp = 0.0;
                    for (j = 1; j <= npno; j++) {
                        const double area = ecom[j].area;
                        sp  += Pm[j] * Pm[j] * area;
                        sdp += dPm[j] * dPm[j] * area;
                    }
                    loc4[1] += sp;
                    loc4[3] += sdp;
                }
                MPI_Allreduce(loc4, glob4, 4, MPI_DOUBLE, MPI_SUM, E->parallel.world);
                E->monitor.vdotv = glob4[0] * inv_vol;
                E->monitor.pdotp = glob4[1] * inv_vol;
                v_norm = sqrt(E->monitor.vdotv);
                p_norm = sqrt(E->monitor.pdotp);
                dvelocity = sqrt((glob4[2] * inv_vol) / (1e-32 + E->monitor.vdotv));
                dpressure = sqrt((glob4[3] * inv_vol) / (1e-32 + E->monitor.pdotp));
                const double div_norm2 = (g[k + 1] * g[k + 1]) * inv_vol;
                E->monitor.incompressibility = sqrt(div_norm2 / (1e-32 + E->monitor.vdotv));
            }

            count++;

            if (!valid) {
                converging = 0;
            } else {
                if (E->control.check_pressure_convergence) {
                    if (dvelocity < imp && dpressure < imp)
                        converging++;
                    else
                        converging = 0;
                } else {
                    if (dvelocity < imp)
                        converging++;
                    else
                        converging = 0;
                }
            }

            if (E->control.print_convergence && E->parallel.me == 0) {
                print_convergence_progress(E, count, time0, v_norm, p_norm,
                                           dvelocity, dpressure, E->monitor.incompressibility);
            }

            if (h_next <= 1e-32)
                break; /* Happy breakdown: exact solution reached in subspace */
        } /* end inner Arnoldi cycle */
    } /* end outer restart loop */

    /* Final true physical divergence verification */
    if (compressible)
        assemble_div_rho_u(E, V, r0, lev);
    else
        assemble_div_u(E, V, r0, lev);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, r0) / (1e-32 + E->monitor.vdotv));

    *steps_max = count;
    return;
}

/* Solve compressible Stokes flow using
 * bi-conjugate gradient stablized (BiCG-stab) iterations
 */

static void solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                   double **V, double **P, double **FF,
                                   double imp, int *steps_max)
{
    /* modern c++: void assemble_div_rho_u() prototype in header */
    /* modern c++: void assemble_del2_u() prototype in header */
    /* modern c++: void assemble_grad_p() prototype in header */
    /* modern c++: void strip_bcs_from_residual() prototype in header */
    /* modern c++: int solve_del2_u() prototype in header */
    /* modern c++: void parallel_process_termination() prototype in header */

    /* modern c++: double global_pdot() prototype in header */
    /* modern c++: global_v_norm2(), global_p_norm2(), global_div_norm2() prototypes in header */
    /* modern c++: double CPU_time0() prototype in header */

    int npno, neq;
    int m, j, count, lev;
    int valid;

    double alpha, beta, omega,inner_imp;
    double r0dotrt, r1dotrt;
    double v_norm, p_norm;
    double dvelocity, dpressure;
    int converging;

    double *F[NCS];
    double *r1[NCS], *r2[NCS], *pt[NCS], *p1[NCS], *p2[NCS];
    double *rt[NCS], *v0[NCS], *s0[NCS], *st[NCS], *t0[NCS];
    double *u0[NCS];
    double *shuffle[NCS];

    double time0, v_res;
    
    inner_imp = imp * E->control.inner_accuracy_scale; /* allow for different innner loop accuracy */

    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    static std::vector<double> ws_F[NCS], ws_r1[NCS], ws_r2[NCS], ws_pt[NCS], ws_p1[NCS], ws_p2[NCS];
    static std::vector<double> ws_rt[NCS], ws_v0[NCS], ws_s0[NCS], ws_st[NCS], ws_t0[NCS], ws_u0[NCS];
    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[m]  = ensure_buf(ws_F[m],  (size_t)neq);
        r1[m] = ensure_buf(ws_r1[m], (size_t)(npno+1));
        r2[m] = ensure_buf(ws_r2[m], (size_t)(npno+1));
        pt[m] = ensure_buf(ws_pt[m], (size_t)(npno+1));
        p1[m] = ensure_buf(ws_p1[m], (size_t)(npno+1));
        p2[m] = ensure_buf(ws_p2[m], (size_t)(npno+1));
        rt[m] = ensure_buf(ws_rt[m], (size_t)(npno+1));
        v0[m] = ensure_buf(ws_v0[m], (size_t)(npno+1));
        s0[m] = ensure_buf(ws_s0[m], (size_t)(npno+1));
        st[m] = ensure_buf(ws_st[m], (size_t)(npno+1));
        t0[m] = ensure_buf(ws_t0[m], (size_t)(npno+1));
        u0[m] = ensure_buf(ws_u0[m], (size_t)neq);
    }

    time0 = CPU_time0();
    count = 0;
    v_res = E->monitor.fdotf;

    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[m][j] = FF[m][j];


    /* calculate the initial velocity residual */
    initial_vel_residual(E, V, P, F, inner_imp*v_res);


    /* initial residual r1 = div(rho_ref*V) */
    assemble_div_rho_u(E, V, r1, lev);

    E->monitor.vdotv = global_v_norm2(E, V);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, r1)
                                        / (1e-32 + E->monitor.vdotv));

    v_norm = sqrt(E->monitor.vdotv);
    p_norm = sqrt(E->monitor.pdotp);
    dvelocity = 1.0;
    dpressure = 1.0;
    converging = 0;


    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0,
                                   v_norm, p_norm,
                                   dvelocity, dpressure,
                                   E->monitor.incompressibility);
    }


    /* initial conjugate residual rt = r1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(j=1; j<=npno; j++)
            rt[m][j] = r1[m][j];


    valid = 1;
    r0dotrt = alpha = omega = 0;

    while( (count < *steps_max) && keep_iterating(E, imp, converging) ) {
        /* require two consecutive converging iterations to quit the while-loop */

        /* r1dotrt = <r1, rt> */
        r1dotrt = global_pdot(E, r1, rt, lev);
        if(r1dotrt == 0.0) {
            /* XXX: can we resume the computation when BiCGstab failed? */
            fprintf(E->fp, "BiCGstab method failed!!\n");
            fprintf(stderr, "BiCGstab method failed!!\n");
            parallel_process_termination();
        }


        /* update search direction */
        if(count == 0)
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p2[m][j] = r1[m][j];
        else {
            /* p2 = r1 + <r1,rt>/<r0,rt> * alpha/omega * (p1 - omega*v0) */
            beta = (r1dotrt / r0dotrt) * (alpha / omega);
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p2[m][j] = r1[m][j] + beta
                        * (p1[m][j] - omega * v0[m][j]);
        }


        /* preconditioner BPI ~= inv(K), pt = BPI*p2 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                pt[m][j] = E->BPI[lev][m][j] * p2[m][j];


        /* solve K*u0 = grad(pt) for u1 */
        assemble_grad_p(E, pt, F, lev);
        valid = solve_del2_u(E, u0, F, inner_imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
        }
        strip_bcs_from_residual(E, u0, lev);


        /* v0 = div(rho_ref*u0) */
        assemble_div_rho_u(E, u0, v0, lev);


        /* alpha = r1dotrt / <rt, v0> */
        alpha = r1dotrt / global_pdot(E, rt, v0, lev);


        /* s0 = r1 - alpha * v0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                s0[m][j] = r1[m][j] - alpha * v0[m][j];


        /* preconditioner BPI ~= inv(K), st = BPI*s0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                st[m][j] = E->BPI[lev][m][j] * s0[m][j];


        /* solve K*u1 = grad(st) for u1 */
        assemble_grad_p(E, st, F, lev);
        valid = solve_del2_u(E, E->u1, F, inner_imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 2\n", stderr);
            fputs("Warning: solver not converging! 2\n", E->fp);
        }
        strip_bcs_from_residual(E, E->u1, lev);


        /* t0 = div(rho_ref * u1) */
        assemble_div_rho_u(E, E->u1, t0, lev);


        /* omega = <t0, s0> / <t0, t0> */
        omega = global_pdot(E, t0, s0, lev) / global_pdot(E, t0, t0, lev);


        /* r2 = s0 - omega * t0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                r2[m][j] = s0[m][j] - omega * t0[m][j];


        /* P = P + alpha * pt + omega * st */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                s0[m][j] = alpha * pt[m][j] + omega * st[m][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                P[m][j] += s0[m][j];


        /* V = V - alpha * u0 - omega * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                F[m][j] = alpha * u0[m][j] + omega * E->u1[m][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[m][j] -= F[m][j];


        /* compute velocity, pressure, delta, and incompressibility residuals in 1 fused MPI_Allreduce */
        assemble_div_rho_u(E, V, t0, lev);
        {
            double norms[5];
            global_stokes_norms(E, V, P, F, s0, t0, norms);
            E->monitor.vdotv = norms[0];
            E->monitor.pdotp = norms[1];
            v_norm = sqrt(norms[0]);
            p_norm = sqrt(norms[1]);
            dvelocity = sqrt(norms[2] / (1e-32 + norms[0]));
            dpressure = sqrt(norms[3] / (1e-32 + norms[1]));
            E->monitor.incompressibility = sqrt(norms[4] / (1e-32 + norms[0]));
        }




        count++;

        if(E->control.print_convergence && E->parallel.me==0) {
            print_convergence_progress(E, count, time0,
                                       v_norm, p_norm,
                                       dvelocity, dpressure,
                                       E->monitor.incompressibility);
        }

	if(!valid){
            /* reset consecutive converging iterations */
            converging = 0;
	}else{
            /* how many consecutive converging iterations? */
            if(E->control.check_pressure_convergence) {
                /* check dv and dp */
                if(dvelocity < imp && dpressure < imp)
                    converging++;
                else
                    converging = 0;
            }else{
                /* check dv only */
                if(dvelocity < imp)
                    converging++;
                else
                    converging = 0;
            }
	  
	}

	/* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[m] = p1[m];
            p1[m] = p2[m];
            p2[m] = shuffle[m];

            shuffle[m] = r1[m];
            r1[m] = r2[m];
            r2[m] = shuffle[m];
        }

        /* shift <r0, rt> = <r1, rt> */
        r0dotrt = r1dotrt;

    } /* end loop for conjugate gradient */


    /* Persistent workspace retained across timesteps */
    *steps_max=count;

    return;

}


/* Solve compressible Stokes flow using
 * conjugate gradient (CG) iterations with an outer iteration
 */

static void solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                     double **V, double **P, double **F,
                                     double imp, int *steps_max)
{
    int m, i;
    int cycles, num_of_loop;
    double relative_err_v, relative_err_p;
    double *old_v[NCS], *old_p[NCS],*diff_v[NCS],*diff_p[NCS];
    double div_res;
    const int npno = E->lmesh.npno;
    const int neq = E->lmesh.neq;
    const int lev = E->mesh.levmax;

    /* modern c++: global_v_norm2(), global_p_norm2() prototypes in header */
    /* modern c++: double global_div_norm2() prototype in header */
    /* modern c++: void assemble_div_rho_u() prototype in header */
    
    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	old_v[m] = (double *)malloc(neq*sizeof(double));
    	diff_v[m] = (double *)malloc(neq*sizeof(double));
    	old_p[m] = (double *)malloc((npno+1)*sizeof(double));
    	diff_p[m] = (double *)malloc((npno+1)*sizeof(double));
    }

    cycles = E->control.p_iterations;

    initial_vel_residual(E, V, P, F,
                         imp * E->control.inner_accuracy_scale * E->monitor.fdotf);

    div_res = 1.0;
    relative_err_v = 1.0;
    relative_err_p = 1.0;
    num_of_loop = 0;

    while((relative_err_v >= imp || relative_err_p >= imp) &&
          (div_res > imp) &&
          (num_of_loop <= E->control.compress_iter_maxstep)) {

        for (m=1;m<=E->sphere.caps_per_proc;m++) {
            for(i=0;i<neq;i++) old_v[m][i] = V[m][i];
            for(i=1;i<=npno;i++) old_p[m][i] = P[m][i];
        }

        if (strcmp(E->control.uzawa, "fcg") == 0)
            solve_Ahat_p_fhat_FCG(E, V, P, F, imp, &cycles);
        else
            solve_Ahat_p_fhat_CG(E, V, P, F, imp, &cycles);

        /* compute norm of div(rho*V) */
        assemble_div_rho_u(E, V, E->u1, lev);
        div_res = sqrt(global_div_norm2(E, E->u1) / (1e-32 + E->monitor.vdotv));

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=0;i<neq;i++) diff_v[m][i] = V[m][i] - old_v[m][i];

        relative_err_v = sqrt( global_v_norm2(E,diff_v) /
                               (1.0e-32 + E->monitor.vdotv) );

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=npno;i++) diff_p[m][i] = P[m][i] - old_p[m][i];

        relative_err_p = sqrt( global_p_norm2(E,diff_p) /
                               (1.0e-32 + E->monitor.pdotp) );

        if(E->parallel.me == 0) {
            fprintf(stderr, "itercg -- div(rho*v)/v=%.2e dv/v=%.2e and dp/p=%.2e loop %d\n\n", div_res, relative_err_v, relative_err_p, num_of_loop);
            fprintf(E->fp, "itercg -- div(rho*v)/v=%.2e dv/v=%.2e and dp/p=%.2e loop %d\n\n", div_res, relative_err_v, relative_err_p, num_of_loop);
        }

        num_of_loop++;

    } /* end of while */

    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	free((void *) old_v[m]);
    	free((void *) old_p[m]);
	free((void *) diff_v[m]);
	free((void *) diff_p[m]);
    }

    return;
}


static void initial_vel_residual(struct All_variables *E,
                                 double **V, double **P, double **F,
                                 double acc)
{
    /* modern c++: void assemble_del2_u() prototype in header */
    /* modern c++: void assemble_grad_p() prototype in header */
    /* modern c++: void strip_bcs_from_residual() prototype in header */
    /* modern c++: int solve_del2_u() prototype in header */

    int neq = E->lmesh.neq;
    int lev = E->mesh.levmax;
    int i, m, valid;

    /* F = F - grad(P) - K*V */
    assemble_grad_p(E, P, E->u1, lev);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[m][i] = F[m][i] - E->u1[m][i];

    assemble_del2_u(E, V, E->u1, lev, 1);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[m][i] = F[m][i] - E->u1[m][i];

    strip_bcs_from_residual(E, F, lev);


    /* solve K*u1 = F for u1 (always use full tolerance-driven solve for initial velocity residual) */
    {
        const int saved_inner_mg = E->control.inner_mg_cycles;
        E->control.inner_mg_cycles = 0;
        valid = solve_del2_u(E, E->u1, F, acc, lev);
        E->control.inner_mg_cycles = saved_inner_mg;
    }
    if(!valid && (E->parallel.me==0)) {
        fputs("Warning: solver not converging! 0\n", stderr);
        fputs("Warning: solver not converging! 0\n", E->fp);
    }
    strip_bcs_from_residual(E, E->u1, lev);


    /* V = V + u1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            V[m][i] += E->u1[m][i];

    return;
}
