/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <sys/file.h>
#include <unistd.h>
#include "global_defs.h"
#include "checkpoints.h"
#include "composition_related.h"

/* Private function prototypes */
static void backup_file(const char *output_file);
static void write_sentinel(FILE *fp);
static void read_sentinel(FILE *fp, int me);

static void general_checkpoint(struct All_variables *E, FILE *fp);
static void tracer_checkpoint(struct All_variables *E, FILE *fp);
static void composition_checkpoint(struct All_variables *E, FILE *fp);
static void energy_checkpoint(struct All_variables *E, FILE *fp);
static void momentum_checkpoint(struct All_variables *E, FILE *fp);

static void read_general_checkpoint(struct All_variables *E, FILE *fp);
static void read_tracer_checkpoint(struct All_variables *E, FILE *fp);
static void read_composition_checkpoint(struct All_variables *E, FILE *fp);
static void read_energy_checkpoint(struct All_variables *E, FILE *fp);
static void read_momentum_checkpoint(struct All_variables *E, FILE *fp);

void myerror(struct All_variables *, char *);

static void prune_old_checkpoints(struct All_variables *E, int current_cycle);

void output_checkpoint(struct All_variables *E)
{
    void parallel_process_sync(struct All_variables *E);
    char output_file[255];
    FILE *fp1;

    sprintf(output_file, "%s.chkpt.%d.%d", E->control.data_file,
            E->parallel.me, E->monitor.solution_cycles);

    fp1 = fopen(output_file, "wb");
    if (fp1 == NULL) {
        fprintf(stderr, "Error: cannot open checkpoint file %s for writing, me=%d\n", output_file, E->parallel.me);
        parallel_process_termination();
    }

    /* checkpoint for general information */
    general_checkpoint(E, fp1);

    /* checkpoint for energy equation */
    energy_checkpoint(E, fp1);

    /* checkpoint for momentum equation */
    momentum_checkpoint(E, fp1);

    /* checkpoint for tracer/composition */
    if(E->control.tracer) {
        tracer_checkpoint(E, fp1);

        if(E->composition.on)
            composition_checkpoint(E, fp1);
    }

    fclose(fp1);

    /* Two-phase commit: wait for all MPI processes to finish writing their checkpoint */
    parallel_process_sync(E);

    /* Rank 0 updates last_good_checkpoint and prunes old checkpoints */
    if (E->parallel.me == 0) {
        char good_file[255];
        FILE *fp_good;

        sprintf(good_file, "%s.last_good_checkpoint", E->control.data_file);
        fp_good = fopen(good_file, "w");
        if (fp_good != NULL) {
            fprintf(fp_good, "cycle %d\n", E->monitor.solution_cycles);
            fprintf(fp_good, "elapsed_time %.12e\n", (double)E->monitor.elapsed_time);
            fprintf(fp_good, "timestep %.12e\n", (double)E->advection.timestep);
            fprintf(fp_good, "nprocs %d\n", E->parallel.nproc);
            fclose(fp_good);
        }

        if (E->control.checkpoints_to_keep > 0) {
            prune_old_checkpoints(E, E->monitor.solution_cycles);
        }
    }

    parallel_process_sync(E);
    return;
}


void read_checkpoint(struct All_variables *E)
{
    void initialize_material(struct All_variables *E);
    void initial_viscosity(struct All_variables *E);

    char output_file[255];
    FILE *fp;

    /* open the checkpoint file */
    snprintf(output_file, 254, "%s.chkpt.%d.%d", E->control.old_P_file,
             E->parallel.me, E->monitor.solution_cycles_init);
    fp = fopen(output_file, "rb");
    if(fp == NULL) {
        fprintf(stderr, "Cannot open file: %s\n", output_file);
        exit(-1);
    }
    if(E->parallel.me == 0)
      fprintf(stderr,"read_checkpoint: restarting from %s\n",output_file);
	
    /* check mesh information in the checkpoint file */
    read_general_checkpoint(E, fp);

    /* init E->mat */
    initialize_material(E);

    /* read energy information in the checkpoint file */
    read_energy_checkpoint(E, fp);

    /* read momentum information in the checkpoint file */
    read_momentum_checkpoint(E, fp);

    /* read tracer/composition/strain information in the checkpoint file */
    if(E->control.tracer) {
      if(E->trace.ic_method_for_flavors == 99){
	if(E->parallel.me == 0)
	  fprintf(stderr,"ic_method_for_flavors = 99 will override checkpoint restart\n");
      }else{
        read_tracer_checkpoint(E, fp);
        if(E->composition.on)
            read_composition_checkpoint(E, fp);
      }
    }

    fclose(fp);

    /* finally, init viscosity */
    initial_viscosity(E);

    /* update next output and checkpoint times based on restored elapsed_time */
    update_next_output_times(E);

    return;
}


static void backup_file(const char *output_file)
{
    char bak_file[255];
    int ierr;

    /* check the existence of output_file */
    if(access(output_file, F_OK) == 0) {
        /* if exist, renamed it to back up */
        sprintf(bak_file, "%s.bak", output_file);
        ierr = rename(output_file, bak_file);
        if(ierr != 0) {
            fprintf(stderr, "Warning, cannot backup checkpoint files\n");
        }
    }

    return;
}


static void write_sentinel(FILE *fp)
{
    int a[4] = {0, 0, 0, 0};

    fwrite(a, sizeof(int), 4, fp);
}


static void read_sentinel(FILE *fp, int me)
{
    int i, a[4];
    int nonzero = 0;

    fread(a, sizeof(int), 4, fp);

    /* check whether a[i] are all zero */
    for(i=0; i<4; i++)
        nonzero |= a[i];

    if(nonzero) {
        fprintf(stderr, "Error in reading checkpoint file: wrong sentinel, "
                "me=%d\n", me);
        exit(-1);
    }

    return;
}


static void general_checkpoint(struct All_variables *E, FILE *fp)
{
    /* write mesh information */
    fwrite(&(E->lmesh.nox), sizeof(int), 1, fp);
    fwrite(&(E->lmesh.noy), sizeof(int), 1, fp);
    fwrite(&(E->lmesh.noz), sizeof(int), 1, fp);
    fwrite(&(E->parallel.nprocx), sizeof(int), 1, fp);
    fwrite(&(E->parallel.nprocy), sizeof(int), 1, fp);
    fwrite(&(E->parallel.nprocz), sizeof(int), 1, fp);
    fwrite(&(E->sphere.caps_per_proc), sizeof(int), 1, fp);

    /* write timing information */
    fwrite(&(E->monitor.solution_cycles), sizeof(int), 1, fp);
    fwrite(&(E->monitor.elapsed_time), sizeof(float), 1, fp);
    fwrite(&(E->advection.timestep), sizeof(float), 1, fp);
    fwrite(&(E->control.start_age), sizeof(float), 1, fp);

    return;
}


static void read_general_checkpoint(struct All_variables *E, FILE *fp)
{
    int tmp[7];
    double dtmp;

    /* read mesh information */
    fread(tmp, sizeof(int), 7, fp);

    if((tmp[0] != E->lmesh.nox) ||
       (tmp[1] != E->lmesh.noy) ||
       (tmp[2] != E->lmesh.noz) ||
       (tmp[3] != E->parallel.nprocx) ||
       (tmp[4] != E->parallel.nprocy) ||
       (tmp[5] != E->parallel.nprocz) ||
       (tmp[6] != E->sphere.caps_per_proc)) {

        fprintf(stderr, "Error in reading checkpoint file: mesh parameters mismatch, me=%d\n",
                E->parallel.me);
        fprintf(stderr, "%d %d %d %d %d %d %d\n",
                tmp[0], tmp[1], tmp[2], tmp[3],
                tmp[4], tmp[5], tmp[6]);
        exit(-1);
    }

    /* read timing information */
    tmp[0] = fread(&(E->monitor.solution_cycles), sizeof(int), 1, fp);
    tmp[0]+= fread(&(E->monitor.elapsed_time), sizeof(float), 1, fp);
    tmp[0]+= fread(&(E->advection.timestep), sizeof(float), 1, fp);
    tmp[0]+= fread(&(E->control.start_age), sizeof(float), 1, fp);
    if(tmp[0] != 4)
      myerror(E,"read_general_checkpoint: header error");

    E->advection.timesteps = E->monitor.solution_cycles;

    return;
}


static void tracer_checkpoint(struct All_variables *E, FILE *fp)
{
    int m, i;

    write_sentinel(fp);

    fwrite(&(E->trace.number_of_basic_quantities), sizeof(int), 1, fp);
    fwrite(&(E->trace.number_of_extra_quantities), sizeof(int), 1, fp);
    fwrite(&(E->trace.nflavors), sizeof(int), 1, fp);
    fwrite(&(E->trace.ilast_tracer_count), sizeof(int), 1, fp);

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        fwrite(&(E->trace.ntracers[m]), sizeof(int), 1, fp);

    /* the 0-th element of basicq/extraq/ielement is not init'd
     * and won't be used when read it. */
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(i=0; i<6; i++) {
            fwrite(E->trace.basicq[m][i], sizeof(double),
                   E->trace.ntracers[m]+1, fp);
        }
        for(i=0; i<E->trace.number_of_extra_quantities; i++) {
            fwrite(E->trace.extraq[m][i], sizeof(double),
                   E->trace.ntracers[m]+1, fp);
        }
        fwrite(E->trace.ielement[m], sizeof(int),
               E->trace.ntracers[m]+1, fp);
    }

    return;
}


static void read_tracer_checkpoint(struct All_variables *E, FILE *fp)
{
    void count_tracers_of_flavors(struct All_variables *E);
    void allocate_tracer_arrays();
    void init_tracer_strain();
    void init_strain();

    int m, i, itmp;

    read_sentinel(fp, E->parallel.me);

    fread(&itmp, sizeof(int), 1, fp);
    if (itmp != E->trace.number_of_basic_quantities) {
        fprintf(stderr, "Error in reading checkpoint file: tracer basicq, me=%d\n",
                E->parallel.me);
        fprintf(stderr, "%d\n", itmp);
        exit(-1);

    }

    fread(&itmp, sizeof(int), 1, fp);
    if (itmp != E->trace.number_of_extra_quantities) {
        fprintf(stderr, "Error in reading checkpoint file: tracer extraq, me=%d\n",
                E->parallel.me);
        fprintf(stderr, "%d\n", itmp);
        exit(-1);

    }

    fread(&itmp, sizeof(int), 1, fp);
    if (itmp != E->trace.nflavors) {
        fprintf(stderr, "Error in reading checkpoint file: tracer nflavors, me=%d\n",
                E->parallel.me);
        fprintf(stderr, "%d\n", itmp);
        exit(-1);

    }

    fread(&itmp, sizeof(int), 1, fp);
    E->trace.ilast_tracer_count = itmp;

    /* # of tracers, allocate memory */
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        fread(&itmp, sizeof(int), 1, fp);
        allocate_tracer_arrays(E, m, itmp);
        E->trace.ntracers[m] = itmp;
    }

    /* read tracer data */
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(i=0; i<6; i++) {
            fread(E->trace.basicq[m][i], sizeof(double),
                  E->trace.ntracers[m]+1, fp);
        }
        for(i=0; i<E->trace.number_of_extra_quantities; i++) {
            fread(E->trace.extraq[m][i], sizeof(double),
                  E->trace.ntracers[m]+1, fp);
        }
        fread(E->trace.ielement[m], sizeof(int),
              E->trace.ntracers[m]+1, fp);
    }

    /* init E->trace.ntracer_flavor */
    count_tracers_of_flavors(E);

    /* In case strain needs to be erased */
    if(E->trace.erase_strain)
       init_tracer_strain(E);

    /* In case strain is tracked */
    /* Average strain from tracers to elements and map on nodes */
    if(E->trace.track_strain)
       init_strain(E);

    return;
}


static void composition_checkpoint(struct All_variables *E, FILE *fp)
{
    int i, m;

    write_sentinel(fp);

    fwrite(&(E->composition.ncomp), sizeof(int), 1, fp);
    fwrite(E->composition.bulk_composition, sizeof(double),
           E->composition.ncomp, fp);
    fwrite(E->composition.initial_bulk_composition, sizeof(double),
           E->composition.ncomp, fp);

    /* the 0-th element of comp_el is not init'd
     * and won't be used when read it. */
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(i=0; i<E->composition.ncomp; i++)
            fwrite(E->composition.comp_el[m][i], sizeof(double),
                   E->lmesh.nel+1, fp);
    }

    return;
}


static void read_composition_checkpoint(struct All_variables *E, FILE *fp)
{
    double tmp;
    int m, i, itmp;

    read_sentinel(fp, E->parallel.me);

    fread(&itmp, sizeof(int), 1, fp);
    if (itmp != E->composition.ncomp) {
        fprintf(stderr, "Error in reading checkpoint file: ncomp, me=%d\n",
                E->parallel.me);
        fprintf(stderr, "%d\n", itmp);
        exit(-1);
    }

    fread(E->composition.bulk_composition, sizeof(double),
          E->composition.ncomp, fp);

    fread(E->composition.initial_bulk_composition, sizeof(double),
          E->composition.ncomp, fp);

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(i=0; i<E->composition.ncomp; i++)
            fread(E->composition.comp_el[m][i], sizeof(double),
                  E->lmesh.nel+1, fp);
    }

    /* init E->composition.comp_node */
    map_composition_to_nodes(E);

    /* preventing uninitialized access */
    E->trace.istat_iempty = 0;

    for (i=0; i<E->composition.ncomp; i++) {
        E->composition.error_fraction[i] = E->composition.bulk_composition[i]
        / E->composition.initial_bulk_composition[i] - 1.0;
    }

    return;
}


static void energy_checkpoint(struct All_variables *E, FILE *fp)
{
    int m;

    write_sentinel(fp);

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        fwrite(E->T[m], sizeof(double), E->lmesh.nno+1, fp);
        fwrite(E->Tdot[m], sizeof(double), E->lmesh.nno+1, fp);
    }

    return;
}


static void read_energy_checkpoint(struct All_variables *E, FILE *fp)
{
    int m;

    read_sentinel(fp, E->parallel.me);

    /* the 0-th element of T/Tdot is not init'd
     * and won't be used when read it. */
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
      if(fread(E->T[m], sizeof(double), E->lmesh.nno+1, fp)!= E->lmesh.nno+1)
	myerror(E,"read_energy_checkpoint: error at T");
      if(fread(E->Tdot[m], sizeof(double), E->lmesh.nno+1, fp)!=E->lmesh.nno+1)
	myerror(E,"read_energy_checkpoint: error at Tdot");
    }

    return;
}


static void momentum_checkpoint(struct All_variables *E, FILE *fp)
{
    int m;
    float junk[2];
    junk[0] = junk[1] = 0;

    write_sentinel(fp);

    /* for backward compatibility */
    fwrite(junk, sizeof(float), 2, fp);

    /* the 0-th element of P/NP/EVI/VI is not init'd
     * and won't be used when read it. */
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        /* Pressure at equation points */
        fwrite(E->P[m], sizeof(double), E->lmesh.npno+1, fp);

        /* velocity at equation points */
        fwrite(E->U[m], sizeof(double), E->lmesh.neq, fp);
    }

    return;
}


static void read_momentum_checkpoint(struct All_variables *E, FILE *fp)
{
    void v_from_vector();
    void p_to_nodes();
    double global_v_norm2(), global_p_norm2();

    int m;
    int lev = E->mesh.levmax;
    float junk[2];

    read_sentinel(fp, E->parallel.me);

    /* for backward compatibility */
    if(fread(junk, sizeof(float), 2, fp)!=2)
      myerror(E,"read_momentum_checkpoint: error at vdotv");

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        /* Pressure at equation points */
      if(fread(E->P[m], sizeof(double), E->lmesh.npno+1, fp) !=  E->lmesh.npno+1)
	myerror(E,"read_momentum_checkpoint: error at P");
        /* velocity at equation points */
      if(fread(E->U[m], sizeof(double), E->lmesh.neq, fp) != E->lmesh.neq)
	myerror(E,"read_momentum_checkpoint: error at U");
    }

    E->monitor.vdotv = global_v_norm2(E, E->U);
    E->monitor.pdotp = global_p_norm2(E, E->P);

    /* update velocity array */
    v_from_vector(E);

    /* init E->NP */
    p_to_nodes(E, E->P, E->NP, lev);

    return;
}


/*=======================================================================
  Functions for two-stage continuation reading directly from binary .chkpt
=========================================================================*/

int get_checkpoint_elapsed_time(struct All_variables *E, const char *output_file, float *elapsed_time)
{
    FILE *fp;
    int tmp[7];
    int cycles;
    float time;

    fp = fopen(output_file, "rb");
    if (fp == NULL)
        return 0;

    if (fread(tmp, sizeof(int), 7, fp) != 7) {
        fclose(fp);
        return 0;
    }

    if (fread(&cycles, sizeof(int), 1, fp) != 1) {
        fclose(fp);
        return 0;
    }

    if (fread(&time, sizeof(float), 1, fp) != 1) {
        fclose(fp);
        return 0;
    }

    *elapsed_time = time;
    fclose(fp);
    return 1;
}


int read_tic_from_checkpoint(struct All_variables *E, const char *output_file)
{
    FILE *fp;
    int tmp[7];
    int m;

    fp = fopen(output_file, "rb");
    if (fp == NULL)
        return 0;

    /* read & check mesh parameters */
    if (fread(tmp, sizeof(int), 7, fp) != 7) {
        fclose(fp);
        return 0;
    }

    if((tmp[0] != E->lmesh.nox) ||
       (tmp[1] != E->lmesh.noy) ||
       (tmp[2] != E->lmesh.noz) ||
       (tmp[3] != E->parallel.nprocx) ||
       (tmp[4] != E->parallel.nprocy) ||
       (tmp[5] != E->parallel.nprocz) ||
       (tmp[6] != E->sphere.caps_per_proc)) {
        fprintf(stderr, "Error in reading checkpoint file: mesh parameters mismatch, me=%d\n",
                E->parallel.me);
        fprintf(stderr, "Expected: %d %d %d %d %d %d %d\n",
                E->lmesh.nox, E->lmesh.noy, E->lmesh.noz,
                E->parallel.nprocx, E->parallel.nprocy, E->parallel.nprocz,
                E->sphere.caps_per_proc);
        fprintf(stderr, "Found:    %d %d %d %d %d %d %d\n",
                tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6]);
        fclose(fp);
        myerror(E, "read_tic_from_checkpoint: mesh mismatch");
    }

    /* skip timing info: 1 int + 3 floats */
    if (fseek(fp, sizeof(int) + 3 * sizeof(float), SEEK_CUR) != 0) {
        fclose(fp);
        myerror(E, "read_tic_from_checkpoint: seek timing error");
    }

    /* read sentinel for energy */
    read_sentinel(fp, E->parallel.me);

    /* read T for each cap, skip Tdot */
    for(m = 1; m <= E->sphere.caps_per_proc; m++) {
        if(fread(E->T[m], sizeof(double), E->lmesh.nno + 1, fp) != E->lmesh.nno + 1) {
            fclose(fp);
            myerror(E, "read_tic_from_checkpoint: error reading T");
        }
        if(fseek(fp, (long)(E->lmesh.nno + 1) * sizeof(double), SEEK_CUR) != 0) {
            fclose(fp);
            myerror(E, "read_tic_from_checkpoint: error seeking past Tdot");
        }
    }

    fclose(fp);
    return 1;
}


int read_tracers_from_checkpoint(struct All_variables *E, const char *output_file)
{
    FILE *fp;
    int tmp[7];
    int check_sentinel[4];
    long offset;

    fp = fopen(output_file, "rb");
    if (fp == NULL)
        return 0;

    /* read & check mesh parameters */
    if (fread(tmp, sizeof(int), 7, fp) != 7) {
        fclose(fp);
        return 0;
    }

    if((tmp[0] != E->lmesh.nox) ||
       (tmp[1] != E->lmesh.noy) ||
       (tmp[2] != E->lmesh.noz) ||
       (tmp[3] != E->parallel.nprocx) ||
       (tmp[4] != E->parallel.nprocy) ||
       (tmp[5] != E->parallel.nprocz) ||
       (tmp[6] != E->sphere.caps_per_proc)) {
        fprintf(stderr, "Error reading checkpoint tracers: mesh mismatch, me=%d\n", E->parallel.me);
        fclose(fp);
        myerror(E, "read_tracers_from_checkpoint: mesh mismatch");
    }

    /* skip timing info: 1 int + 3 floats */
    if (fseek(fp, sizeof(int) + 3 * sizeof(float), SEEK_CUR) != 0) {
        fclose(fp);
        myerror(E, "read_tracers_from_checkpoint: seek timing error");
    }

    /* skip energy checkpoint: sentinel + 2 * (nno+1) doubles per cap */
    read_sentinel(fp, E->parallel.me);
    offset = (long)E->sphere.caps_per_proc * 2 * (E->lmesh.nno + 1) * sizeof(double);
    if (fseek(fp, offset, SEEK_CUR) != 0) {
        fclose(fp);
        myerror(E, "read_tracers_from_checkpoint: seek energy error");
    }

    /* skip momentum checkpoint: sentinel + 2 floats + ((npno+1) + neq) doubles per cap */
    read_sentinel(fp, E->parallel.me);
    offset = 2 * sizeof(float) + (long)E->sphere.caps_per_proc * ((E->lmesh.npno + 1) + E->lmesh.neq) * sizeof(double);
    if (fseek(fp, offset, SEEK_CUR) != 0) {
        fclose(fp);
        myerror(E, "read_tracers_from_checkpoint: seek momentum error");
    }

    /* verify tracer section exists (sentinel check preview) */
    if (fread(check_sentinel, sizeof(int), 4, fp) != 4) {
        fclose(fp);
        return 0;
    }
    fseek(fp, -4 * (long)sizeof(int), SEEK_CUR);

    /* read tracer checkpoint */
    if (E->trace.ic_method_for_flavors == 99) {
        if (E->parallel.me == 0)
            fprintf(stderr, "ic_method_for_flavors = 99 will override checkpoint restart\n");
    } else {
        read_tracer_checkpoint(E, fp);
        if (E->composition.on)
            read_composition_checkpoint(E, fp);
    }

    fclose(fp);
    return 1;
}


/*=======================================================================
  Sliding Window Checkpoint Pruning & Auto-Restart Implementation
=========================================================================*/

static void prune_old_checkpoints(struct All_variables *E, int current_cycle)
{
    char hist_file[255];
    FILE *fp;
    int cycles[1000];
    int count = 0;
    int i, r;
    char path[255];

    if (E->control.checkpoints_to_keep <= 0)
        return;

    snprintf(hist_file, sizeof(hist_file), "%s.checkpoint_history.txt", E->control.data_file);

    /* 1. Read existing cycle history */
    fp = fopen(hist_file, "r");
    if (fp != NULL) {
        while (count < 990 && fscanf(fp, "%d", &cycles[count]) == 1) {
            count++;
        }
        fclose(fp);
    }

    /* 2. Append current_cycle if not already the last entry */
    if (count == 0 || cycles[count - 1] != current_cycle) {
        cycles[count++] = current_cycle;
    }

    /* 3. If count exceeds checkpoints_to_keep, prune the oldest entries */
    while (count > E->control.checkpoints_to_keep) {
        int old_cycle = cycles[0];

        /* Delete checkpoint files for all ranks at old_cycle */
        for (r = 0; r < E->parallel.nproc; r++) {
            snprintf(path, sizeof(path), "%s.chkpt.%d.%d", E->control.data_file, r, old_cycle);
            remove(path);
        }

        fprintf(stderr, "Checkpoints: pruned old checkpoint files for cycle %d (keeping last %d checkpoints)\n",
                old_cycle, E->control.checkpoints_to_keep);

        /* Shift history array left */
        for (i = 0; i < count - 1; i++) {
            cycles[i] = cycles[i + 1];
        }
        count--;
    }

    /* 4. Write back updated history */
    fp = fopen(hist_file, "w");
    if (fp != NULL) {
        for (i = 0; i < count; i++) {
            fprintf(fp, "%d\n", cycles[i]);
        }
        fclose(fp);
    }
}


int determine_auto_restart_checkpoint(struct All_variables *E)
{
    char good_file[255];
    FILE *fp;
    int cycle = -1;
    char line[100];
    int found = 0;

    if (E->parallel.me == 0) {
        /* 1. Check last_good_checkpoint file */
        snprintf(good_file, sizeof(good_file), "%s.last_good_checkpoint", E->control.old_P_file);
        fp = fopen(good_file, "r");
        if (fp != NULL) {
            while (fgets(line, sizeof(line), fp)) {
                if (sscanf(line, "cycle %d", &cycle) == 1) {
                    found = 1;
                    break;
                }
            }
            fclose(fp);
        }

        /* 2. Fallback: check checkpoint_history.txt */
        if (!found) {
            char hist_file[255];
            snprintf(hist_file, sizeof(hist_file), "%s.checkpoint_history.txt", E->control.old_P_file);
            fp = fopen(hist_file, "r");
            if (fp != NULL) {
                int last_c = -1, c;
                while (fscanf(fp, "%d", &c) == 1) {
                    last_c = c;
                }
                fclose(fp);
                if (last_c >= 0) {
                    cycle = last_c;
                    found = 1;
                }
            }
        }

        /* 3. If cycle identified, verify checkpoint rank 0 file exists */
        if (found && cycle >= 0) {
            char chkpt_test[255];
            snprintf(chkpt_test, sizeof(chkpt_test), "%s.chkpt.0.%d", E->control.old_P_file, cycle);
            fp = fopen(chkpt_test, "rb");
            if (fp == NULL) {
                fprintf(stderr, "Auto-restart: candidate cycle %d found but %s is missing!\n",
                        cycle, chkpt_test);
                found = 0;
                cycle = -1;
            } else {
                fclose(fp);
            }
        }
    }

    /* Broadcast decision to all ranks */
    MPI_Bcast(&found, 1, MPI_INT, 0, E->parallel.world);
    MPI_Bcast(&cycle, 1, MPI_INT, 0, E->parallel.world);

    if (found && cycle >= 0) {
        E->control.restart = 1;
        E->monitor.solution_cycles_init = cycle;
        if (E->parallel.me == 0) {
            fprintf(stderr, "Auto-restart: successfully detected checkpoint at cycle %d. Resuming computation.\n", cycle);
            fprintf(E->fp, "Auto-restart: successfully detected checkpoint at cycle %d. Resuming computation.\n", cycle);
        }
        return 1;
    } else {
        E->control.restart = 0;
        E->monitor.solution_cycles_init = 0;
        if (E->parallel.me == 0) {
            fprintf(stderr, "Auto-restart: no valid checkpoint found. Starting initial run from cycle 0.\n");
            fprintf(E->fp, "Auto-restart: no valid checkpoint found. Starting initial run from cycle 0.\n");
        }
        return 0;
    }
}


void update_next_output_times(struct All_variables *E)
{
    if (E->control.storage_interval_time > 0.0) {
        double interval = E->control.storage_interval_time;
        int k = (int)(E->monitor.elapsed_time / interval);
        E->control.next_storage_time = (k + 1) * interval;
        if (E->control.next_storage_time <= E->monitor.elapsed_time + 1.0e-4 * interval)
            E->control.next_storage_time += interval;
    } else {
        E->control.next_storage_time = 1.0e20;
    }

    if (E->control.checkpoint_interval_time > 0.0) {
        double interval = E->control.checkpoint_interval_time;
        int k = (int)(E->monitor.elapsed_time / interval);
        E->control.next_checkpoint_time = (k + 1) * interval;
        if (E->control.next_checkpoint_time <= E->monitor.elapsed_time + 1.0e-4 * interval)
            E->control.next_checkpoint_time += interval;
    } else {
        E->control.next_checkpoint_time = 1.0e20;
    }
}
