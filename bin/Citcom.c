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
#include <mpi.h>

#include <math.h>
#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"
#include "citcom_init.h"
#include "interuption.h"
#include "output.h"
#include "parallel_related.h"
#include "checkpoints.h"

extern int Emergency_stop;

void solver_init(struct All_variables *E);

int main(argc,argv)
     int argc;
     char **argv;

{	/* Functions called by main*/
  void general_stokes_solver();
  void general_stokes_solver_pseudo_surf();
  void global_default_values();
  void read_instructions();
  void initial_setup();
  void initial_conditions();
  void post_processing();
  void read_velocity_boundary_from_file();
  void read_rayleigh_from_file();
  void read_mat_from_file();
  void read_temperature_boundary_from_file();

  void output_finalize();
  void tracer_advection();
  void heat_flux();

  float cpu_time_on_vp_it;

  int cpu_total_seconds,k,need_init_sol;
  double CPU_time0(),time,initial_time,start_time;
  int last_saved_checkpoint_cycle;
  int do_field_output, do_checkpoint;

  struct All_variables *E;
  MPI_Comm world;

  MPI_Init(&argc,&argv); /* added here to allow command-line input */

  if (argc < 2)   {
    fprintf(stderr,"Usage: %s PARAMETERFILE\n", argv[0]);
    parallel_process_termination();
  }



  /* this section reads input, allocates memory, and set some initial values;
   * replaced by CitcomS.Controller.initialize() and
   * CitcomS.Solver.initialize() in Pyre. */
  world = MPI_COMM_WORLD;
  E = citcom_init(&world); /* allocate global E and do initializaion here */

  /* define common aliases for full/regional functions */
  solver_init(E);

  start_time = time = CPU_time0();

  /* Global interuption handling routine defined once here */
  set_signal();

  /* default values for various parameters */
  global_default_values(E);

  /* read input parameters from file */
  read_instructions(E, argv[1]);

  /* create mesh, setup solvers etc. */
  initial_setup(E);

  cpu_time_on_vp_it = CPU_time0();
  initial_time = cpu_time_on_vp_it - time;
  if (E->parallel.me == 0)  {
    fprintf(stderr,"Input parameters taken from file '%s'\n",argv[1]);
    fprintf(stderr,"Initialization complete after %g seconds\n\n",initial_time);
    fprintf(E->fp,"Initialization complete after %g seconds\n\n",initial_time);
    fflush(E->fp);
  }

  /* write all config parameters to a file named pidXXXXXXXXX */
  //print_all_config_parameters(E);



  /* this section sets the initial condition;
   * replaced by CitcomS.Controller.launch() ->
   * CitcomS.Solver.launch() in Pyre. */
  if (E->control.post_p) {
      /* the initial condition is from previous checkpoint */
      read_checkpoint(E);

      /* the program will finish after post_processing */
      post_processing(E);
      (E->problem_output)(E, E->monitor.solution_cycles);

      citcom_finalize(E, 0);
  }

  if (E->control.restart) {
      /* the initial condition is from previous checkpoint */
      read_checkpoint(E);
      if(E->control.tracer && (E->trace.ic_method_for_flavors == 99)){
	/* 
	   if ggrd tracer input is selected, this will override
	   existing tracers, or allow addition of tracers after
	   restart from a thermal model
	 
	*/
	initialize_tracers(E);
        if (E->composition.on)
	  init_composition(E);
	need_init_sol = 1;
      }else
	need_init_sol = 0;
  }
  else {
      /* regular init, or read T from file only */

      initial_conditions(E);
      need_init_sol = 1;	/*  */
  }
  if(need_init_sol){
    /* find first solution */
      if(E->control.pseudo_free_surf) {
          if(E->mesh.topvbc == 2)
              general_stokes_solver_pseudo_surf(E);
          else
              assert(0);
      }
      else
          general_stokes_solver(E);
  }

  /* stop the computation if only computes stokes' problem */
  if (E->control.stokes)  {
    if(E->control.tracer==1)
      tracer_advection(E);
    (E->problem_output)(E, E->monitor.solution_cycles);
    citcom_finalize(E, 0);
  }


  (E->problem_output)(E, E->monitor.solution_cycles);

  /* information about simulation time and wall clock time */
  output_time(E, E->monitor.solution_cycles);

  E->control.last_checkpoint_walltime = CPU_time0();
  last_saved_checkpoint_cycle = E->monitor.solution_cycles;

  if(!E->control.restart)	/* if we have not restarted, print new
				   checkpoint, else leave as is to
				   allow reusing directories */
    output_checkpoint(E);


 

  /* this section advances the time step;
   * replaced by CitcomS.Controller.march() in Pyre. */
  while ( E->control.keep_going   &&  (Emergency_stop == 0) ) {


    
    /* The next few lines of code were replaced by
     * pyCitcom_PG_timestep_solve() in Pyre version.
     * If you modify here, make sure its Pyre counterpart
     * is modified as well */
    E->monitor.solution_cycles++;
    if(E->monitor.solution_cycles>E->control.print_convergence)
      E->control.print_convergence=1;

    (E->next_buoyancy_field)(E);
    /* */


    if(((E->advection.total_timesteps < E->advection.max_total_timesteps) &&
	(E->advection.timesteps < E->advection.max_timesteps)) ||
       (E->advection.total_timesteps < E->advection.min_timesteps) )
      E->control.keep_going = 1;
    else
      E->control.keep_going = 0;

    cpu_total_seconds = CPU_time0()-start_time;
    if (cpu_total_seconds > E->control.record_all_until)  {
      E->control.keep_going = 0;
    }

    if (E->monitor.T_interior > E->monitor.T_interior_max_for_exit)  {
      fprintf(E->fp,"quit due to maxT = %.4e sub_iteration%d\n",E->monitor.T_interior,E->advection.last_sub_iterations);
      parallel_process_termination();
    }

    if(E->control.tracer==1)
      tracer_advection(E);
    general_stokes_solver(E);

    if(E->output.write_q_files)
      if ((E->monitor.solution_cycles % E->output.write_q_files)==0)
	heat_flux(E);
    /* 1. Field output triggering (Step-based or Model Time/Age-based) */
    do_field_output = 0;
    if (E->control.record_every > 0 && (E->monitor.solution_cycles % E->control.record_every) == 0) {
        do_field_output = 1;
    }
    if (E->control.storage_interval_time > 0.0 &&
        E->monitor.elapsed_time >= E->control.next_storage_time - 1.0e-4 * E->control.storage_interval_time) {
        do_field_output = 1;
        while (E->control.next_storage_time <= E->monitor.elapsed_time + 1.0e-4 * E->control.storage_interval_time) {
            E->control.next_storage_time += E->control.storage_interval_time;
        }
    }
    if (do_field_output) {
	(E->problem_output)(E, E->monitor.solution_cycles);
    }

    /* information about simulation time and wall clock time */
    output_time(E, E->monitor.solution_cycles);

    /* 2. Checkpoint triggering (Step-based, Model Time/Age-based, or Wall-clock-based) */
    do_checkpoint = 0;
    if (E->control.checkpoint_frequency > 0 && (E->monitor.solution_cycles % E->control.checkpoint_frequency) == 0) {
        do_checkpoint = 1;
    }
    if (E->control.checkpoint_interval_time > 0.0 &&
        E->monitor.elapsed_time >= E->control.next_checkpoint_time - 1.0e-4 * E->control.checkpoint_interval_time) {
        do_checkpoint = 1;
        while (E->control.next_checkpoint_time <= E->monitor.elapsed_time + 1.0e-4 * E->control.checkpoint_interval_time) {
            E->control.next_checkpoint_time += E->control.checkpoint_interval_time;
        }
    }
    if (E->control.checkpoint_walltime > 0.0 &&
        (CPU_time0() - E->control.last_checkpoint_walltime >= E->control.checkpoint_walltime)) {
        do_checkpoint = 1;
    }

    if (do_checkpoint &&
        ((!E->control.restart) || (E->monitor.solution_cycles != E->monitor.solution_cycles_init))) {
        output_checkpoint(E);
        E->control.last_checkpoint_walltime = CPU_time0();
        last_saved_checkpoint_cycle = E->monitor.solution_cycles;
    }
    /* updating time-dependent material group
     * if mat_control is 0, the material group has already been
     * initialized in initial_conditions() */
    if(E->control.mat_control==1)
      read_mat_from_file(E);

#ifdef USE_GGRD
    /* updating local rayleigh number (based on Netcdf grds, the
       rayleigh number may be modified laterally in the surface
       layers) */
    /* no counterpart in pyre */
    if(E->control.ggrd.ray_control)
      read_rayleigh_from_file(E);
    if(E->control.ggrd.vtop_control)
      ggrd_read_vtop_from_file(E, TRUE);

#endif

    /* updating plate velocity boundary condition */
    if(E->control.vbcs_file==1)
      read_velocity_boundary_from_file(E);

    /* updating plate temperature boundary condition */
    if(E->control.tbcs_file)
      read_temperature_boundary_from_file(E);   


    if (E->parallel.me == 0)  {
      fprintf(E->fp,"CPU total = %g & CPU = %g for step %d time = %.4e dt = %.4e  maxT = %.4e sub_iteration%d\n",
	      CPU_time0()-start_time,CPU_time0()-time,E->monitor.solution_cycles,E->monitor.elapsed_time,
	      E->advection.timestep,E->monitor.T_interior,E->advection.last_sub_iterations);

      time = CPU_time0();
    }

  }


  /* Checkpoint on Termination */
  if (E->control.checkpoint_on_termination) {
    if (last_saved_checkpoint_cycle != E->monitor.solution_cycles) {
      if (E->parallel.me == 0)
        fprintf(stderr, "Checkpoint on termination: saving final checkpoint at cycle %d\n", E->monitor.solution_cycles);
      output_checkpoint(E);
      last_saved_checkpoint_cycle = E->monitor.solution_cycles;
    }
  }

  /* this section prints time accounting;
   * no counterpart in pyre */
  if (E->parallel.me == 0)  {
    fprintf(stderr,"cycles=%d\n",E->monitor.solution_cycles);
    cpu_time_on_vp_it=CPU_time0()-cpu_time_on_vp_it;
    fprintf(stderr,"Average cpu time taken for velocity step = %f\n",
	    cpu_time_on_vp_it/((float)(E->monitor.solution_cycles-E->control.restart)));
    fprintf(E->fp,"Initialization overhead = %f\n",initial_time);
    fprintf(E->fp,"Average cpu time taken for velocity step = %f\n",
	    cpu_time_on_vp_it/((float)(E->monitor.solution_cycles-E->control.restart)));
  }
  citcom_finalize(E, 0);
  return(0);

}
