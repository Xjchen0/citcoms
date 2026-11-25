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
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#ifdef USE_GGRD
#include "ggrd_handling.h"
#endif

/*=======================================================================
  Calculate ages (MY) for opening input files -> material, ages, velocities
  Open these files, read in results, and average if necessary
=========================================================================*/

void regional_read_input_files_for_timesteps(struct All_variables *E, int action, int output)
{
    /* modern c++: float find_age_in_MY() prototype in header */

    FILE *fp1, *fp2;
    float age, newage1, newage2;
    char output_file1[1024],output_file2[1024];
    float *TB1, *TB2, *VB1[5],*VB2[5], inputage1, inputage2;
    float *ST1, *ST2, *SS1, *SS2;
    int nox,noz,noy,nnn,nox1,noz1,noy1;
    int i,ii,ll,mm,j,k,n,nodeg,nodel,node;
    int intage, pos_age;
    int nodea;
    int nn, el;

    const int dims=E->mesh.nsd;

    int elx,ely,elz,elg,emax;
    float *VIP1,*VIP2;
    int *LL1, *LL2;

    int llayer;
    /* modern c++: int layers() prototype in header */

    /*if( E->parallel.me == 0)  
        fprintf(stderr, "\nINSIDE regional_read_input_files_for_timesteps   action=%d\n",action); */

    nox=E->mesh.nox;
    noy=E->mesh.noy;
    noz=E->mesh.noz;

    nox1=E->lmesh.nox;
    noz1=E->lmesh.noz;
    noy1=E->lmesh.noy;


    elx=E->lmesh.elx;
    elz=E->lmesh.elz;
    ely=E->lmesh.ely;

    emax=E->mesh.elx*E->mesh.elz*E->mesh.ely;

    age=find_age_in_MY(E);

    const int exact_age_action = (action == 1 || action == 2 || action == 5 || action == 6 || action == 7);
    if (age < 0.0f || (exact_age_action && age <= 1.0e-5f)) { /* age <= 0 -> use age=0 */
      intage = 0;
      newage2 = newage1 = 0.0f;
      pos_age = 0;
    }
    else {
      float rounded_age = floorf(age + 0.5f);
      if (exact_age_action && fabsf(age - rounded_age) < 1.0e-5f) {
        /* Exact integer age: read single file, avoiding crash from non-existent (max_age + 1) file. */
        intage = (int)rounded_age;
        newage2 = newage1 = rounded_age;
        pos_age = 0;
      } else {
        intage = (int)age;
        newage1 = 1.0f * intage;
        newage2 = 1.0f * intage + 1.0f;
        pos_age = 1;
      }
    }

    int need_vel_read = 1;
    if (action == 1) {
      need_vel_read = (E->parallel.me_loc[3] == E->parallel.nprocz - 1) ? 1 : 0;
      if (!need_vel_read && E->mesh.toplayerbc > 0) {
        for (int kz = 1; kz <= noz1; kz++) {
          if (E->SX[E->mesh.levmax][1][3][kz] >= E->mesh.toplayerbc_r) {
            need_vel_read = 1;
            break;
          }
        }
      } else if (!need_vel_read && E->mesh.toplayerbc < 0) {
        int gk = E->mesh.noz + E->mesh.toplayerbc;
        int lk = gk - E->lmesh.nzs + 1;
        if (lk >= 1 && lk <= noz1)
          need_vel_read = 1;
      }
    }

    switch (action) { /* set up files to open */
    case 1:  /* read velocity boundary conditions */
#ifdef USE_GGRD
      if(!E->control.ggrd.vtop_control){	/* regular input */
#endif
      if (!need_vel_read)
        break;
      snprintf(output_file1, sizeof(output_file1), "%s%0.0f", E->control.velocity_boundary_file, newage1);
      snprintf(output_file2, sizeof(output_file2), "%s%0.0f", E->control.velocity_boundary_file, newage2);
      fp1=fopen(output_file1,"r");
	if (fp1 == NULL) {
          fprintf(stderr, "(Problem_related #4) Cannot open %s\n", output_file1);
          if (E->fp) fprintf(E->fp,"(Problem_related #4) Cannot open %s\n",output_file1);
          parallel_process_termination();
	}
      if (pos_age) {
        fp2=fopen(output_file2,"r");
	 if (fp2 == NULL) {
          fprintf(stderr, "(Problem_related #5) Cannot open %s\n", output_file2);
          if (E->fp) fprintf(E->fp,"(Problem_related #5) Cannot open %s\n",output_file2);
          fclose(fp1);
          parallel_process_termination();
	 }
      }
      if((E->parallel.me==0) && (output==1))   {
         fprintf(E->fp,"Velocity: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
         fprintf(E->fp,"Velocity: File1 = %s\n",output_file1);
        if (pos_age)
           fprintf(E->fp,"Velocity: File2 = %s\n",output_file2);
        else
           fprintf(E->fp,"Velocity: File2 = No file inputted (exact/non-positive age)\n");
      }
#ifdef USE_GGRD
      }	
#endif
      break;

      case 2:  /* read ages for lithosphere temperature assimilation */
#ifdef USE_GGRD
      if(!E->control.ggrd.age_control){	/* regular input */
#endif
        sprintf(output_file1,"%s%0.0f",E->control.lith_age_file,newage1);
        sprintf(output_file2,"%s%0.0f",E->control.lith_age_file,newage2);
        fp1=fopen(output_file1,"r");
        if (fp1 == NULL) {
          fprintf(E->fp,"(Problem_related #6) Cannot open %s\n",output_file1);
          exit(8);
        }
        if (pos_age) {
          fp2=fopen(output_file2,"r");
          if (fp2 == NULL) {
            fprintf(E->fp,"(Problem_related #7) Cannot open %s\n",output_file2);            exit(8);
          }
        }
        if((E->parallel.me==0) && (output==1))   {
          fprintf(E->fp,"Age: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
          fprintf(E->fp,"Age: File1 = %s\n",output_file1);
          if (pos_age)
            fprintf(E->fp,"Age: File2 = %s\n",output_file2);
          else
            fprintf(E->fp,"Age: File2 = No file inputted (negative age)\n");
        }
#ifdef USE_GGRD
      }	
#endif
        break;

      case 3:  /* read element materials */
#ifdef USE_GGRD
	if(E->control.ggrd.mat_control == 0 ){
#endif
        sprintf(output_file1,"%s%0.0f.0",E->control.mat_file,newage1);
        sprintf(output_file2,"%s%0.0f.0",E->control.mat_file,newage2);
        fp1=fopen(output_file1,"r");
        if (fp1 == NULL) {
          fprintf(E->fp,"(Problem_related #8) Cannot open %s\n",output_file1);
          exit(8);
        }
        if (pos_age) {
          fp2=fopen(output_file2,"r");
          if (fp2 == NULL) {
            fprintf(E->fp,"(Problem_related #9) Cannot open %s\n",output_file2);
            exit(8);
          }
        }
        if((E->parallel.me==0) && (output==1))   {
          fprintf(E->fp,"Mat: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
          fprintf(E->fp,"Mat: File1 = %s\n",output_file1);
          if (pos_age)
            fprintf(E->fp,"Mat: File2 = %s\n",output_file2);
          else
            fprintf(E->fp,"Mat: File2 = No file inputted (negative age)\n");
        }

#ifdef USE_GGRD
	}
#endif
	break;

	/* mode 4 is rayleigh control for GGRD, see below */

      case 5:  /* read temperature boundary conditions, top surface */
        sprintf(output_file1,"%s%0.0f",E->control.temperature_boundary_file,newage1);
        sprintf(output_file2,"%s%0.0f",E->control.temperature_boundary_file,newage2);
        fp1=fopen(output_file1,"r");
	  if (fp1 == NULL) {
            fprintf(E->fp,"(Problem_related #10) Cannot open %s\n",output_file1);
            exit(8);
	  }
        if (pos_age) {
          fp2=fopen(output_file2,"r");
	   if (fp2 == NULL) {
            fprintf(E->fp,"(Problem_related #11) Cannot open %s\n",output_file2);
            exit(8);
	   }
        }
        if((E->parallel.me==0) && (output==1))   {
           fprintf(E->fp,"Surface Temperature: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
           fprintf(E->fp,"Surface Temperature: File1 = %s\n",output_file1);
          if (pos_age)
             fprintf(E->fp,"Surface Temperature: File2 = %s\n",output_file2);
          else
             fprintf(E->fp,"Surface Temperature: File2 = No file inputted (negative age)\n");
        }
      break;

      case 6:  /* EarthByte: read temperature and stencil for slab assimilation */
        snprintf(output_file1, sizeof(output_file1), "%s%0.0f", E->control.slab_assim_file, newage1);
        snprintf(output_file2, sizeof(output_file2), "%s%0.0f", E->control.slab_assim_file, newage2);
        fp1 = fopen(output_file1, "r");
        if (fp1 == NULL) {
          fprintf(stderr, "(Problem_related #12) Cannot open %s\n", output_file1);
          if (E->fp) fprintf(E->fp, "(Problem_related #12) Cannot open %s\n", output_file1);
          parallel_process_termination();
        }
        if (pos_age) {
          fp2 = fopen(output_file2, "r");
          if (fp2 == NULL) {
            fprintf(stderr, "(Problem_related #13) Cannot open %s\n", output_file2);
            if (E->fp) fprintf(E->fp, "(Problem_related #13) Cannot open %s\n", output_file2);
            fclose(fp1);
            parallel_process_termination();
          }
        }
        if ((E->parallel.me == 0) && (output == 1)) {
          fprintf(E->fp, "Slab Assimilation: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",
                  E->control.start_age, E->monitor.elapsed_time, age);
          fprintf(E->fp, "Slab Assimilation: File1 = %s\n", output_file1);
          if (pos_age)
            fprintf(E->fp, "Slab Assimilation: File2 = %s\n", output_file2);
          else
            fprintf(E->fp, "Slab Assimilation: File2 = No file inputted (exact/non-positive age)\n");
        }
        break;

      case 7:  /* EarthByte: read internal velocity and stencil */
        snprintf(output_file1, sizeof(output_file1), "%s%0.0f", E->control.velocity_internal_file, newage1);
        snprintf(output_file2, sizeof(output_file2), "%s%0.0f", E->control.velocity_internal_file, newage2);
        fp1 = fopen(output_file1, "r");
        if (fp1 == NULL) {
          fprintf(stderr, "(Problem_related #14) Cannot open %s\n", output_file1);
          if (E->fp) fprintf(E->fp, "(Problem_related #14) Cannot open %s\n", output_file1);
          parallel_process_termination();
        }
        if (pos_age) {
          fp2 = fopen(output_file2, "r");
          if (fp2 == NULL) {
            fprintf(stderr, "(Problem_related #15) Cannot open %s\n", output_file2);
            if (E->fp) fprintf(E->fp, "(Problem_related #15) Cannot open %s\n", output_file2);
            fclose(fp1);
            parallel_process_termination();
          }
        }
        if ((E->parallel.me == 0) && (output == 1)) {
          fprintf(E->fp, "Internal Velocity: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",
                  E->control.start_age, E->monitor.elapsed_time, age);
          fprintf(E->fp, "Internal Velocity: File1 = %s\n", output_file1);
          if (pos_age)
            fprintf(E->fp, "Internal Velocity: File2 = %s\n", output_file2);
          else
            fprintf(E->fp, "Internal Velocity: File2 = No file inputted (exact/non-positive age)\n");
        }
        break;

    } /* end switch */



    switch (action) { /* Read the contents of files and average */

    case 1:  /* velocity boundary conditions */
#ifdef USE_GGRD
      if(!E->control.ggrd.vtop_control){ /* grd control is called from boundary condition file */
#endif
      if (!need_vel_read)
        break;
      nnn=nox*noy;
      for(i=1;i<=2;i++)  {
        VB1[i]=(float*) malloc ((nnn+1)*sizeof(float));
        VB2[i]=(float*) malloc ((nnn+1)*sizeof(float));
      }
      for(i=1;i<=nnn;i++)   {
         if(fscanf(fp1,"%f %f",&(VB1[1][i]),&(VB1[2][i])) != 2) {
           fprintf(stderr,"Error while reading file '%s'\n",output_file1);
           fclose(fp1);
           if (pos_age) fclose(fp2);
           parallel_process_termination();
         }
         VB1[1][i]=E->data.timedir*VB1[1][i];
         VB1[2][i]=E->data.timedir*VB1[2][i];
         if (pos_age) {
             if(fscanf(fp2,"%f %f",&(VB2[1][i]),&(VB2[2][i])) != 2) {
                 fprintf(stderr,"Error while reading file '%s'\n",output_file2);
                 fclose(fp1);
                 fclose(fp2);
                 parallel_process_termination();
             }
             VB2[1][i]=E->data.timedir*VB2[1][i];
             VB2[2][i]=E->data.timedir*VB2[2][i];
         }
      }
      fclose(fp1);
      if (pos_age) fclose(fp2);

      {
        const int is_top_rank = (E->parallel.me_loc[3] == E->parallel.nprocz - 1);
        const int internal_lk = (E->mesh.toplayerbc < 0) ? (E->mesh.noz + E->mesh.toplayerbc - E->lmesh.nzs + 1) : -1;
        for (int kz = 1; kz <= noz1; kz++) {
          const int is_top_layer = (is_top_rank && kz == noz1);
          int assign_layer = is_top_layer;
          if (!assign_layer && E->mesh.toplayerbc > 0 && E->SX[E->mesh.levmax][1][3][kz] >= E->mesh.toplayerbc_r)
            assign_layer = 1;
          else if (!assign_layer && E->mesh.toplayerbc < 0 && kz == internal_lk)
            assign_layer = 1;
          if (!assign_layer)
            continue;
          for(k=1;k<=noy1;k++)
             for(i=1;i<=nox1;i++)    {
                nodeg = E->lmesh.nxs+i-1 + (E->lmesh.nys+k-2)*nox;
                nodel = (k-1)*nox1*noz1 + (i-1)*noz1 + kz;
		if (pos_age) {
                    E->sphere.cap[1].VB[1][nodel] = (VB1[1][nodeg] + (VB2[1][nodeg]-VB1[1][nodeg])/(newage2-newage1)*(age-newage1))*E->data.scalev;
                    E->sphere.cap[1].VB[2][nodel] = (VB1[2][nodeg] + (VB2[2][nodeg]-VB1[2][nodeg])/(newage2-newage1)*(age-newage1))*E->data.scalev;
		}
		else {
                    E->sphere.cap[1].VB[1][nodel] = VB1[1][nodeg]*E->data.scalev;
                    E->sphere.cap[1].VB[2][nodel] = VB1[2][nodeg]*E->data.scalev;
		}
		if (is_top_layer)
                    E->sphere.cap[1].VB[3][nodel] = 0.0;
             }
        }
      }
      for(i=1;i<=2;i++) {
          free ((void *) VB1[i]);
          free ((void *) VB2[i]);
      }


#ifdef USE_GGRD
      } /* end of branch if allowing for ggrd handling */
#endif
      break;

      case 2:  /* ages for lithosphere temperature assimilation */
        for(i=1;i<=noy;i++)
          for(j=1;j<=nox;j++) {
            node=j+(i-1)*nox;
            if(fscanf(fp1,"%f",&inputage1) != 1) {
              fprintf(stderr,"Error while reading file '%s'\n",output_file1);
              fclose(fp1);
              if (pos_age) fclose(fp2);
              parallel_process_termination();
            }
            if (pos_age) { /* positive ages - we must interpolate */
              if(fscanf(fp2,"%f",&inputage2) != 1) {
                fprintf(stderr,"Error while reading file '%s'\n",output_file2);
                fclose(fp1);
                fclose(fp2);
                parallel_process_termination();
              }
              /* EarthByte: preserve stencil masking if either endpoint is masked */
              if (inputage1 <= E->control.lith_age_stencil_value ||
                  inputage2 <= E->control.lith_age_stencil_value) {
                inputage1 = min(inputage1, inputage2);
                inputage2 = inputage1;
              }
              E->age_t[node] = (inputage1 + (inputage2-inputage1)/(newage2-newage1)*(age-newage1))/E->data.scalet;
            }
            else { /* exact or non-positive age - convert Myr to dimensionless time */
              E->age_t[node] = inputage1 / E->data.scalet;
            }
          }
        fclose(fp1);
        if (pos_age) fclose(fp2);
        break;

      case 3:  /* read element materials */
#ifdef USE_GGRD
	if(E->control.ggrd.mat_control != 0){
	  ggrd_read_mat_from_file(E, 1);
	}else{
#endif

        VIP1 = (float*) malloc ((emax+1)*sizeof(float));
        VIP2 = (float*) malloc ((emax+1)*sizeof(float));
        LL1 = (int*) malloc ((emax+1)*sizeof(int));
        LL2 = (int*) malloc ((emax+1)*sizeof(int));

          for (el=1; el<=elx*ely*elz; el++)  {
            nodea = E->ien[1][el].node[2];
            llayer = layers(E,1,nodea);
            if (llayer)  { /* for layers:1-lithosphere,2-upper, 3-trans, and 4-lower mantle */
              E->mat[1][el] = llayer;
            }
          }
          for(i=1;i<=emax;i++)  {
              if(fscanf(fp1,"%d %d %f", &nn,&(LL1[i]),&(VIP1[i])) != 3) {
                  fprintf(stderr,"Error while reading file '%s'\n",output_file1);
                  exit(8);
              }
              if (pos_age) {
                  if(fscanf(fp2,"%d %d %f", &nn,&(LL2[i]),&(VIP2[i])) != 3) {
                      fprintf(stderr,"Error while reading file '%s'\n",output_file2);
                      exit(8);
                  }
              }
          }

          fclose(fp1);
          if (pos_age) fclose(fp2);

          for (k=1;k<=ely;k++)   {
            for (i=1;i<=elx;i++)   {
              for (j=1;j<=elz;j++)  {
                el = j + (i-1)*E->lmesh.elz + (k-1)*E->lmesh.elz*E->lmesh.elx;
                elg = E->lmesh.ezs+j + (E->lmesh.exs+i-1)*E->mesh.elz + (E->lmesh.eys+k-1)*E->mesh.elz*E->mesh.elx;
                if (pos_age) { /* positive ages - we must interpolate */
                    E->VIP[1][el] = VIP1[elg]+(VIP2[elg]-VIP1[elg])/(newage2-newage1)*(age-newage1);
                }
                else { /* negative ages - don't do the interpolation */
                    E->VIP[1][el] = VIP1[elg];
                }

                /* E->mat[1][el] = LL1[elg]; */ /*use the mat numbers base on radius*/

              }     /* end for j  */
            }     /*  end for i */
          }     /*  end for k  */

         free ((void *) VIP1);
         free ((void *) VIP2);
         free ((void *) LL1);
         free ((void *) LL2);
#ifdef USE_GGRD
	} /* end of branch if allowing for ggrd handling */
#endif
	break;
    case 4:			/* material control */
#ifdef USE_GGRD
      if(E->control.ggrd.ray_control)
	ggrd_read_ray_from_file(E, 1);
#else
      myerror(E,"input_from_files: mode 4 only for GGRD");
#endif
      break;

    case 5:  /* read temperature boundary conditions, top surface */
      nnn=nox*noy;
      TB1=(float*) malloc ((nnn+1)*sizeof(float));
      TB2=(float*) malloc ((nnn+1)*sizeof(float));

      for(i=1;i<=nnn;i++)   {
        if(fscanf(fp1,"%f",&(TB1[i])) != 1) {
          fprintf(stderr,"Error while reading file '%s'\n",output_file1);
          exit(8);
        }
         /* if( E->parallel.me == 0)  
        fprintf(stderr, "\nINSIDE regional_read_input_files_for_timesteps TB1=%f %d\n",TB1[i],i); */
         if (pos_age) {
             if(fscanf(fp2,"%f",&(TB2[i])) != 1) {
                 fprintf(stderr,"Error while reading file '%s'\n",output_file2);
                 exit(8);
             }
         }
      }
      fclose(fp1);
      if (pos_age) fclose(fp2);

      if(E->parallel.me_loc[3]==E->parallel.nprocz-1 )  {
          for(k=1;k<=noy1;k++)
             for(i=1;i<=nox1;i++)    {
                nodeg = E->lmesh.nxs+i-1 + (E->lmesh.nys+k-2)*nox;
                nodel = (k-1)*nox1*noz1 + (i-1)*noz1+noz1;
		if (pos_age) { /* positive ages - we must interpolate */
                    E->sphere.cap[1].TB[1][nodel] = (TB1[nodeg] + (TB2[nodeg]-TB1[nodeg])/(newage2-newage1)*(age-newage1));
                    E->sphere.cap[1].TB[2][nodel] = (TB1[nodeg] + (TB2[nodeg]-TB1[nodeg])/(newage2-newage1)*(age-newage1));
                    E->sphere.cap[1].TB[3][nodel] = (TB1[nodeg] + (TB2[nodeg]-TB1[nodeg])/(newage2-newage1)*(age-newage1));
		}
		else { /* negative ages - don't do the interpolation */
                    E->sphere.cap[1].TB[1][nodel] = TB1[nodeg];
                    E->sphere.cap[1].TB[2][nodel] = TB1[nodeg];
                    E->sphere.cap[1].TB[3][nodel] = TB1[nodeg];
		}
             }
      }   /* end of E->parallel.me_loc[3]==E->parallel.nprocz-1   */
      free ((void *) TB1);
      free ((void *) TB2);

      break;

      case 6:  /* EarthByte: read temperature and stencil for slab assimilation */
        nnn = nox * noy * noz;
        ST1 = (float*) malloc((nnn + 1) * sizeof(float));
        ST2 = (float*) malloc((nnn + 1) * sizeof(float));
        SS1 = (float*) malloc((nnn + 1) * sizeof(float));
        SS2 = (float*) malloc((nnn + 1) * sizeof(float));

        for (i = 1; i <= nnn; i++) {
          if (fscanf(fp1, "%g %g", &(ST1[i]), &(SS1[i])) != 2) {
            fprintf(stderr, "Error while reading file '%s'\n", output_file1);
            fclose(fp1);
            if (pos_age) fclose(fp2);
            parallel_process_termination();
          }
          if (pos_age) {
            if (fscanf(fp2, "%g %g", &(ST2[i]), &(SS2[i])) != 2) {
              fprintf(stderr, "Error while reading file '%s'\n", output_file2);
              fclose(fp1);
              fclose(fp2);
              parallel_process_termination();
            }
          }
        }
        fclose(fp1);
        if (pos_age) fclose(fp2);

        for (j = 1; j <= noy1; j++)
          for (i = 1; i <= nox1; i++)
            for (k = 1; k <= noz1; k++) {
              nodel = k + (i - 1) * noz1 + (j - 1) * nox1 * noz1;
              nodeg = (E->lmesh.nzs + k - 1) + (E->lmesh.nxs + i - 2) * noz + (E->lmesh.nys + j - 2) * noz * nox;

              if (pos_age) {
                E->sphere.cap[1].slab_temp[nodel] = ST1[nodeg] + (ST2[nodeg] - ST1[nodeg]) / (newage2 - newage1) * (age - newage1);
                E->sphere.cap[1].slab_sten[nodel] = SS1[nodeg] + (SS2[nodeg] - SS1[nodeg]) / (newage2 - newage1) * (age - newage1);
              } else {
                E->sphere.cap[1].slab_temp[nodel] = ST1[nodeg];
                E->sphere.cap[1].slab_sten[nodel] = SS1[nodeg];
              }
            }

        free((void *) ST1);
        free((void *) ST2);
        free((void *) SS1);
        free((void *) SS2);
        break;

      case 7:  /* EarthByte: read internal velocity and stencil */
        nnn = nox * noy * noz;
        for (i = 1; i <= dims + 1; i++) {
          VB1[i] = (float*) malloc((nnn + 1) * sizeof(float));
          VB2[i] = (float*) malloc((nnn + 1) * sizeof(float));
        }

        for (i = 1; i <= nnn; i++) {
          if (fscanf(fp1, "%f %f %f %f", &(VB1[1][i]), &(VB1[2][i]), &(VB1[3][i]), &(VB1[4][i])) != 4) {
            fprintf(stderr, "Error while reading file '%s'\n", output_file1);
            fclose(fp1);
            if (pos_age) fclose(fp2);
            parallel_process_termination();
          }
          VB1[1][i] *= E->data.timedir;
          VB1[2][i] *= E->data.timedir;
          VB1[3][i] *= E->data.timedir;
          if (pos_age) {
            if (fscanf(fp2, "%f %f %f %f", &(VB2[1][i]), &(VB2[2][i]), &(VB2[3][i]), &(VB2[4][i])) != 4) {
              fprintf(stderr, "Error while reading file '%s'\n", output_file2);
              fclose(fp1);
              fclose(fp2);
              parallel_process_termination();
            }
            VB2[1][i] *= E->data.timedir;
            VB2[2][i] *= E->data.timedir;
            VB2[3][i] *= E->data.timedir;
          }
        }
        fclose(fp1);
        if (pos_age) fclose(fp2);

        for (j = 1; j <= noy1; j++)
          for (i = 1; i <= nox1; i++)
            for (k = 1; k <= noz1; k++) {
              nodel = k + (i - 1) * noz1 + (j - 1) * nox1 * noz1;
              nodeg = (E->lmesh.nzs + k - 1) + (E->lmesh.nxs + i - 2) * noz + (E->lmesh.nys + j - 2) * noz * nox;

              if ((int)VB1[4][nodeg] != 2) {
                if (pos_age && (fabs(age - newage1) > fabs(age - newage2))) {
                  E->sphere.cap[1].VB[1][nodel] = VB2[1][nodeg] * E->data.scalev;
                  E->sphere.cap[1].VB[2][nodel] = VB2[2][nodeg] * E->data.scalev;
                  E->sphere.cap[1].VB[3][nodel] = VB2[3][nodeg] * E->data.scalev;
                  E->sphere.cap[1].slab_sten2[nodel] = (int)VB2[4][nodeg];
                } else {
                  E->sphere.cap[1].VB[1][nodel] = VB1[1][nodeg] * E->data.scalev;
                  E->sphere.cap[1].VB[2][nodel] = VB1[2][nodeg] * E->data.scalev;
                  E->sphere.cap[1].VB[3][nodel] = VB1[3][nodeg] * E->data.scalev;
                  E->sphere.cap[1].slab_sten2[nodel] = (int)VB1[4][nodeg];
                }
              } else {
                E->sphere.cap[1].slab_sten2[nodel] = (int)VB1[4][nodeg];
              }
            }

        for (i = 1; i <= dims + 1; i++) {
          free((void *) VB1[i]);
          free((void *) VB2[i]);
        }
        break;

    } /* end switch */

   return;
}
