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

/* Routine to process the output of the finite element cycles
   and to turn them into a coherent suite of files  */


#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "output.h"
#include "string.h"

#if defined(USE_GZDIR) || defined(USE_ZLIB)
#include "zlib.h"
#endif

#define CHUNK 16384

static void write_binary_array(int nn, float* array, FILE * f);
static void write_ascii_array(int nn, int perLine, float *array, FILE *fp);

static void vts_file_header(struct All_variables *E, FILE *fp)
{

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <StructuredGrid WholeExtent=\"%s\">\n"
        "    <Piece Extent=\"%s\">\n";

    char extent[64], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
             E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz,
             E->lmesh.exs, E->lmesh.exs + E->lmesh.elx,
             E->lmesh.eys, E->lmesh.eys + E->lmesh.ely);

    snprintf(header, 1024, format, extent, extent);

    fputs(header, fp);

    return;
}


static void vts_file_trailer(struct All_variables *E, FILE *fp)
{
    const char trailer[] =
        "    </Piece>\n"
        "  </StructuredGrid>\n"
        "</VTKFile>\n";

    fputs(trailer, fp);

    return;
}


static void vtk_point_data_header(struct All_variables *E, FILE *fp)
{
    fputs("      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n", fp);
    return;
}


static void vtk_point_data_trailer(struct All_variables *E, FILE *fp)
{
    fputs("      </PointData>\n", fp);
    return;
}


static void vtk_cell_data_header(struct All_variables *E, FILE *fp)
{
    fputs("      <CellData>\n", fp);
    return;
}


static void vtk_cell_data_trailer(struct All_variables *E, FILE *fp)
{
    fputs("      </CellData>\n", fp);
    return;
}


static void vtk_output_temp(struct All_variables *E, FILE *fp)
{
    int i;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floattemp = (float*)malloc(nodes*sizeof(float));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\">\n", E->output.vtk_format);

    for(i=0;i < nodes;i++)
        floattemp[i] =  (float) *(E->T[1]+i+1);

    if (strcmp(E->output.vtk_format,"binary") == 0) {
        write_binary_array(nodes,floattemp,fp);
    } else {
        write_ascii_array(nodes,1,floattemp,fp);
    }
    fputs("        </DataArray>\n", fp);
    free(floattemp);
    return;
}


static void vtk_output_velo(struct All_variables *E, FILE *fp)
{
    int i, j;
    int nodes=E->sphere.caps_per_proc*E->lmesh.nno;
    double sint, sinf, cost, cosf;
    float *V[4];
    const int lev = E->mesh.levmax;
    float* floatvel = (float*)malloc(nodes*3*sizeof(float));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        V[1] = E->sphere.cap[j].V[1];
        V[2] = E->sphere.cap[j].V[2];
        V[3] = E->sphere.cap[j].V[3];

        for(i=1; i<=E->lmesh.nno; i++) {
            sint = E->SinCos[lev][j][0][i];
            sinf = E->SinCos[lev][j][1][i];
            cost = E->SinCos[lev][j][2][i];
            cosf = E->SinCos[lev][j][3][i];

            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+0] = (float)(V[1][i]*cost*cosf - V[2][i]*sinf + V[3][i]*sint*cosf);
            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+1] = (float)(V[1][i]*cost*sinf + V[2][i]*cosf + V[3][i]*sint*sinf);
            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+2] = (float)(-V[1][i]*sint + V[3][i]*cost);
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes*3,floatvel,fp);
    else
        write_ascii_array(nodes*3,3,floatvel,fp);
    fputs("        </DataArray>\n", fp);

    free(floatvel);
    return;
}


static void vtk_output_visc(struct All_variables *E, FILE *fp)
{
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    int lev = E->mesh.levmax;

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0) {
            write_binary_array(nodes,&E->VI[lev][1][1],fp);
        } else {
            write_ascii_array(nodes,1,&E->VI[lev][1][1],fp);
        }

    fputs("        </DataArray>\n", fp);
    return;
}


static void vtk_output_coord(struct All_variables *E, FILE *fp)
{
    /* Output Cartesian coordinates as most VTK visualization softwares
       assume it. */
    int i, j;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floatpos = (float*)malloc(nodes*3*sizeof(float));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++){
                floatpos[((j-1)*E->lmesh.nno+i-1)*3] = (float)(E->x[j][1][i]);
	        floatpos[((j-1)*E->lmesh.nno+i-1)*3+1]=(float)(E->x[j][2][i]);
	        floatpos[((j-1)*E->lmesh.nno+i-1)*3+2]=(float)(E->x[j][3][i]);
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes*3,floatpos,fp);
    else
        write_ascii_array(nodes*3,3,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(floatpos);
    return;
}

static void vtk_output_stress(struct All_variables *E, FILE *fp)
{
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
   /* for stress computation */
    /* modern c++: void allocate_STD_mem() prototype in header */
    /* modern c++: void compute_nodal_stress() prototype in header */
    /* modern c++: void free_STD_mem() prototype in header */
    float *SXX[NCS],*SYY[NCS],*SXY[NCS],*SXZ[NCS],*SZY[NCS],*SZZ[NCS];
    float *divv[NCS],*vorv[NCS];

    /* those are sorted like stt spp srr stp str srp  */
    allocate_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    compute_nodal_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    free_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\">\n", E->output.vtk_format);

    if (strcmp(E->output.vtk_format, "binary") == 0) {
        write_binary_array(nodes*6,&E->gstress[1][1],fp);
    } else {
        write_ascii_array(nodes*6,6,&E->gstress[1][1],fp);
    }

    fputs("        </DataArray>\n", fp);
    return;
}

static void vtk_output_comp_nd(struct All_variables *E, FILE *fp)
{
    int i, j, k;
    char name[255];
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floatcompo = (float*)malloc(nodes*sizeof(float));

    for(k=0;k<E->composition.ncomp;k++) {
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"composition%d\" format=\"%s\">\n", k+1, E->output.vtk_format);

        for(j=1; j<=E->sphere.caps_per_proc; j++) {
            for(i=1; i<=E->lmesh.nno; i++) {
                floatcompo[(j-1)*E->lmesh.nno+i-1] = (float) (E->composition.comp_node[j][k][i]);
	    }
        }

        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(nodes,floatcompo,fp);
        else
            write_ascii_array(nodes,1,floatcompo,fp);
        fputs("        </DataArray>\n", fp);
    }
    free(floatcompo);
    return;
}


static void vtk_output_surf(struct All_variables *E,  FILE *fp, int cycles)
{
    int i, j, k;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    char output_file[255];
    float* floattopo = (float*)malloc(nodes*sizeof(float));

    if((E->output.write_q_files == 0) || (cycles == 0) ||
      (cycles % E->output.write_q_files)!=0)
        heat_flux(E);
  /* else, the heat flux will have been computed already */

    if(E->control.use_cbf_topo){
        get_CBF_topo(E,E->slice.tpg,E->slice.tpgb);
    }
    else{
        get_STD_topo(E,E->slice.tpg,E->slice.tpgb,E->slice.divg,E->slice.vort,cycles);
    }

    fprintf(fp,"        <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1;j<=E->sphere.caps_per_proc;j++){
        for(i=1;i<=E->lmesh.nsf;i++){
            for(k=1;k<=E->lmesh.noz;k++){
                floattopo[(j-1)*E->lmesh.nno + (i-1)*E->lmesh.noz + k-1] = 0.0;
            }

            if (E->parallel.me_loc[3]==E->parallel.nprocz-1) {

                /* choose either STD topo or pseudo-free-surf topo */
                if(E->control.pseudo_free_surf)
                floattopo[(j-1)*E->lmesh.nno + i*E->lmesh.noz-1] = E->slice.freesurf[j][i];
                else
                floattopo[(j-1)*E->lmesh.nno + i*E->lmesh.noz-1] = E->slice.tpg[j][i];

            }
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes,floattopo,fp);
    else
        write_ascii_array(nodes,1,floattopo,fp);

    fputs("        </DataArray>\n", fp);
    free(floattopo);
    return;
}


static void write_vtm(struct All_variables *E, int cycles)
{
    FILE *fp;
    char vtm_file[255];
    int n;

    const char header[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <vtkMultiBlockDataSet>\n";

    snprintf(vtm_file, 255, "%s.%d.vtm",
             E->control.data_file, cycles);
    fp = output_open(vtm_file, "w");
    fputs(header, fp);

    for(n=0; n<E->parallel.nproc; n++) {
        fprintf(fp, "    <DataSet index=\"%d\" file=\"%s.proc%d.%d.vts\"/>\n",
                n, E->control.data_prefix, n, cycles);
    }
    fputs("  </vtkMultiBlockDataSet>\n",fp);
    fputs("</VTKFile>",fp);

    fclose(fp);
}

static void write_vtm_caps(struct All_variables *E, int cycles)
{
    FILE *fp;
    char vtm_file[255];
    int c;

    const char header[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <vtkMultiBlockDataSet>\n";

    snprintf(vtm_file, 255, "%s.%d.vtm",
             E->control.data_file, cycles);
    fp = output_open(vtm_file, "w");
    fputs(header, fp);

    for(c = 0; c < E->sphere.caps; c++) {
        fprintf(fp, "    <DataSet index=\"%d\" file=\"%s.cap%d.%d.vts\"/>\n",
                c, E->control.data_prefix, c, cycles);
    }
    fputs("  </vtkMultiBlockDataSet>\n", fp);
    fputs("</VTKFile>", fp);

    fclose(fp);
}

static void write_visit_caps(struct All_variables *E, int cycles)
{
    FILE *fp;
    char visit_file[255];
    int c;

    const char header[] = "!NBLOCKS %d\n";

    snprintf(visit_file, 255, "%s.%d.visit",
             E->control.data_file, cycles);
    fp = output_open(visit_file, "w");
    fprintf(fp, header, E->sphere.caps);

    for(c = 0; c < E->sphere.caps; c++) {
        fprintf(fp, "%s.cap%d.%d.vts\n",
                E->control.data_prefix, c, cycles);
    }
    fclose(fp);
}

static void write_visit(struct All_variables *E, int cycles)
{
    FILE *fp;
    char visit_file[255];
    int n;

    const char header[] = "!NBLOCKS %d\n";

    snprintf(visit_file, 255, "%s.%d.visit",
             E->control.data_file, cycles);
    fp = output_open(visit_file, "w");
    fprintf(fp, header, E->parallel.nproc);

    for(n=0; n<E->parallel.nproc; n++) {
        fprintf(fp, "%s.proc%d.%d.vts\n",
                E->control.data_prefix, n, cycles);
    }
    fclose(fp);
}

static void write_pvts(struct All_variables *E, int cycles)
{
    FILE *fp;
    char pvts_file[255];
    int i,j,k;
    snprintf(pvts_file, 255, "%s.%d.pvts",
             E->control.data_file,cycles);
    fp = output_open(pvts_file, "w");

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <PStructuredGrid WholeExtent=\"%s\" GhostLevel=\"#\">\n"
        "    <PPointData Scalars=\"temperature\" Vectors=\"velocity\">\n"
        "      <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\"/>\n"
        "      <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\"/>\n"
        "      <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\"/>\n";

    char extent[64], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
        E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz*E->parallel.nprocz,
        E->lmesh.exs, E->lmesh.exs + E->lmesh.elx*E->parallel.nprocx,
        E->lmesh.eys, E->lmesh.eys + E->lmesh.ely*E->parallel.nprocy);

    snprintf(header, 1024, format, extent, E->output.vtk_format,
             E->output.vtk_format, E->output.vtk_format);
    fputs(header, fp);

    if (E->output.stress){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\"/>\n", E->output.vtk_format);
    }
    if (E->output.comp_nd && E->composition.on){
        for (int c = 0; c < E->composition.ncomp; c++) {
            fprintf(fp,"      <DataArray type=\"Float32\" Name=\"composition%d\" format=\"%s\"/>\n", c + 1, E->output.vtk_format);
        }
    }
    if (E->output.surf){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\"/>\n", E->output.vtk_format);
    }

    fputs("    </PPointData>\n \n"
    "    <PCellData>\n"
    "    </PCellData>\n \n"
    "    <PPoints>\n"
    "      <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"binary\" />\n"
    "    </PPoints>\n", fp);

    for(i=0; i < E->parallel.nprocy;i++){
        for(j=0; j < E->parallel.nprocx;j++){
            for(k=0; k < E->parallel.nprocz;k++){
                fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s.proc%d.%d.vts\"/>\n",
                    (k%E->parallel.nprocz)*E->lmesh.elz,
                    (k%E->parallel.nprocz+1)*E->lmesh.elz,
                    (j%E->parallel.nprocx)*E->lmesh.elx, (j%E->parallel.nprocx+1)*E->lmesh.elx,
                    (i%E->parallel.nprocy)*E->lmesh.ely, (i%E->parallel.nprocy+1)*E->lmesh.ely,
                    E->control.data_prefix,
                    i*E->parallel.nprocx*E->parallel.nprocz+j*E->parallel.nprocz+k, cycles);
            }
        }
    }

    fputs("  </PStructuredGrid>\n",fp);
    fputs("</VTKFile>",fp);

    fclose(fp);
}


static void write_pvd(struct All_variables *E, int cycles)
{
    FILE *fp_hist;
    FILE *fp_pvd;
    char hist_file[512];
    char pvd_file[512];
    char line[512];
    char target_file[512];
    double time_val;
    int max_entries = 10000;
    int count = 0;
    int i;
    int *hist_cycles;
    double *hist_times;
    char (*hist_files)[512];

    if (strcmp(E->output.vtk_output_mode, "single") == 0) {
        snprintf(target_file, sizeof(target_file), "%s.%d.vtu", E->control.data_prefix, cycles);
    } else if (strcmp(E->output.vtk_output_mode, "cap") == 0) {
        if (E->sphere.caps == 12) {
            snprintf(target_file, sizeof(target_file), "%s.%d.vtm", E->control.data_prefix, cycles);
        } else {
            snprintf(target_file, sizeof(target_file), "%s.cap0.%d.vts", E->control.data_prefix, cycles);
        }
    } else {
        if (E->sphere.caps == 12) {
            snprintf(target_file, sizeof(target_file), "%s.%d.vtm", E->control.data_prefix, cycles);
        } else {
            snprintf(target_file, sizeof(target_file), "%s.%d.pvts", E->control.data_prefix, cycles);
        }
    }

    if (E->data.scalet > 0.0) {
        time_val = (double)(E->monitor.elapsed_time * E->data.scalet);
    } else {
        time_val = (double)E->monitor.elapsed_time;
    }

    snprintf(hist_file, sizeof(hist_file), "%s.vtk_history.txt", E->control.data_file);
    snprintf(pvd_file, sizeof(pvd_file), "%s.pvd", E->control.data_file);

    hist_cycles = (int *)malloc(max_entries * sizeof(int));
    hist_times = (double *)malloc(max_entries * sizeof(double));
    hist_files = (char (*)[512])malloc(max_entries * 512);

    if (!hist_cycles || !hist_times || !hist_files) {
        if (hist_cycles) free(hist_cycles);
        if (hist_times) free(hist_times);
        if (hist_files) free(hist_files);
        return;
    }

    /* Read existing history and filter entries >= cycles (restart protection) */
    fp_hist = fopen(hist_file, "r");
    if (fp_hist) {
        while (fgets(line, sizeof(line), fp_hist) && count < max_entries - 1) {
            int c;
            double t;
            char f[512];
            if (sscanf(line, "%d %lf %511s", &c, &t, f) == 3) {
                if (c < cycles) {
                    hist_cycles[count] = c;
                    hist_times[count] = t;
                    strncpy(hist_files[count], f, 511);
                    hist_files[count][511] = '\0';
                    count++;
                }
            }
        }
        fclose(fp_hist);
    }

    /* Add current step */
    if (count < max_entries) {
        hist_cycles[count] = cycles;
        hist_times[count] = time_val;
        strncpy(hist_files[count], target_file, 511);
        hist_files[count][511] = '\0';
        count++;
    }

    /* Rewrite history file */
    fp_hist = fopen(hist_file, "w");
    if (fp_hist) {
        for (i = 0; i < count; i++) {
            fprintf(fp_hist, "%d %.8e %s\n", hist_cycles[i], hist_times[i], hist_files[i]);
        }
        fclose(fp_hist);
    }

    /* Write PVD collection file */
    fp_pvd = fopen(pvd_file, "w");
    if (fp_pvd) {
        fprintf(fp_pvd, "<?xml version=\"1.0\"?>\n");
        fprintf(fp_pvd, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(fp_pvd, "  <Collection>\n");
        for (i = 0; i < count; i++) {
            fprintf(fp_pvd, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\" file=\"%s\"/>\n",
                    hist_times[i], hist_files[i]);
        }
        fprintf(fp_pvd, "  </Collection>\n");
        fprintf(fp_pvd, "</VTKFile>\n");
        fclose(fp_pvd);
    }

    free(hist_cycles);
    free(hist_times);
    free(hist_files);
}

static void write_ascii_array(int nn, int perLine, float *array, FILE *fp)
{
    int i;

    switch (perLine) {
    case 1:
        for(i=0; i<nn; i++)
            fprintf(fp, "%.4e\n", array[i]);
        break;
    case 3:
        for(i=0; i < nn/3; i++)
            fprintf(fp,"%.4e %.4e %.4e\n",array[3*i],array[3*i+1],array[3*i+2]);
        break;
    case 6:
        for(i=0; i < nn/6; i++)
            fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e\n",
                    array[6*i],array[6*i+1],array[6*i+2],
                    array[6*i+3],array[6*i+4],array[6*i+5]);
        break;
    }
    return;
}

static void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray)
{
    /* simple float to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union FloatToUnsignedChars
        {
            float input;
            unsigned char output[4];
        } floattransform;

    for (i=0; i<nn; i++){
        floattransform.input=floatarray[i];
        chararray[4*i]=floattransform.output[0];
        chararray[4*i+1]=floattransform.output[1];
        chararray[4*i+2]=floattransform.output[2];
        chararray[4*i+3]=floattransform.output[3];
    }
    return;
}

static void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray)
{
    /* simple int - to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union IntToUnsignedChars
        {
            int input;
            unsigned char output[4];
        } inttransform;

    for (i=0; i<nn; i++){
        inttransform.input=intarray[i];
        chararray[4*i]=inttransform.output[0];
        chararray[4*i+1]=inttransform.output[1];
        chararray[4*i+2]=inttransform.output[2];
        chararray[4*i+3]=inttransform.output[3];
        }
}


static void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2)
/* function to compress "in" to "out" reducing size from nn to nn2 */
{
#if defined(USE_GZDIR) || defined(USE_ZLIB)
    int ntemp=0;

    /* in and out of z-stream */
    unsigned char inz[CHUNK];
    unsigned char outz[CHUNK];

    /* compression level */
    int level = Z_DEFAULT_COMPRESSION;
    int ret,flush;
    int i,j,k;

    /* zlib compression stream */
    z_stream strm;

    /* hope compressed data will be <= uncompressed */
    *out = (unsigned char *)malloc(sizeof(unsigned char)*nn);

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    /* zlib init */
    ret = deflateInit(&strm, level);
    if (ret == Z_OK){
        i=0;     // position in "in" array
        do{
            j=0; // position in "inz"
            do{
                inz[j++]=in[i++];
            } while((j<CHUNK) && (i<nn)); // stopps if "inz"-buffer is full or "in" array empty
            strm.avail_in=j;              // set number of input chars

            flush = (i==nn) ? Z_FINISH : Z_NO_FLUSH; // done?
            strm.next_in = inz;           // set input buffer

            do{
                strm.avail_out = CHUNK;   // set number of max output chars
                strm.next_out = outz;     // set output buffer

                /* zlib compress */
                ret = deflate(&strm, flush);
                assert(ret != Z_STREAM_ERROR);

                /* zlib changed strm.avail_out=CHUNK
                 to the number of chars we can NOT use
                 in outz */

                for (k=0;k<CHUNK-strm.avail_out;k++){
                    (*out)[ntemp+k]=outz[k];
                }

                /* increase position in "out" */
                ntemp+=(CHUNK-strm.avail_out);
            }while(strm.avail_out==0);
            assert(strm.avail_in == 0);

        }while (flush != Z_FINISH);
    }
    else{fprintf(stderr,"Error during compression init\n");}

    // now we know how short "out" should be!
    *nn2=ntemp;
    *out = (unsigned char *)realloc(*out, sizeof(unsigned char)*ntemp);

    (void)deflateEnd(&strm);
#endif

    return;
}

static void base64(unsigned char * in, int nn, unsigned char* out)
{
    /*takes *in*-array and "in"-length-"nn" and fills "out"-array
    with base64(in) "out" needs to be big enough!!!
    length(out) >= 4* |^ nn/3.0 ^| */
    char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int len;
    int i;

    for (i=0; i < nn; i+=3){

        len = (3 < nn-i ? 3 : nn-i);
        if (len >= 3){
        /* normal base64 encoding */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[ ((in[i+1] & 0x0f) << 2) | ((in[i+2] & 0xc0) >> 6)];
            out[i/3*4+3] = cb64[ in[i+2] & 0x3f ];
        } else if (len == 2){
        /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[((in[i+1] & 0x0f) << 2)];
            out[i/3*4+3] = (unsigned char) '=';
        } else if (len == 1){
        /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) ];
            out[i/3*4+2] = (unsigned char) '=';
            out[i/3*4+3] = (unsigned char) '=';
        }
    }
}


static void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out)
{
    /* writing vtk compatible zlib compressed base64 encoded data to "out" */
    int i;
    unsigned char * b64head;
    int b64bodylength;
    unsigned char * b64body;
    /* header of data */
    unsigned char * charhead = (unsigned char *)malloc(sizeof(unsigned char)*16);
    /* - consists of "1" (number of pieces) */
    /* - original datalength in byte */
    /* - original datalength in byte */
    /* - new datalength after z-lib compression */
    int * headInts= (int *)malloc(sizeof(int)*4);
    headInts[0]=1;
    headInts[1]=orinn;
    headInts[2]=orinn;
    headInts[3]=nn;
    // transform to unsigned char
    IntToUnsignedChar(headInts,4,charhead);

    // base64: 16byte -> 24byte
    b64head = (unsigned char *)malloc(sizeof(unsigned char)*24);
    // fills b64head
    base64(charhead, 16, b64head);

    // base64 data
    b64bodylength = 4*ceil((double) nn/3.0);
    b64body = (unsigned char *)malloc(sizeof(unsigned char)*b64bodylength);
    // writes base64 data to b64body
    base64(in,nn,b64body);

    // combines header and body
    for (i=0; i<24 ; i++){
        out[i]=b64head[i];
    }

    for (i=0; i<b64bodylength ; i++){
        out[24+i]=b64body[i];
    }

    if(b64body){free(b64body);}
    if(b64head){free(b64head);}
    if(headInts){free(headInts);}
    if(charhead){free(charhead);}
}

static void write_binary_array(int nn, float* array, FILE * f)
{
    /* writes vtk-data array of floats and performs zip and base64 encoding */
    int chararraylength=4*nn;	/* nn floats -> 4*nn unsigned chars */
    unsigned char * chararray = (unsigned char *)malloc(chararraylength * sizeof(unsigned char));
    int compressedarraylength = 0;
    unsigned char * compressedarray;
    unsigned char ** pointertocompressedarray= &compressedarray;
    int base64plusheadlength;
    unsigned char * base64plusheadarray;

    FloatToUnsignedChar(array,nn,chararray);

    /* compression routine */
    zlibcompress(chararray,chararraylength,pointertocompressedarray,&compressedarraylength);

    /* special header for zip compressed and bas64 encoded data
    header needs 4 int32 = 16 byte -> 24 byte due to base64 (4*16/3) */
    base64plusheadlength = 24 + 4*ceil((double) compressedarraylength/3.0);
    base64plusheadarray = (unsigned char *)malloc(sizeof(unsigned char)* base64plusheadlength);

    /* fills base64plusheadarray with everything ready for simple writing */
    base64plushead(compressedarray,compressedarraylength, chararraylength, base64plusheadarray);

    fwrite(base64plusheadarray,sizeof(unsigned char),base64plusheadlength,f);
    fprintf(f,"\n");
    free(chararray);
    free(base64plusheadarray);
    free(compressedarray);
}

static void write_binary_uint8_array(int nn, const unsigned char* array, FILE * f)
{
    int compressedarraylength = 0;
    unsigned char * compressedarray = NULL;
    int base64plusheadlength;
    unsigned char * base64plusheadarray;

    zlibcompress((unsigned char*)array, nn, &compressedarray, &compressedarraylength);

    base64plusheadlength = 24 + 4 * (int)ceil((double)compressedarraylength / 3.0);
    base64plusheadarray = (unsigned char *)malloc(sizeof(unsigned char) * base64plusheadlength);

    base64plushead(compressedarray, compressedarraylength, nn, base64plusheadarray);

    fwrite(base64plusheadarray, sizeof(unsigned char), base64plusheadlength, f);
    fprintf(f, "\n");

    free(base64plusheadarray);
    free(compressedarray);
}

static void write_binary_int_array(int nn, int* array, FILE * f)
{
    int chararraylength = 4 * nn;
    unsigned char * chararray = (unsigned char *)malloc(chararraylength * sizeof(unsigned char));
    int compressedarraylength = 0;
    unsigned char * compressedarray = NULL;
    int base64plusheadlength;
    unsigned char * base64plusheadarray;

    IntToUnsignedChar(array, nn, chararray);
    zlibcompress(chararray, chararraylength, &compressedarray, &compressedarraylength);

    base64plusheadlength = 24 + 4 * (int)ceil((double)compressedarraylength / 3.0);
    base64plusheadarray = (unsigned char *)malloc(sizeof(unsigned char) * base64plusheadlength);

    base64plushead(compressedarray, compressedarraylength, chararraylength, base64plusheadarray);

    fwrite(base64plusheadarray, sizeof(unsigned char), base64plusheadlength, f);
    fprintf(f, "\n");

    free(chararray);
    free(base64plusheadarray);
    free(compressedarray);
}

static void write_ascii_int_array(int nn, int perLine, const int *array, FILE *fp)
{
    int i;
    for (i = 0; i < nn; i++) {
        fprintf(fp, "%d%s", array[i], ((i + 1) % perLine == 0 || i == nn - 1) ? "\n" : " ");
    }
}

/**********************************************************************/

static void vtk_output_proc(struct All_variables *E, int cycles)
{
    char output_file[255];
    FILE *fp;
    snprintf(output_file, 255, "%s.proc%d.%d.vts",
             E->control.data_file, E->parallel.me, cycles);
    fp = output_open(output_file, "w");

    /* first, write volume data to vts file */
    vts_file_header(E, fp);

    /* write node-based field */
    vtk_point_data_header(E, fp);
    vtk_output_temp(E, fp);

    vtk_output_velo(E, fp);

    vtk_output_visc(E, fp);

    if (E->output.stress)
        vtk_output_stress(E, fp);

    if (E->output.comp_nd && E->composition.on)
        vtk_output_comp_nd(E, fp);

    if (E->output.surf)
        vtk_output_surf(E, fp, cycles);

    vtk_point_data_trailer(E, fp);

    /* write element-based field */
    vtk_cell_data_header(E, fp);
    /* TODO: comp_el, heating */
    vtk_cell_data_trailer(E, fp);

    /* write coordinate */
    vtk_output_coord(E, fp);

    vts_file_trailer(E, fp);
    fclose(fp);

    /* then, write other type of data */

    if (E->parallel.me == 0) {
        if (E->sphere.caps == 12) {
            /* VTK multiblock file */
            write_vtm(E, cycles);
            /* LLNL VisIt multiblock file */
            write_visit(E, cycles);
        }
        else
            write_pvts(E, cycles);

        /* ParaView PVD collection file */
        write_pvd(E, cycles);
    }
    return;
}

static void vtk_output_cap(struct All_variables *E, int cycles)
{
    int i, j, k;
    int nproc_cap = E->parallel.nprocx * E->parallel.nprocy * E->parallel.nprocz;
    int cap_id = E->parallel.me / nproc_cap;
    int rank_in_cap = E->parallel.me % nproc_cap;

    MPI_Comm cap_comm;
    MPI_Comm_split(E->parallel.world, cap_id, rank_in_cap, &cap_comm);

    int local_nno = E->lmesh.nno;
    int nox = E->lmesh.nox;
    int noy = E->lmesh.noy;
    int noz = E->lmesh.noz;

    int NOX = E->mesh.nox;
    int NOY = E->mesh.noy;
    int NOZ = E->mesh.noz;
    int NNO_cap = NOX * NOY * NOZ;

    const int has_comp_nd = (E->output.comp_nd && E->composition.on && E->composition.ncomp > 0) ? 1 : 0;
    const int ncomp = has_comp_nd ? E->composition.ncomp : 0;

    /* 1. Extract local fields */
    float *local_coord = (float *)malloc(local_nno * 3 * sizeof(float));
    float *local_temp = (float *)malloc(local_nno * sizeof(float));
    float *local_velo = (float *)malloc(local_nno * 3 * sizeof(float));
    float *local_visc = (float *)malloc(local_nno * sizeof(float));
    float *local_stress = NULL;
    float **local_comp = NULL;
    float *local_surf = NULL;

    for (i = 0; i < local_nno; i++) {
        local_coord[i * 3 + 0] = (float)E->x[1][1][i + 1];
        local_coord[i * 3 + 1] = (float)E->x[1][2][i + 1];
        local_coord[i * 3 + 2] = (float)E->x[1][3][i + 1];
        local_temp[i] = (float)*(E->T[1] + i + 1);
        local_visc[i] = E->VI[E->mesh.levmax][1][i + 1];
    }

    const int lev = E->mesh.levmax;
    for (i = 1; i <= local_nno; i++) {
        double sint = E->SinCos[lev][1][0][i];
        double sinf = E->SinCos[lev][1][1][i];
        double cost = E->SinCos[lev][1][2][i];
        double cosf = E->SinCos[lev][1][3][i];
        double v1 = E->sphere.cap[1].V[1][i];
        double v2 = E->sphere.cap[1].V[2][i];
        double v3 = E->sphere.cap[1].V[3][i];
        local_velo[(i - 1) * 3 + 0] = (float)(v1 * cost * cosf - v2 * sinf + v3 * sint * cosf);
        local_velo[(i - 1) * 3 + 1] = (float)(v1 * cost * sinf + v2 * cosf + v3 * sint * sinf);
        local_velo[(i - 1) * 3 + 2] = (float)(-v1 * sint + v3 * cost);
    }

    if (E->output.stress) {
        float *SXX[NCS], *SYY[NCS], *SZZ[NCS], *SXY[NCS], *SXZ[NCS], *SZY[NCS];
        float *divv[NCS], *vorv[NCS];
        allocate_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
        compute_nodal_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
        free_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);

        local_stress = (float *)malloc(local_nno * 6 * sizeof(float));
        for (i = 0; i < local_nno * 6; i++) {
            local_stress[i] = E->gstress[1][i + 1];
        }
    }

    if (has_comp_nd) {
        local_comp = (float **)malloc(ncomp * sizeof(float *));
        for (int c = 0; c < ncomp; c++) {
            local_comp[c] = (float *)malloc(local_nno * sizeof(float));
            for (i = 0; i < local_nno; i++) {
                local_comp[c][i] = (float)E->composition.comp_node[1][c][i + 1];
            }
        }
    }

    if (E->output.surf) {
        if ((E->output.write_q_files == 0) || (cycles == 0) ||
            (cycles % E->output.write_q_files) != 0)
            heat_flux(E);
        if (E->control.use_cbf_topo) {
            get_CBF_topo(E, E->slice.tpg, E->slice.tpgb);
        } else {
            get_STD_topo(E, E->slice.tpg, E->slice.tpgb, E->slice.divg, E->slice.vort, cycles);
        }
        local_surf = (float *)calloc(local_nno, sizeof(float));
        if (E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
            for (i = 1; i <= E->lmesh.nsf; i++) {
                local_surf[i * noz - 1] = E->control.pseudo_free_surf ? E->slice.freesurf[1][i] : E->slice.tpg[1][i];
            }
        }
    }

    /* 2. Gather to cap leader */
    float *cap_coord = NULL;
    float *cap_temp = NULL;
    float *cap_velo = NULL;
    float *cap_visc = NULL;
    float *cap_stress = NULL;
    float **cap_comp = NULL;
    float *cap_surf = NULL;

    if (rank_in_cap == 0) {
        cap_coord = (float *)malloc(NNO_cap * 3 * sizeof(float));
        cap_temp = (float *)malloc(NNO_cap * sizeof(float));
        cap_velo = (float *)malloc(NNO_cap * 3 * sizeof(float));
        cap_visc = (float *)malloc(NNO_cap * sizeof(float));
        if (E->output.stress) {
            cap_stress = (float *)malloc(NNO_cap * 6 * sizeof(float));
        }
        if (has_comp_nd) {
            cap_comp = (float **)malloc(ncomp * sizeof(float *));
            for (int c = 0; c < ncomp; c++) {
                cap_comp[c] = (float *)malloc(NNO_cap * sizeof(float));
            }
        }
        if (E->output.surf) {
            cap_surf = (float *)calloc(NNO_cap, sizeof(float));
        }

        /* Insert rank 0's own data */
        int xs = 0, ys = 0, zs = 0;
        for (j = 0; j < noy; j++) {
            for (i = 0; i < nox; i++) {
                for (k = 0; k < noz; k++) {
                    int l_idx = k + i * noz + j * noz * nox;
                    int gx = xs + i;
                    int gy = ys + j;
                    int gz = zs + k;
                    int c_idx = gz + gx * NOZ + gy * NOZ * NOX;
                    cap_temp[c_idx] = local_temp[l_idx];
                    cap_visc[c_idx] = local_visc[l_idx];
                    cap_coord[3 * c_idx + 0] = local_coord[3 * l_idx + 0];
                    cap_coord[3 * c_idx + 1] = local_coord[3 * l_idx + 1];
                    cap_coord[3 * c_idx + 2] = local_coord[3 * l_idx + 2];
                    cap_velo[3 * c_idx + 0] = local_velo[3 * l_idx + 0];
                    cap_velo[3 * c_idx + 1] = local_velo[3 * l_idx + 1];
                    cap_velo[3 * c_idx + 2] = local_velo[3 * l_idx + 2];
                    if (cap_stress && local_stress) {
                        for (int c = 0; c < 6; c++)
                            cap_stress[6 * c_idx + c] = local_stress[6 * l_idx + c];
                    }
                    if (cap_comp && local_comp) {
                        for (int c = 0; c < ncomp; c++)
                            cap_comp[c][c_idx] = local_comp[c][l_idx];
                    }
                    if (cap_surf && local_surf) {
                        cap_surf[c_idx] = local_surf[l_idx];
                    }
                }
            }
        }

        /* Receive from workers */
        float *recv_buf = (float *)malloc(local_nno * (E->output.stress ? 6 : 3) * sizeof(float));
        for (int r = 1; r < nproc_cap; r++) {
            int r_pz = r % E->parallel.nprocz;
            int r_px = ((r - r_pz) / E->parallel.nprocz) % E->parallel.nprocx;
            int r_py = (r - r_pz) / (E->parallel.nprocz * E->parallel.nprocx);
            int r_xs = r_px * (nox - 1);
            int r_ys = r_py * (noy - 1);
            int r_zs = r_pz * (noz - 1);

            MPI_Status st;
            /* Recv temp */
            MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 101, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_temp[c_idx] = recv_buf[l_idx];
                    }
                }
            }

            /* Recv visc */
            MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 102, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_visc[c_idx] = recv_buf[l_idx];
                    }
                }
            }

            /* Recv coord */
            MPI_Recv(recv_buf, local_nno * 3, MPI_FLOAT, r, 103, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_coord[3 * c_idx + 0] = recv_buf[3 * l_idx + 0];
                        cap_coord[3 * c_idx + 1] = recv_buf[3 * l_idx + 1];
                        cap_coord[3 * c_idx + 2] = recv_buf[3 * l_idx + 2];
                    }
                }
            }

            /* Recv velo */
            MPI_Recv(recv_buf, local_nno * 3, MPI_FLOAT, r, 104, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_velo[3 * c_idx + 0] = recv_buf[3 * l_idx + 0];
                        cap_velo[3 * c_idx + 1] = recv_buf[3 * l_idx + 1];
                        cap_velo[3 * c_idx + 2] = recv_buf[3 * l_idx + 2];
                    }
                }
            }

            /* Recv stress */
            if (E->output.stress) {
                MPI_Recv(recv_buf, local_nno * 6, MPI_FLOAT, r, 105, cap_comm, &st);
                for (j = 0; j < noy; j++) {
                    for (i = 0; i < nox; i++) {
                        for (k = 0; k < noz; k++) {
                            int l_idx = k + i * noz + j * noz * nox;
                            int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                            for (int c = 0; c < 6; c++)
                                cap_stress[6 * c_idx + c] = recv_buf[6 * l_idx + c];
                        }
                    }
                }
            }

            /* Recv comp_nd */
            if (has_comp_nd && cap_comp) {
                for (int c = 0; c < ncomp; c++) {
                    MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 110 + c, cap_comm, &st);
                    for (j = 0; j < noy; j++) {
                        for (i = 0; i < nox; i++) {
                            for (k = 0; k < noz; k++) {
                                int l_idx = k + i * noz + j * noz * nox;
                                int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                                cap_comp[c][c_idx] = recv_buf[l_idx];
                            }
                        }
                    }
                }
            }

            /* Recv surface */
            if (E->output.surf && cap_surf) {
                MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 106, cap_comm, &st);
                for (j = 0; j < noy; j++) {
                    for (i = 0; i < nox; i++) {
                        for (k = 0; k < noz; k++) {
                            int l_idx = k + i * noz + j * noz * nox;
                            int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                            cap_surf[c_idx] = recv_buf[l_idx];
                        }
                    }
                }
            }
        }
        free(recv_buf);

        /* Write cap VTS file */
        char output_file[255];
        snprintf(output_file, 255, "%s.cap%d.%d.vts", E->control.data_file, cap_id, cycles);
        FILE *fp = output_open(output_file, "w");

        char extent[128];
        snprintf(extent, sizeof(extent), "0 %d 0 %d 0 %d", NOZ - 1, NOX - 1, NOY - 1);
        fprintf(fp, "<?xml version=\"1.0\"?>\n");
        fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n");
        fprintf(fp, "  <StructuredGrid WholeExtent=\"%s\">\n", extent);
        fprintf(fp, "    <Piece Extent=\"%s\">\n", extent);
        fprintf(fp, "      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n");

        /* Temperature */
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(NNO_cap, cap_temp, fp);
        else
            write_ascii_array(NNO_cap, 1, cap_temp, fp);
        fprintf(fp, "        </DataArray>\n");

        /* Velocity */
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(NNO_cap * 3, cap_velo, fp);
        else
            write_ascii_array(NNO_cap * 3, 3, cap_velo, fp);
        fprintf(fp, "        </DataArray>\n");

        /* Viscosity */
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(NNO_cap, cap_visc, fp);
        else
            write_ascii_array(NNO_cap, 1, cap_visc, fp);
        fprintf(fp, "        </DataArray>\n");

        /* Stress */
        if (E->output.stress && cap_stress) {
            fprintf(fp, "        <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\">\n", E->output.vtk_format);
            if (strcmp(E->output.vtk_format, "binary") == 0)
                write_binary_array(NNO_cap * 6, cap_stress, fp);
            else
                write_ascii_array(NNO_cap * 6, 6, cap_stress, fp);
            fprintf(fp, "        </DataArray>\n");
        }

        /* Nodal compositions (comp_nd) */
        if (has_comp_nd && cap_comp) {
            for (int c = 0; c < ncomp; c++) {
                fprintf(fp, "        <DataArray type=\"Float32\" Name=\"composition%d\" format=\"%s\">\n",
                        c + 1, E->output.vtk_format);
                if (strcmp(E->output.vtk_format, "binary") == 0)
                    write_binary_array(NNO_cap, cap_comp[c], fp);
                else
                    write_ascii_array(NNO_cap, 1, cap_comp[c], fp);
                fprintf(fp, "        </DataArray>\n");
            }
        }

        /* Surface topography */
        if (E->output.surf && cap_surf) {
            fprintf(fp, "        <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\">\n", E->output.vtk_format);
            if (strcmp(E->output.vtk_format, "binary") == 0)
                write_binary_array(NNO_cap, cap_surf, fp);
            else
                write_ascii_array(NNO_cap, 1, cap_surf, fp);
            fprintf(fp, "        </DataArray>\n");
        }

        fprintf(fp, "      </PointData>\n");
        fprintf(fp, "      <CellData>\n      </CellData>\n");
        fprintf(fp, "      <Points>\n");
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(NNO_cap * 3, cap_coord, fp);
        else
            write_ascii_array(NNO_cap * 3, 3, cap_coord, fp);
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "      </Points>\n");
        fprintf(fp, "    </Piece>\n");
        fprintf(fp, "  </StructuredGrid>\n");
        fprintf(fp, "</VTKFile>\n");
        fclose(fp);

        free(cap_coord);
        free(cap_temp);
        free(cap_velo);
        free(cap_visc);
        if (cap_stress) free(cap_stress);
        if (cap_comp) {
            for (int c = 0; c < ncomp; c++) free(cap_comp[c]);
            free(cap_comp);
        }
        if (cap_surf) free(cap_surf);

    } else {
        /* Worker sends local buffers */
        MPI_Send(local_temp, local_nno, MPI_FLOAT, 0, 101, cap_comm);
        MPI_Send(local_visc, local_nno, MPI_FLOAT, 0, 102, cap_comm);
        MPI_Send(local_coord, local_nno * 3, MPI_FLOAT, 0, 103, cap_comm);
        MPI_Send(local_velo, local_nno * 3, MPI_FLOAT, 0, 104, cap_comm);
        if (E->output.stress && local_stress) {
            MPI_Send(local_stress, local_nno * 6, MPI_FLOAT, 0, 105, cap_comm);
        }
        if (has_comp_nd && local_comp) {
            for (int c = 0; c < ncomp; c++) {
                MPI_Send(local_comp[c], local_nno, MPI_FLOAT, 0, 110 + c, cap_comm);
            }
        }
        if (E->output.surf && local_surf) {
            MPI_Send(local_surf, local_nno, MPI_FLOAT, 0, 106, cap_comm);
        }
    }

    free(local_coord);
    free(local_temp);
    free(local_velo);
    free(local_visc);
    if (local_stress) free(local_stress);
    if (local_comp) {
        for (int c = 0; c < ncomp; c++) free(local_comp[c]);
        free(local_comp);
    }
    if (local_surf) free(local_surf);

    MPI_Comm_free(&cap_comm);

    /* Global rank 0 writes index files */
    if (E->parallel.me == 0) {
        if (E->sphere.caps == 12) {
            write_vtm_caps(E, cycles);
            write_visit_caps(E, cycles);
        }
        write_pvd(E, cycles);
    }
}

static void vtk_output_single(struct All_variables *E, int cycles)
{
    int i, j, k;
    int nproc_cap = E->parallel.nprocx * E->parallel.nprocy * E->parallel.nprocz;
    int cap_id = E->parallel.me / nproc_cap;
    int rank_in_cap = E->parallel.me % nproc_cap;

    MPI_Comm cap_comm;
    MPI_Comm_split(E->parallel.world, cap_id, rank_in_cap, &cap_comm);

    int local_nno = E->lmesh.nno;
    int nox = E->lmesh.nox;
    int noy = E->lmesh.noy;
    int noz = E->lmesh.noz;

    int NOX = E->mesh.nox;
    int NOY = E->mesh.noy;
    int NOZ = E->mesh.noz;
    int NNO_cap = NOX * NOY * NOZ;
    int caps = E->sphere.caps;

    const int has_comp_nd = (E->output.comp_nd && E->composition.on && E->composition.ncomp > 0) ? 1 : 0;
    const int ncomp = has_comp_nd ? E->composition.ncomp : 0;

    /* 1. Extract local fields */
    float *local_coord = (float *)malloc(local_nno * 3 * sizeof(float));
    float *local_temp = (float *)malloc(local_nno * sizeof(float));
    float *local_velo = (float *)malloc(local_nno * 3 * sizeof(float));
    float *local_visc = (float *)malloc(local_nno * sizeof(float));
    float *local_stress = NULL;
    float **local_comp = NULL;
    float *local_surf = NULL;

    for (i = 0; i < local_nno; i++) {
        local_coord[i * 3 + 0] = (float)E->x[1][1][i + 1];
        local_coord[i * 3 + 1] = (float)E->x[1][2][i + 1];
        local_coord[i * 3 + 2] = (float)E->x[1][3][i + 1];
        local_temp[i] = (float)*(E->T[1] + i + 1);
        local_visc[i] = E->VI[E->mesh.levmax][1][i + 1];
    }

    const int lev = E->mesh.levmax;
    for (i = 1; i <= local_nno; i++) {
        double sint = E->SinCos[lev][1][0][i];
        double sinf = E->SinCos[lev][1][1][i];
        double cost = E->SinCos[lev][1][2][i];
        double cosf = E->SinCos[lev][1][3][i];
        double v1 = E->sphere.cap[1].V[1][i];
        double v2 = E->sphere.cap[1].V[2][i];
        double v3 = E->sphere.cap[1].V[3][i];
        local_velo[(i - 1) * 3 + 0] = (float)(v1 * cost * cosf - v2 * sinf + v3 * sint * cosf);
        local_velo[(i - 1) * 3 + 1] = (float)(v1 * cost * sinf + v2 * cosf + v3 * sint * sinf);
        local_velo[(i - 1) * 3 + 2] = (float)(-v1 * sint + v3 * cost);
    }

    if (E->output.stress) {
        float *SXX[NCS], *SYY[NCS], *SZZ[NCS], *SXY[NCS], *SXZ[NCS], *SZY[NCS];
        float *divv[NCS], *vorv[NCS];
        allocate_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
        compute_nodal_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
        free_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);

        local_stress = (float *)malloc(local_nno * 6 * sizeof(float));
        for (i = 0; i < local_nno * 6; i++) {
            local_stress[i] = E->gstress[1][i + 1];
        }
    }

    if (has_comp_nd) {
        local_comp = (float **)malloc(ncomp * sizeof(float *));
        for (int c = 0; c < ncomp; c++) {
            local_comp[c] = (float *)malloc(local_nno * sizeof(float));
            for (i = 0; i < local_nno; i++) {
                local_comp[c][i] = (float)E->composition.comp_node[1][c][i + 1];
            }
        }
    }

    if (E->output.surf) {
        if ((E->output.write_q_files == 0) || (cycles == 0) ||
            (cycles % E->output.write_q_files) != 0)
            heat_flux(E);
        if (E->control.use_cbf_topo) {
            get_CBF_topo(E, E->slice.tpg, E->slice.tpgb);
        } else {
            get_STD_topo(E, E->slice.tpg, E->slice.tpgb, E->slice.divg, E->slice.vort, cycles);
        }
        local_surf = (float *)calloc(local_nno, sizeof(float));
        if (E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
            for (i = 1; i <= E->lmesh.nsf; i++) {
                local_surf[i * noz - 1] = E->control.pseudo_free_surf ? E->slice.freesurf[1][i] : E->slice.tpg[1][i];
            }
        }
    }

    /* 2. Assemble on cap leader */
    float *cap_coord = NULL;
    float *cap_temp = NULL;
    float *cap_velo = NULL;
    float *cap_visc = NULL;
    float *cap_stress = NULL;
    float **cap_comp = NULL;
    float *cap_surf = NULL;

    if (rank_in_cap == 0) {
        cap_coord = (float *)malloc(NNO_cap * 3 * sizeof(float));
        cap_temp = (float *)malloc(NNO_cap * sizeof(float));
        cap_velo = (float *)malloc(NNO_cap * 3 * sizeof(float));
        cap_visc = (float *)malloc(NNO_cap * sizeof(float));
        if (E->output.stress) {
            cap_stress = (float *)malloc(NNO_cap * 6 * sizeof(float));
        }
        if (has_comp_nd) {
            cap_comp = (float **)malloc(ncomp * sizeof(float *));
            for (int c = 0; c < ncomp; c++) {
                cap_comp[c] = (float *)malloc(NNO_cap * sizeof(float));
            }
        }
        if (E->output.surf) {
            cap_surf = (float *)calloc(NNO_cap, sizeof(float));
        }

        int xs = 0, ys = 0, zs = 0;
        for (j = 0; j < noy; j++) {
            for (i = 0; i < nox; i++) {
                for (k = 0; k < noz; k++) {
                    int l_idx = k + i * noz + j * noz * nox;
                    int gx = xs + i;
                    int gy = ys + j;
                    int gz = zs + k;
                    int c_idx = gz + gx * NOZ + gy * NOZ * NOX;
                    cap_temp[c_idx] = local_temp[l_idx];
                    cap_visc[c_idx] = local_visc[l_idx];
                    cap_coord[3 * c_idx + 0] = local_coord[3 * l_idx + 0];
                    cap_coord[3 * c_idx + 1] = local_coord[3 * l_idx + 1];
                    cap_coord[3 * c_idx + 2] = local_coord[3 * l_idx + 2];
                    cap_velo[3 * c_idx + 0] = local_velo[3 * l_idx + 0];
                    cap_velo[3 * c_idx + 1] = local_velo[3 * l_idx + 1];
                    cap_velo[3 * c_idx + 2] = local_velo[3 * l_idx + 2];
                    if (cap_stress && local_stress) {
                        for (int c = 0; c < 6; c++)
                            cap_stress[6 * c_idx + c] = local_stress[6 * l_idx + c];
                    }
                    if (cap_comp && local_comp) {
                        for (int c = 0; c < ncomp; c++)
                            cap_comp[c][c_idx] = local_comp[c][l_idx];
                    }
                    if (cap_surf && local_surf) {
                        cap_surf[c_idx] = local_surf[l_idx];
                    }
                }
            }
        }

        float *recv_buf = (float *)malloc(local_nno * (E->output.stress ? 6 : 3) * sizeof(float));
        for (int r = 1; r < nproc_cap; r++) {
            int r_pz = r % E->parallel.nprocz;
            int r_px = ((r - r_pz) / E->parallel.nprocz) % E->parallel.nprocx;
            int r_py = (r - r_pz) / (E->parallel.nprocz * E->parallel.nprocx);
            int r_xs = r_px * (nox - 1);
            int r_ys = r_py * (noy - 1);
            int r_zs = r_pz * (noz - 1);

            MPI_Status st;
            MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 201, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_temp[c_idx] = recv_buf[l_idx];
                    }
                }
            }

            MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 202, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_visc[c_idx] = recv_buf[l_idx];
                    }
                }
            }

            MPI_Recv(recv_buf, local_nno * 3, MPI_FLOAT, r, 203, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_coord[3 * c_idx + 0] = recv_buf[3 * l_idx + 0];
                        cap_coord[3 * c_idx + 1] = recv_buf[3 * l_idx + 1];
                        cap_coord[3 * c_idx + 2] = recv_buf[3 * l_idx + 2];
                    }
                }
            }

            MPI_Recv(recv_buf, local_nno * 3, MPI_FLOAT, r, 204, cap_comm, &st);
            for (j = 0; j < noy; j++) {
                for (i = 0; i < nox; i++) {
                    for (k = 0; k < noz; k++) {
                        int l_idx = k + i * noz + j * noz * nox;
                        int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                        cap_velo[3 * c_idx + 0] = recv_buf[3 * l_idx + 0];
                        cap_velo[3 * c_idx + 1] = recv_buf[3 * l_idx + 1];
                        cap_velo[3 * c_idx + 2] = recv_buf[3 * l_idx + 2];
                    }
                }
            }

            if (E->output.stress) {
                MPI_Recv(recv_buf, local_nno * 6, MPI_FLOAT, r, 205, cap_comm, &st);
                for (j = 0; j < noy; j++) {
                    for (i = 0; i < nox; i++) {
                        for (k = 0; k < noz; k++) {
                            int l_idx = k + i * noz + j * noz * nox;
                            int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                            for (int c = 0; c < 6; c++)
                                cap_stress[6 * c_idx + c] = recv_buf[6 * l_idx + c];
                        }
                    }
                }
            }

            if (has_comp_nd && cap_comp) {
                for (int c = 0; c < ncomp; c++) {
                    MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 210 + c, cap_comm, &st);
                    for (j = 0; j < noy; j++) {
                        for (i = 0; i < nox; i++) {
                            for (k = 0; k < noz; k++) {
                                int l_idx = k + i * noz + j * noz * nox;
                                int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                                cap_comp[c][c_idx] = recv_buf[l_idx];
                            }
                        }
                    }
                }
            }

            if (E->output.surf && cap_surf) {
                MPI_Recv(recv_buf, local_nno, MPI_FLOAT, r, 206, cap_comm, &st);
                for (j = 0; j < noy; j++) {
                    for (i = 0; i < nox; i++) {
                        for (k = 0; k < noz; k++) {
                            int l_idx = k + i * noz + j * noz * nox;
                            int c_idx = (r_zs + k) + (r_xs + i) * NOZ + (r_ys + j) * NOZ * NOX;
                            cap_surf[c_idx] = recv_buf[l_idx];
                        }
                    }
                }
            }
        }
        free(recv_buf);
    } else {
        MPI_Send(local_temp, local_nno, MPI_FLOAT, 0, 201, cap_comm);
        MPI_Send(local_visc, local_nno, MPI_FLOAT, 0, 202, cap_comm);
        MPI_Send(local_coord, local_nno * 3, MPI_FLOAT, 0, 203, cap_comm);
        MPI_Send(local_velo, local_nno * 3, MPI_FLOAT, 0, 204, cap_comm);
        if (E->output.stress && local_stress) {
            MPI_Send(local_stress, local_nno * 6, MPI_FLOAT, 0, 205, cap_comm);
        }
        if (has_comp_nd && local_comp) {
            for (int c = 0; c < ncomp; c++) {
                MPI_Send(local_comp[c], local_nno, MPI_FLOAT, 0, 210 + c, cap_comm);
            }
        }
        if (E->output.surf && local_surf) {
            MPI_Send(local_surf, local_nno, MPI_FLOAT, 0, 206, cap_comm);
        }
    }

    free(local_coord);
    free(local_temp);
    free(local_velo);
    free(local_visc);
    if (local_stress) free(local_stress);
    if (local_comp) {
        for (int c = 0; c < ncomp; c++) free(local_comp[c]);
        free(local_comp);
    }
    if (local_surf) free(local_surf);
    MPI_Comm_free(&cap_comm);

    /* 3. Gather all caps to global rank 0 and write monolithic .vtu */
    int total_points = caps * NNO_cap;
    int elx = NOX - 1;
    int ely = NOY - 1;
    int elz = NOZ - 1;
    int cells_per_cap = elx * ely * elz;
    int total_cells = caps * cells_per_cap;

    float *all_coord = NULL;
    float *all_temp = NULL;
    float *all_velo = NULL;
    float *all_visc = NULL;
    float *all_stress = NULL;
    float **all_comp = NULL;
    float *all_surf = NULL;

    if (E->parallel.me == 0) {
        all_coord = (float *)malloc(total_points * 3 * sizeof(float));
        all_temp = (float *)malloc(total_points * sizeof(float));
        all_velo = (float *)malloc(total_points * 3 * sizeof(float));
        all_visc = (float *)malloc(total_points * sizeof(float));
        if (E->output.stress) {
            all_stress = (float *)malloc(total_points * 6 * sizeof(float));
        }
        if (has_comp_nd) {
            all_comp = (float **)malloc(ncomp * sizeof(float *));
            for (int c = 0; c < ncomp; c++) {
                all_comp[c] = (float *)malloc(total_points * sizeof(float));
            }
        }
        if (E->output.surf) {
            all_surf = (float *)malloc(total_points * sizeof(float));
        }

        /* Copy cap 0 */
        memcpy(all_coord, cap_coord, NNO_cap * 3 * sizeof(float));
        memcpy(all_temp, cap_temp, NNO_cap * sizeof(float));
        memcpy(all_velo, cap_velo, NNO_cap * 3 * sizeof(float));
        memcpy(all_visc, cap_visc, NNO_cap * sizeof(float));
        if (all_stress && cap_stress) {
            memcpy(all_stress, cap_stress, NNO_cap * 6 * sizeof(float));
        }
        if (all_comp && cap_comp) {
            for (int c = 0; c < ncomp; c++)
                memcpy(all_comp[c], cap_comp[c], NNO_cap * sizeof(float));
        }
        if (all_surf && cap_surf) {
            memcpy(all_surf, cap_surf, NNO_cap * sizeof(float));
        }

        /* Receive remaining caps from cap leaders */
        for (int c = 1; c < caps; c++) {
            int leader_rank = c * nproc_cap;
            MPI_Status st;
            MPI_Recv(all_temp + c * NNO_cap, NNO_cap, MPI_FLOAT, leader_rank, 301, E->parallel.world, &st);
            MPI_Recv(all_visc + c * NNO_cap, NNO_cap, MPI_FLOAT, leader_rank, 302, E->parallel.world, &st);
            MPI_Recv(all_coord + c * NNO_cap * 3, NNO_cap * 3, MPI_FLOAT, leader_rank, 303, E->parallel.world, &st);
            MPI_Recv(all_velo + c * NNO_cap * 3, NNO_cap * 3, MPI_FLOAT, leader_rank, 304, E->parallel.world, &st);
            if (all_stress) {
                MPI_Recv(all_stress + c * NNO_cap * 6, NNO_cap * 6, MPI_FLOAT, leader_rank, 305, E->parallel.world, &st);
            }
            if (all_comp) {
                for (int comp_i = 0; comp_i < ncomp; comp_i++) {
                    MPI_Recv(all_comp[comp_i] + c * NNO_cap, NNO_cap, MPI_FLOAT, leader_rank, 310 + comp_i, E->parallel.world, &st);
                }
            }
            if (all_surf) {
                MPI_Recv(all_surf + c * NNO_cap, NNO_cap, MPI_FLOAT, leader_rank, 306, E->parallel.world, &st);
            }
        }
    } else if (rank_in_cap == 0) {
        /* Other cap leaders send to global rank 0 */
        MPI_Send(cap_temp, NNO_cap, MPI_FLOAT, 0, 301, E->parallel.world);
        MPI_Send(cap_visc, NNO_cap, MPI_FLOAT, 0, 302, E->parallel.world);
        MPI_Send(cap_coord, NNO_cap * 3, MPI_FLOAT, 0, 303, E->parallel.world);
        MPI_Send(cap_velo, NNO_cap * 3, MPI_FLOAT, 0, 304, E->parallel.world);
        if (E->output.stress && cap_stress) {
            MPI_Send(cap_stress, NNO_cap * 6, MPI_FLOAT, 0, 305, E->parallel.world);
        }
        if (has_comp_nd && cap_comp) {
            for (int comp_i = 0; comp_i < ncomp; comp_i++) {
                MPI_Send(cap_comp[comp_i], NNO_cap, MPI_FLOAT, 0, 310 + comp_i, E->parallel.world);
            }
        }
        if (E->output.surf && cap_surf) {
            MPI_Send(cap_surf, NNO_cap, MPI_FLOAT, 0, 306, E->parallel.world);
        }
    }

    if (cap_coord) free(cap_coord);
    if (cap_temp) free(cap_temp);
    if (cap_velo) free(cap_velo);
    if (cap_visc) free(cap_visc);
    if (cap_stress) free(cap_stress);
    if (cap_comp) {
        for (int c = 0; c < ncomp; c++) free(cap_comp[c]);
        free(cap_comp);
    }
    if (cap_surf) free(cap_surf);

    /* 4. Write monolithic .vtu file on global rank 0 */
    if (E->parallel.me == 0) {
        char vtu_file[255];
        snprintf(vtu_file, 255, "%s.%d.vtu", E->control.data_file, cycles);
        FILE *fp = output_open(vtu_file, "w");

        fprintf(fp, "<?xml version=\"1.0\"?>\n");
        fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n");
        fprintf(fp, "  <UnstructuredGrid>\n");
        fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", total_points, total_cells);

        /* PointData */
        fprintf(fp, "      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n");
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(total_points, all_temp, fp);
        else
            write_ascii_array(total_points, 1, all_temp, fp);
        fprintf(fp, "        </DataArray>\n");

        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(total_points * 3, all_velo, fp);
        else
            write_ascii_array(total_points * 3, 3, all_velo, fp);
        fprintf(fp, "        </DataArray>\n");

        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(total_points, all_visc, fp);
        else
            write_ascii_array(total_points, 1, all_visc, fp);
        fprintf(fp, "        </DataArray>\n");

        if (E->output.stress && all_stress) {
            fprintf(fp, "        <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\">\n", E->output.vtk_format);
            if (strcmp(E->output.vtk_format, "binary") == 0)
                write_binary_array(total_points * 6, all_stress, fp);
            else
                write_ascii_array(total_points * 6, 6, all_stress, fp);
            fprintf(fp, "        </DataArray>\n");
        }

        if (has_comp_nd && all_comp) {
            for (int c = 0; c < ncomp; c++) {
                fprintf(fp, "        <DataArray type=\"Float32\" Name=\"composition%d\" format=\"%s\">\n",
                        c + 1, E->output.vtk_format);
                if (strcmp(E->output.vtk_format, "binary") == 0)
                    write_binary_array(total_points, all_comp[c], fp);
                else
                    write_ascii_array(total_points, 1, all_comp[c], fp);
                fprintf(fp, "        </DataArray>\n");
            }
        }

        if (E->output.surf && all_surf) {
            fprintf(fp, "        <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\">\n", E->output.vtk_format);
            if (strcmp(E->output.vtk_format, "binary") == 0)
                write_binary_array(total_points, all_surf, fp);
            else
                write_ascii_array(total_points, 1, all_surf, fp);
            fprintf(fp, "        </DataArray>\n");
        }
        fprintf(fp, "      </PointData>\n");

        /* Points */
        fprintf(fp, "      <Points>\n");
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(total_points * 3, all_coord, fp);
        else
            write_ascii_array(total_points * 3, 3, all_coord, fp);
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "      </Points>\n");

        /* Cells */
        int *connectivity = (int *)malloc(total_cells * 8 * sizeof(int));
        int *offsets = (int *)malloc(total_cells * sizeof(int));
        unsigned char *cell_types = (unsigned char *)malloc(total_cells * sizeof(unsigned char));

        int cell_cnt = 0;
        for (int c = 0; c < caps; c++) {
            int cap_base = c * NNO_cap;
            for (j = 0; j < ely; j++) {
                for (i = 0; i < elx; i++) {
                    for (k = 0; k < elz; k++) {
                        /* Hexahedron node connectivity: 8 vertices */
                        int n0 = k + i * NOZ + j * NOZ * NOX + cap_base;
                        int n1 = (k + 1) + i * NOZ + j * NOZ * NOX + cap_base;
                        int n2 = (k + 1) + (i + 1) * NOZ + j * NOZ * NOX + cap_base;
                        int n3 = k + (i + 1) * NOZ + j * NOZ * NOX + cap_base;
                        int n4 = k + i * NOZ + (j + 1) * NOZ * NOX + cap_base;
                        int n5 = (k + 1) + i * NOZ + (j + 1) * NOZ * NOX + cap_base;
                        int n6 = (k + 1) + (i + 1) * NOZ + (j + 1) * NOZ * NOX + cap_base;
                        int n7 = k + (i + 1) * NOZ + (j + 1) * NOZ * NOX + cap_base;

                        connectivity[cell_cnt * 8 + 0] = n0;
                        connectivity[cell_cnt * 8 + 1] = n1;
                        connectivity[cell_cnt * 8 + 2] = n2;
                        connectivity[cell_cnt * 8 + 3] = n3;
                        connectivity[cell_cnt * 8 + 4] = n4;
                        connectivity[cell_cnt * 8 + 5] = n5;
                        connectivity[cell_cnt * 8 + 6] = n6;
                        connectivity[cell_cnt * 8 + 7] = n7;

                        offsets[cell_cnt] = (cell_cnt + 1) * 8;
                        cell_types[cell_cnt] = 12; /* VTK_HEXAHEDRON */
                        cell_cnt++;
                    }
                }
            }
        }

        fprintf(fp, "      <Cells>\n");
        /* Connectivity */
        fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_int_array(total_cells * 8, connectivity, fp);
        else
            write_ascii_int_array(total_cells * 8, 8, connectivity, fp);
        fprintf(fp, "        </DataArray>\n");

        /* Offsets */
        fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_int_array(total_cells, offsets, fp);
        else
            write_ascii_int_array(total_cells, 8, offsets, fp);
        fprintf(fp, "        </DataArray>\n");

        /* Types */
        fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_uint8_array(total_cells, cell_types, fp);
        else {
            for (i = 0; i < total_cells; i++)
                fprintf(fp, "%d%s", cell_types[i], ((i + 1) % 16 == 0 || i == total_cells - 1) ? "\n" : " ");
        }
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "      </Cells>\n");

        fprintf(fp, "    </Piece>\n");
        fprintf(fp, "  </UnstructuredGrid>\n");
        fprintf(fp, "</VTKFile>\n");
        fclose(fp);

        free(connectivity);
        free(offsets);
        free(cell_types);

        free(all_coord);
        free(all_temp);
        free(all_velo);
        free(all_visc);
        if (all_stress) free(all_stress);
        if (all_comp) {
            for (int c = 0; c < ncomp; c++) free(all_comp[c]);
            free(all_comp);
        }
        if (all_surf) free(all_surf);

        /* Write PVD */
        write_pvd(E, cycles);
    }
}

void vtk_output(struct All_variables *E, int cycles)
{
    if (strcmp(E->output.vtk_output_mode, "proc") == 0) {
        vtk_output_proc(E, cycles);
    } else if (strcmp(E->output.vtk_output_mode, "single") == 0) {
        vtk_output_single(E, cycles);
    } else {
        /* default: "cap" */
        vtk_output_cap(E, cycles);
    }
}

