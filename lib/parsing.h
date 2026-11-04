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
#pragma once

void setup_parser(struct All_variables *E, const char *filename);
void shutdown_parser(struct All_variables *E);

int input_string(const char *name, char *value, const char *Default, int m);
int input_boolean(const char *name, int *value, const char *interpret, int m);
int input_int(const char *name, int *value, const char *interpret, int m);
int input_float(const char *name, float *value, const char *interpret, int m);
int input_double(const char *name, double *value, const char *interpret, int m);
int input_int_vector(const char *name, int number, int *value, int m);
int input_char_vector(const char *name, int number, char *value, int m);
int input_float_vector(const char *name, int number, float *value, int m);
int input_double_vector(const char *name, int number, double *value, int m);
