#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>

#include "VerletList.h"
#include "main.h"

using std::vector;

int VerletList(int N, double r_cut, double Lx,
			   vector<double> &x, vector<double> &y,
			   vector<double> &x_save, vector<double> &y_save,
			   vector<vector <int>> &VL, int VL_count_old, int first_call)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list_sq = r_factor * r_factor * r_cut_sq;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		n, m;

	int 		VL_count = 0;

	if (first_call == 0){
		double dr_max = 0.0;
		for (n = 0; n < N; n++){
			dx = x[n] - x_save[n];
			dy = y[n] - y_save[n];
			dx -= std::round(dx / Lx) * Lx;
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return VL_count_old;
	    }
	}

	double 		x1, y1;

	for (n = 0; n < N - 1; n++){
	    for (m = n + 1; m < N; m++){
			dx = x[m] - x[n];
			dy = y[m] - y[n];
			dx -= std::round(dx / Lx) * Lx;
			dr = dx * dx + dy * dy;
			if (dr < r_list_sq){
				VL[VL_count][0] = n;
				VL[VL_count][1] = m;
				VL_count++;
			}
		}
	}	

	for (n = 0; n < N; n++){
		x_save[n] = x[n];
		y_save[n] = y[n];
	}
	
	return VL_count;
}

int VerletList(int N, int N_up, double r_cut, double Lx,
			   vector<double> &xp, vector<double> &yp,
			   vector<double> &xp_save, vector<double> &yp_save,
			   vector<vector <int>> &VL, int VL_count_old, int first_call)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list_sq = r_factor * r_factor * r_cut_sq;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		n, m;

	int 		VL_count = 0;

	if (first_call == 0){
		double dr_max = 0.0;
		for (n = 0; n < N_up; n++){
			dx = xp_save[n] - xp[n];
			dy = yp_save[n] - yp[n];
			dx -= std::round(dx / Lx) * Lx;
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return VL_count_old;
	    }
	}

	double 		x1, y1;

	for (n = 0; n < N_up; n++){
	    for (m = n + 1; m < N; m++){
			dx = xp[m] - xp[n];
			dy = yp[m] - yp[n];
			dx -= std::round(dx / Lx) * Lx;
			dr = dx * dx + dy * dy;
			if (dr < r_list_sq){
				VL[VL_count][0] = n;
				VL[VL_count][1] = m;
				VL_count++;
			}
		}
	}	

	for (n = 0; n < N; n++){
		xp_save[n] = xp[n];
		yp_save[n] = yp[n];
	}

	return VL_count;
}