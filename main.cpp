#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <vector>

#include "main.h"
#include "MD.h"

using std::vector;

void Initialize(int N, double Dl, vector<double> &x, vector<double> &y, double Lx, int pid)
{
	std::mt19937 gen(pid); //Standard mersenne_twister_engine seeded with pid
	std::uniform_real_distribution<> dis(0.0, 1.0);
	Dl *= 1.01; // prevent particle from touching
	int				Nx = std::floor(Lx / Dl); // column number for the lattice sites
	int 			Ny = std::ceil((double)N / (double)Nx); // row number for the lattice sites
	vector<int> 	idx(Nx * Ny, 0);
	for (int i = 0; i < Nx * Ny; i++)
		idx[i] = i;
	shuffle(idx.begin(), idx.end(), std::default_random_engine(pid));
	int 			nx, ny;
	for (int i = 0; i < Nx * Ny; i++){
		if (idx[i] < N){
			ny = i / Nx;
			nx = i % Nx;
			x[idx[i]] = Dl * (double)nx + 0.001 * Dl * (dis(gen) - 0.5);
			y[idx[i]] = Dl * (double)ny + Dl / 2.0 + 0.001 * Dl * (dis(gen) - 0.5);
		}
	}

	return;
}


int main(int argc, char *argv[])
{
    int 			N = (int)atoi(argv[1]); // number of particles
    int          	mu_c = (int)atoi(argv[2]); // friction coefficient
	int				mu_e = (int)atoi(argv[3]); // friction coefficient
	int 			v0_c = (int)atoi(argv[4]);
	int 			v0_e = (int)atoi(argv[5]);
	int 			pid = (int)atoi(argv[6]); // random seed number

	double 			mu;
	if (mu_c == 0)
		mu = 0.0;
	else
		mu = (double)mu_c * std::pow(10.0, mu_e);

    int 			i, j, k, l;
    char 			dongdir[200] = "/Users/wangd/Documents/Yale/FrictionVibrationModes/RiverBed/Programs/Flow/Files/";
	
	// initialization
	double 			Dl = 1.4; // diameter of large particles
	double 			Ds = 1.0; // diameter of small particles
	vector<double> 	D(N, 0.0), D_sq(N, 0.0), M(N, 0.0), I(N, 0.0); // diameter, half big and half small
													 // mass and inertia
	for (i = 0; i < N / 2; i++){
		D[i] = Dl;
		D_sq[i] = Dl * Dl;
		M[i] = Dl * Dl;
		I[i] = 0.125 * M[i] * D_sq[i];
	}
	for (i = N / 2; i < N; i++){
		D[i] = Ds;
		D_sq[i] = Ds * Ds;
		M[i] = Ds * Ds;
		I[i] = 0.125 * M[i] * D_sq[i];
	}

	int 			N_height = 10; // Average 10 layers of particles in height
	double 			Lx = 0.5 * (Dl + Ds) * (double)N / (double)N_height; // box length, periodic

	// std::mt19937 gen(pid); //Standard mersenne_twister_engine seeded with pid
	// std::uniform_real_distribution<> dis(0.0, 1.0); // to randomly draw number from a uniform distribution within (0, 1)
	// generate initial configuration by random placing particles on lattice sites
	FILE*			otptfl;
	vector<double> 	x(N, 0.0), y(N, 0.0);
	char 			cmd[300], initfile[300], sedfile[300], rundir[300], flowdir[300];
	sprintf(rundir, "%sN_%04d/Mu_%dE%d/%05d/", dongdir, N, mu_c, mu_e, pid);
	sprintf(cmd, "mkdir -p %s", rundir);
	std::system(cmd);
	sprintf(initfile, "%sPos_Initial_%05d.txt", rundir, pid); // file name to save particle positions
	sprintf(sedfile, "%sPos_Sediment_%05d.txt", rundir, pid);
	sprintf(flowdir, "%sv0_%dE%d", rundir, v0_c, v0_e);
	sprintf(cmd, "mkdir -p %s", flowdir);
	std::system(cmd);
	
	std::ifstream posfile(initfile);
    if (posfile.good()){
		double Lx_temp;
		posfile >> Lx_temp;
		for (i = 0; i < N; i++){
			posfile >> x[i];
			posfile >> y[i];
		}
		posfile.close();
    }
	else{
		posfile.close();
		Initialize(N, Dl, x, y, Lx, pid);

		otptfl = fopen(initfile, "w");

		fprintf(otptfl, "%.32e\n", Lx);
		for (i = 0; i < N; i++){
			fprintf(otptfl, "%.32e\n", x[i]);
			fprintf(otptfl, "%.32e\n", y[i]);
		}

		fclose(otptfl);
	}

	vector<vector <double>>	Ut(N, vector<double> (N, 0.0));
	vector<vector <int>> 	Touch(N, vector<int> (N, 0));

	// particle properties
	double 			K = 1.0; // spring constant
	double 			Kt = K / 3.0; // tangential spring constant
	double 			g = K / 1000.0; // gravity
	double 			phi_t = 0.5;
	double 			en = 0.9; // restitution
	double 			gamma_v = -4.0 * std::log10(en) * K / PI; // particle collision dissipation
	// fluid properties
	double 			rho_f = 0.5; // density of fluid
	double 			nu = 0.1; // viscosity of fluid
	double 			v0 = (double)v0_c * std::pow(10.0, v0_e);
	double 			B1 = 3.0 * PI * rho_f * nu;
	double 			B2 = 0.005 * PI * rho_f;
	double 			b = 5.0;
	// sediment force threshold
    double 			Fthresh = std::pow(10.0, -12) * K;  // force balance threshold
	double 			Eval, P;
	int 			Nt_flow = 500000;

	if (mu_c == 0){
		std::ifstream pos_sed(sedfile);
		if (pos_sed.good()){
			double Lx_temp;
			pos_sed >> Lx_temp;
			for (i = 0; i < N; i++){
				pos_sed >> x[i];
				pos_sed >> y[i];
			}
		}
		else {
			Sediment(N, Lx, K, g, D, M, x, y, Fthresh);
			for (i = 0; i < N; i++)
				x[i] -= std::floor(x[i] / Lx) * Lx;

			otptfl = fopen(sedfile, "w");
			fprintf(otptfl, "%.32e\n", Lx);
			for (i = 0; i < N; i++){
				fprintf(otptfl, "%.32e\n", x[i]);
				fprintf(otptfl, "%.32e\n", y[i]);
			}
			fclose(otptfl);
		}
		pos_sed.close();
		
		// particles at the bottom are fixed
		// rearrangle particle indices so that bottom particles are at the end of the array
		vector<double> 	D_new(N, 0.0), D_sq_new(N, 0.0), M_new(N, 0.0), I_new(N, 0.0);
		vector<double> 	xp_new(N, 0.0), yp_new(N, 0.0);
		int 			N_bottom = 0, N_up = 0;
		for (i = 0; i < N; i++){
			if (y[i] < 0.5 * D[i]){
				N_bottom++;
				j = N - N_bottom;
				xp_new[j] = x[i];
				yp_new[j] = y[i];
				D_new[j] = D[i];
				D_sq_new[j] = D_sq[i];
				M_new[j] = M[i];
				I_new[j] = I[i];
			}
			else {
				xp_new[N_up] = x[i];
				yp_new[N_up] = y[i];
				D_new[N_up] = D[i];
				D_sq_new[N_up] = D_sq[i];
				M_new[N_up] = M[i];
				I_new[N_up] = I[i];
				N_up++;
			}
		}
		
		Flow(N, N_up, Lx, K, g, Dl, D_new, D_sq_new, M_new, xp_new, yp_new,
            	phi_t, gamma_v, b, B1, B2, v0, Nt_flow, flowdir);
	}
	else{
		std::ifstream pos_sed(sedfile);
		if (pos_sed.good()){
			double Lx_temp;
			pos_sed >> Lx_temp;
			for (i = 0; i < N; i++){
				pos_sed >> x[i];
				pos_sed >> y[i];
			}
			for (i = 0; i < N - 1; i++)
				for (j = i + 1; j < N; j++)
					pos_sed >> Ut[i][j];
			for (i = 0; i < N - 1; i++)
				for (j = i + 1; j < N; j++)
					pos_sed >> Touch[i][j];
		}
		else {
			Sediment(N, Lx, K, Kt, g, D, M, I, x, y, Ut, Touch, mu, Fthresh);
			for (i = 0; i < N; i++)
				x[i] -= std::floor(x[i] / Lx) * Lx;

			otptfl = fopen(sedfile, "w");
			fprintf(otptfl, "%.32e\n", Lx);
			for (i = 0; i < N; i++){
				fprintf(otptfl, "%.32e\n", x[i]);
				fprintf(otptfl, "%.32e\n", y[i]);
			}
			for (i = 0; i < N - 1; i++)
				for (j = i + 1; j < N; j++)
					fprintf(otptfl, "%.32e\n", Ut[i][j]);
			for (i = 0; i < N - 1; i++)
				for (j = i + 1; j < N; j++)
					fprintf(otptfl, "%d\n", Touch[i][j]);
			fclose(otptfl);
		}
		pos_sed.close();

		// particles at the bottom are fixed
		// rearrangle particle indices so that bottom particles are at the end of the array
		vector<double> 	D_new(N, 0.0), D_sq_new(N, 0.0), M_new(N, 0.0), I_new(N, 0.0);
		vector<double> 	xp_new(N, 0.0), yp_new(N, 0.0);
		vector<vector <double>> Ut_new(N, vector<double> (N, 0.0));
		vector<vector <int>> 	Touch_new(N, vector<int> (N, 0));
		vector<int>    	idx_reshuffle(N, 0);
		int 			N_bottom = 0, N_up = 0;
		for (i = 0; i < N; i++){
			if (y[i] < 0.5 * D[i]){
				N_bottom++;
				j = N - N_bottom;
				idx_reshuffle[i] = j;
				D_new[j] = D[i];
				D_sq_new[j] = D_sq[i];
				M_new[j] = M[i];
				I_new[j] = I[i];
			}
			else {
				idx_reshuffle[i] = N_up;
				D_new[N_up] = D[i];
				D_sq_new[N_up] = D_sq[i];
				M_new[N_up] = M[i];
				I_new[N_up] = I[i];
				N_up++;
			}
		}
		for (i = 0; i < N; i++){
			j = idx_reshuffle[i];
			xp_new[j] = x[i];
			yp_new[j] = y[i];
			for (k = i + 1; k < N; k++){
				l = idx_reshuffle[k];
				Ut_new[j][l] = Ut[i][k];
				Ut_new[l][j] = Ut[i][k];
				Touch_new[j][l] = Touch[i][k];
				Touch_new[l][j] = Touch[i][k];
			}
		}
		
		Flow(N, N_up, Lx, K, Kt, g, Dl, D_new, D_sq_new, M_new, I_new, xp_new, yp_new,
            	Ut_new, Touch_new, mu, phi_t, gamma_v, b, B1, B2, v0, Nt_flow, flowdir);
	}

    return 0;

}
