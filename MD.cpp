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

#include "MD.h"
#include "VerletList.h"
#include "Energy.h"
#include "main.h"

using std::vector;

double Sediment(int N, double Lx, double K, double g,
                vector<double> &D, vector<double> &M,
                vector<double> &x, vector<double> &y, double Fthresh)
{
    FILE*       otptfl;
	int 		n, nt = 0;
	vector<double>	Fx(N, 0.0), Fy(N, 0.0);
	vector<double>	Vx(N, 0.0), Vy(N, 0.0);
	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL(10 * N, vector<int> (2, -1));
	int 		VL_count;
	vector<double>	 		x_save(N, 0.0), y_save(N, 0.0);
	double 		r_cut = D[0];

    int         Nt = 100000000;
    int         Ncoll = 20;
    double      dt = PI / std::sqrt(K) / Ncoll; // time step for MD simulation
	double 		dt_half = 0.5 * dt;
	double 		dt_sq_half = 0.5 * dt * dt;
	double 		B = 0.001 * K;
	double 		B_denorm = 1.0 + B * dt_half;

	// Verlet list initialization
	VL_count = VerletList(N, r_cut, Lx, x, y, x_save, y_save, VL, VL_count, first_call);
	// for (n = 0; n < VL_count; n++)
	//  	printf("VL n: %d  VL m: %d\n", VL[n][0], VL[n][1]);
	//printf("VL_count: %5d\n", VL_count);

	double Eval, Acc_max, Acc_abs, Uk;
	Eval = Energy_Disk_VL(N, Lx, K, g, D, M, x, y, Fx, Fy, VL, VL_count);
	Acc_max = std::abs(Fx[0]);
	for (n = 1; n < N; n++){
		Acc_abs = std::abs(Fx[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	for (n = 0; n < N; n++){
		Acc_abs = std::abs(Fy[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	if (Acc_max < Fthresh)
		return Eval;

	first_call = 0;
	for (nt = 1; nt < Nt; nt++){
		// First step in Velocity Verlet
		for (n = 0; n < N; n++){
		 	Vx[n] += dt_half * Fx[n] / M[n];
			Vy[n] += dt_half * Fy[n] / M[n];
		 	x[n] += dt * Vx[n];
			y[n] += dt * Vy[n];
		}
		VL_count = VerletList(N, r_cut, Lx, x, y, x_save, y_save, VL, VL_count, first_call);
		Eval = Energy_Disk_VL(N, Lx, K, g, D, M, x, y, Fx, Fy, VL, VL_count);
        /*
        Uk = 0.0;
        for (n = 0; n < Nall; n++)
            Uk += Vel[n] * Vel[n];
        Uk *= 0.5;
        */
		Acc_max = std::abs(Fx[0]);
		for (n = 1; n < N; n++){
			Acc_abs = std::abs(Fx[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		for (n = 0; n < N; n++){
			Acc_abs = std::abs(Fy[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		if ((Acc_max < Fthresh) || (nt > Nt))
			break;

		// Second step in Velocity Verlet
		for (n = 0; n < N; n++){
		 	Fx[n] = (Fx[n] - B * Vx[n]) / B_denorm;
			Fy[n] = (Fy[n] - B * Vy[n]) / B_denorm;
		 	// F[n] -= Bn * Vel[n];
		 	Vx[n] += dt_half * Fx[n] / M[n];
			Vy[n] += dt_half * Fy[n] / M[n];
		}
        /*
		if (nt % 1000 == 0){
            printf("nt: %d  Max Force: %.4e\n", nt, Acc_max / K);

            otptfl = fopen("./Sed_Intermediate.txt", "a");
            for (n = 0; n < Nall; n++)
                fprintf(otptfl, "%.32e\n", pos[n]);
            fclose(otptfl);
        }
        */
        // if (nt % 10000 == 0) printf("nt: %d  Max Force: %.5e\n", nt, Acc_max / K);
	}
	printf("Iteration Number: %d  Max Force: %.4e  Kinetic Energy: %.4e\n", nt, Acc_max / K, Uk / N);
	return Eval;
}

void Flow(int N, int N_up, double Lx, double K, double g,
            double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
            vector<double> &xp, vector<double> &yp,
            double phi_t, double gamma_v, double b,
            double B1, double B2, double v0, int Nt, char *rundir)
{
    FILE*       otptfl_pos;
    FILE *      otptfl_vel;
    char        posname[300], velname[300];
	int 		n, nt = 0;
	vector<double>	Fx(N_up, 0.0), Fy(N_up, 0.0);
	vector<double>	Vx(N, 0.0), Vy(N, 0.0);
	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL(40 * N_up, vector<int> (2, -1));
	int 		VL_count;
	vector<double>	 		xp_save(N, 0.0), yp_save(N, 0.0);
	double 		r_cut = 2.0 * Dl; // larger cut off distance for local phi calculation

    int         save_step = 1000; // save position and velocity every save_step steps
    int         Ncoll = 20;
    double      dt = PI / std::sqrt(K) / Ncoll; // time step for MD simulation
	double 		dt_half = 0.5 * dt;
	double 		dt_sq_half = 0.5 * dt * dt;

    sprintf(posname, "%s/Pos_Flow_All.txt", rundir);
	sprintf(velname, "%s/Vel_Flow_All.txt", rundir);

    otptfl_pos = fopen(posname, "w");
    for (n = 0; n < N_up; n++){
        fprintf(otptfl_pos, "%.16e\n", xp[n]);
        fprintf(otptfl_pos, "%.16e\n", yp[n]);
    }
    fclose(otptfl_pos);

    
    otptfl_vel = fopen(velname, "w");
    for (n = 0; n < N_up; n++){
        fprintf(otptfl_vel, "%.16e\n", Vx[n]);
        fprintf(otptfl_vel, "%.16e\n", Vy[n]);
    }
    fclose(otptfl_vel);

	// Verlet list initialization
	VL_count = VerletList(N, N_up, r_cut, Lx, xp, yp, xp_save, yp_save, VL, VL_count, first_call);

	double Eval, Uk;
	Force_Disk_VL(N, N_up, Lx, K, g, Dl, D, D_sq, M, xp, yp,
                    Vx, Vy, phi_t, gamma_v, b, B1, B2, v0,
                    Fx, Fy, VL, VL_count);

	first_call = 0;
	for (nt = 1; nt <= Nt; nt++){
		// First step in Velocity Verlet
		for (n = 0; n < N_up; n++){
		 	Vx[n] += dt_half * Fx[n] / M[n];
            Vy[n] += dt_half * Fy[n] / M[n];
		 	xp[n] += dt * Vx[n];
            yp[n] += dt * Vy[n];
		}
		VL_count = VerletList(N, N_up, r_cut, Lx, xp, yp, xp_save, yp_save, VL, VL_count, first_call);
		Force_Disk_VL(N, N_up, Lx, K, g, Dl, D, D_sq, M, xp, yp,
                        Vx, Vy, phi_t, gamma_v, b, B1, B2, v0,
                        Fx, Fy, VL, VL_count);

		// Second step in Velocity Verlet
		for (n = 0; n < N_up; n++){
		 	Vx[n] += dt_half * Fx[n] / M[n];
		    Vy[n] += dt_half * Fy[n] / M[n];
        }
        
        if (nt % save_step == 0){
            otptfl_pos = fopen(posname, "a");
            for (n = 0; n < N_up; n++){
                fprintf(otptfl_pos, "%.16e\n", xp[n]);
                fprintf(otptfl_pos, "%.16e\n", yp[n]);
            }
            fclose(otptfl_pos);

            otptfl_vel = fopen(velname, "a");
            for (n = 0; n < N_up; n++){
                fprintf(otptfl_vel, "%.16e\n", Vx[n]);
                fprintf(otptfl_vel, "%.16e\n", Vy[n]);
            }
            fclose(otptfl_vel);
        }
	}
	return;
}


double Sediment(int N, double Lx, double K, double Kt, double g,
                vector<double> &D, vector<double> &M, vector<double> &I,
                vector<double> &x, vector<double> &y, vector<vector <double>> &Ut,
                vector<vector <int>> &Touch, double mu, double Fthresh)
{
    FILE*       otptfl;
	int 		n, nt = 0;
	vector<double>	Fx(N, 0.0), Fy(N, 0.0);
	vector<double> 	T(N, 0.0);
	vector<double>	Vx(N, 0.0), Vy(N, 0.0);
	vector<double> 	W(N, 0.0);
	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL(10 * N, vector<int> (2, -1));
	int 		VL_count;
	vector<double>	 		x_save(N, 0.0), y_save(N, 0.0);
	double 		r_cut = D[0];

    int         Nt = 100000000;
    int         Ncoll = 20;
    double      dt = PI / std::sqrt(K) / Ncoll; // time step for MD simulation
	double 		dt_half = 0.5 * dt;
	double 		dt_sq_half = 0.5 * dt * dt;
	double 		Bn = 0.001 * K;
	double 		Bt = 0.001 * Kt;
	double 		Bn_denorm = 1.0 + Bn * dt_half;
	double 		Bt_denorm = 1.0 + Bt * dt_half;

	// Verlet list initialization
	VL_count = VerletList(N, r_cut, Lx, x, y, x_save, y_save, VL, VL_count, first_call);
	// for (n = 0; n < VL_count; n++)
	//  	printf("VL n: %d  VL m: %d\n", VL[n][0], VL[n][1]);
	//printf("VL_count: %5d\n", VL_count);

	double Eval, Acc_max, Acc_abs, Uk;
	Eval = Energy_Disk_VL(N, Lx, K, Kt, g, D, M, x, y, Ut, Touch,
                          Vx, Vy, W, dt, mu, Fx, Fy, T, VL, VL_count);
	Acc_max = std::abs(Fx[0]);
	for (n = 1; n < N; n++){
		Acc_abs = std::abs(Fx[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	for (n = 0; n < N; n++){
		Acc_abs = std::abs(Fy[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	for (n = 0; n < N; n++){
		Acc_abs = std::abs(T[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	if (Acc_max < Fthresh)
		return Eval;

	first_call = 0;
	for (nt = 1; nt < Nt; nt++){
		// First step in Velocity Verlet
		for (n = 0; n < N; n++){
			W[n] += dt_half * T[n] / I[n];
		 	Vx[n] += dt_half * Fx[n] / M[n];
			Vy[n] += dt_half * Fy[n] / M[n];
		 	x[n] += dt * Vx[n];
			y[n] += dt * Vy[n];
		}
			
		VL_count = VerletList(N, r_cut, Lx, x, y, x_save, y_save, VL, VL_count, first_call);
		Eval = Energy_Disk_VL(N, Lx, K, Kt, g, D, M, x, y, Ut, Touch,
                              Vx, Vy, W, dt, mu, Fx, Fy, T, VL, VL_count);
        /*
        Uk = 0.0;
        for (n = 0; n < Nall; n++)
            Uk += Vel[n] * Vel[n];
        for (n = 0; n < Nall; n++)
            Uk += W[n] * W[n];
        Uk *= 0.5;
        */
		Acc_max = std::abs(Fx[0]);
		for (n = 1; n < N; n++){
			Acc_abs = std::abs(Fx[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		for (n = 0; n < N; n++){
			Acc_abs = std::abs(Fy[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		for (n = 0; n < N; n++){
			Acc_abs = std::abs(T[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		if ((Acc_max < Fthresh) || (nt > Nt))
			break;

		// Second step in Velocity Verlet
		for (n = 0; n < N; n++){
		 	Fx[n] = (Fx[n] - Bn * Vx[n]) / Bn_denorm;
			Fy[n] = (Fy[n] - Bn * Vy[n]) / Bn_denorm;
		 	// F[n] -= Bn * Vel[n];
		 	Vx[n] += dt_half * Fx[n] / M[n];
			Vy[n] += dt_half * Fy[n] / M[n];
			T[n] = (T[n] - Bt * W[n]) / Bt_denorm;
			// T[n] -= Bt * W[n];
			W[n] += dt_half * T[n] / I[n];
		}
        /*
		if (nt % 1000 == 0){
            printf("nt: %d  Max Force: %.4e\n", nt, Acc_max / K);
            
            otptfl = fopen("./Sed_Intermediate.txt", "a");
            for (n = 0; n < Nall; n++)
                fprintf(otptfl, "%.32e\n", pos[n]);
            fclose(otptfl);
        }
        */
        if (nt % 10000 == 0) printf("nt: %d  Max Force: %.5e\n", nt, Acc_max / K);
	}
	printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
	return Eval;
}


void Flow(int N, int N_up, double Lx, double K, double Kt, double g,
            double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
            vector<double> &I, vector<double> &xp, vector<double> &yp,
            vector<vector <double>> &Ut, vector<vector <int>> &Touch,
            double mu, double phi_t, double gamma_v, double b,
            double B1, double B2, double v0, int Nt, char *rundir)
{
    FILE*       otptfl_pos;
    FILE *      otptfl_vel;
    char        posname[300], velname[300];
	int 		n, nt = 0;
	vector<double>	Fx(N_up, 0.0), Fy(N_up, 0.0), T(N_up, 0.0);
	vector<double>	Vx(N, 0.0), Vy(N, 0.0), W(N, 0.0);
	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL(40 * N_up, vector<int> (2, -1));
	int 		VL_count;
	vector<double>	 		xp_save(N, 0.0), yp_save(N, 0.0);
	double 		r_cut = 2.0 * Dl; // larger cut off distance for local phi calculation

    int         save_step = 1000; // save position and velocity every save_step steps
    int         Ncoll = 20;
    double      dt = PI / std::sqrt(K) / Ncoll; // time step for MD simulation
	double 		dt_half = 0.5 * dt;
	double 		dt_sq_half = 0.5 * dt * dt;

    sprintf(posname, "%s/Pos_Flow_All.txt", rundir);
	sprintf(velname, "%s/Vel_Flow_All.txt", rundir);

    otptfl_pos = fopen(posname, "w");
    for (n = 0; n < N_up; n++){
        fprintf(otptfl_pos, "%.16e\n", xp[n]);
        fprintf(otptfl_pos, "%.16e\n", yp[n]);
    }
    fclose(otptfl_pos);
    
    otptfl_vel = fopen(velname, "w");
    for (n = 0; n < N_up; n++){
        fprintf(otptfl_vel, "%.16e\n", Vx[n]);
        fprintf(otptfl_vel, "%.16e\n", Vy[n]);
    }
    fclose(otptfl_vel);

	// Verlet list initialization
	VL_count = VerletList(N, N_up, r_cut, Lx, xp, yp, xp_save, yp_save, VL, VL_count, first_call);

	double Eval, Uk;
	Force_Disk_VL(N, N_up, Lx, K, Kt, g, Dl, D, D_sq, M, xp, yp, Ut, Touch,
                    Vx, Vy, W, dt, mu, phi_t, gamma_v, b, B1, B2, v0,
                    Fx, Fy, T, VL, VL_count);

	first_call = 0;
	for (nt = 1; nt <= Nt; nt++){
		// First step in Velocity Verlet
		for (n = 0; n < N_up; n++){
		 	Vx[n] += dt_half * Fx[n] / M[n];
            Vy[n] += dt_half * Fy[n] / M[n];
            W[n] += dt_half * T[n] / I[n];
		 	xp[n] += dt * Vx[n];
            yp[n] += dt * Vy[n];
		}
		VL_count = VerletList(N, N_up, r_cut, Lx, xp, yp, xp_save, yp_save, VL, VL_count, first_call);
		Force_Disk_VL(N, N_up, Lx, K, Kt, g, Dl, D, D_sq, M, xp, yp, Ut, Touch,
                        Vx, Vy, W, dt, mu, phi_t, gamma_v, b, B1, B2, v0,
                        Fx, Fy, T, VL, VL_count);

		// Second step in Velocity Verlet
		for (n = 0; n < N_up; n++){
		 	Vx[n] += dt_half * Fx[n] / M[n];
		    Vy[n] += dt_half * Fy[n] / M[n];
			W[n] += dt_half * T[n] / I[n];
        }
        if (nt % save_step == 0){
            otptfl_pos = fopen(posname, "a");
            for (n = 0; n < N_up; n++){
                fprintf(otptfl_pos, "%.16e\n", xp[n]);
                fprintf(otptfl_pos, "%.16e\n", yp[n]);
            }
            fclose(otptfl_pos);

            otptfl_vel = fopen(velname, "a");
            for (n = 0; n < N_up; n++){
                fprintf(otptfl_vel, "%.16e\n", Vx[n]);
                fprintf(otptfl_vel, "%.16e\n", Vy[n]);
            }
            fclose(otptfl_vel);
        }
	}
	return;
}