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

#include "Energy.h"
#include "main.h"

using std::vector;

double sign(double x){
    return (x >= 0.0) ? 1.0 : -1.0;
}

double min(double x, double y){
    return (x >= y) ? y : x;
}

double Energy_Disk_VL(int N, double Lx, double K, double g,
                      vector<double> &D, vector<double> &M,
                      vector<double> &x, vector<double> &y,
                      vector<double> &Fx, vector<double> &Fy,
                      vector<vector <int>> &VL, int VL_count)
{
    int 		n, m;
    double 		Eval = 0.0;

	for (n = 0; n < N; n++){
        Fx[n] = 0.0;
        Fy[n] = 0.0;
    }

	// Disk-disk Force
    double 		Dnm, dnm_sq, dnm, dd;
    double 		dx, dy, F, dFx, dFy;
    int 		vl_idx;

	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		n = VL[vl_idx][0];
		m = VL[vl_idx][1];
        Dnm = 0.5 * (D[n] + D[m]);
		dx = x[m] - x[n];
        dx -= std::round(dx / Lx) * Lx;
		if (std::abs(dx) < Dnm){
			dy = y[m] - y[n];
            if (std::abs(dy) < Dnm){
                dnm_sq = dx * dx + dy * dy;
                dnm = std::sqrt(dnm_sq);
                if (dnm < Dnm){
                    dd = Dnm - dnm;
                    Eval += 0.5 * K * dd * dd;

                    F = K * dd / dnm;
                    dFx = F * dx;
                    dFy = F * dy;
                    Fx[n] -= dFx;
                    Fy[n] -= dFy;
                    Fx[m] += dFx; // 3rd law
                    Fy[m] += dFy;
                }
            }
		}
	}

    for (n = 0; n < N; n++){
        Dnm = D[n] / 2.0;
        if (y[n] < Dnm) // force from bottom wall
            Fy[n] += K * (Dnm - y[n]);
        Fy[n] -= M[n] * g; // gravity
    }

	return Eval;
}

void Force_Disk_VL(int N, int N_up, double Lx, double K, double g,
                    double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
                    vector<double> &xp, vector<double> &yp,
                    vector<double> &Vx, vector<double> &Vy,
                    double phi_t, double gamma_v, double b,
                    double B1, double B2, double v0,
                    vector<double> &Fx, vector<double> &Fy,
                    vector<vector <int>> &VL, int VL_count)
{
    int 		n, m;

	for (n = 0; n < N_up; n++){
        Fx[n] = 0.0;
        Fy[n] = 0.0;
    }

	// Disk-disk Force
    vector<double>  phi(N_up, 0.0);
    double      D_phi, D_phi_n, D_phi_m, D_phi_n_sq, D_phi_m_sq;
    double      theta, theta_small, d1, d2;
    double      Dnm, dnm_sq, dnm, dd;
    double 		dx, dy, dFx, dFy, Fn, Fd;
    double      dvx, dvy, dv, dw, vt;
    int 		vl_idx;

    // particle contact & dissipation forces && local packing fraction
	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		n = VL[vl_idx][0];
		m = VL[vl_idx][1];
        Dnm = 0.5 * (D[n] + D[m]);
        D_phi = Dnm + Dl;
        D_phi_n = 0.5 * D[n] + Dl;
        D_phi_m = 0.5 * D[m] + Dl;
        D_phi_n_sq = D_phi_n * D_phi_n;
        D_phi_m_sq = D_phi_m * D_phi_m;
		dx = xp[m] - xp[n];
        dx -= std::round(dx / Lx) * Lx;
        dy = yp[m] - yp[n];
        dnm_sq = dx * dx + dy * dy;
        dnm = std::sqrt(dnm_sq);
        if (dnm < Dnm){ // contact force and local packing fraction
            Fn = K * (Dnm / dnm - 1.0);

            dvx = Vx[n] - Vx[m];
            dvy = Vy[n] - Vy[m];
            Fd = gamma_v / std::sqrt(0.5 * (M[n] + M[m])) * M[n] * M[m] / (M[n] + M[m]) *
                    (dvx * dx + dvy * dy) / dnm_sq; // dissipation force

            dFx = (Fn + Fd) * dx;
            dFy = (Fn + Fd) * dy;
            Fx[n] -= dFx;
            Fy[n] -= dFy;

            phi[n] += 0.25 * PI * D_sq[m];
            
            if (m < N_up){
                Fx[m] += dFx; // 3rd law
                Fy[m] += dFy;
                phi[m] += 0.25 * PI * D_sq[n];
            }
        }
        else if (dnm < D_phi){ // local packing fraction only
            if (dnm < D_phi_n - 0.5 * D[m]) // particle m fully inside of the probe circle
                phi[n] += 0.25 * PI * D_sq[m];
            else if (dnm < std::sqrt(D_phi_n_sq - 0.25 * D_sq[m])){ // particle m partially inside of the probe circle    
                d2 = (D_phi_n_sq - 0.25 * D_sq[m] - dnm_sq) / 2.0 / dnm;
                d1 = std::sqrt(0.25 * D_sq[m] - d2 * d2);
                theta = std::acos(d2 / 0.5 / D[m]);
                theta_small = std::asin(d1 / D_phi_n);
                phi[n] += theta_small * D_phi_n_sq + (d2 - dnm) * d1 +
                            0.25 * (PI - theta) * D_sq[m];
            }
            else {
                d2 = (0.25 * D_sq[m] + dnm_sq - D_phi_n_sq) / 2.0 / dnm;
                d1 = std::sqrt(0.25 * D_sq[m] - d2 * d2);
                theta = std::acos(d2 / 0.5 / D[m]);
                theta_small = std::asin(d1 / D_phi_n);
                phi[n] += theta_small * D_phi_n_sq - dnm * d1 + 0.25 * theta * D_sq[m];
            }
            if (m < N_up){
                if (dnm < D_phi_m - 0.5 * D[n])  // particle n fully inside of the probe circle
                    phi[m] += 0.25 * PI * D_sq[n];
                else if (dnm < std::sqrt(D_phi_m_sq - 0.25 * D_sq[n])){ // particle n partially inside of the probe circle
                    d2 = (D_phi_m_sq - 0.25 * D_sq[n] - dnm_sq) / 2.0 / dnm;
                    d1 = std::sqrt(0.25 * D_sq[n] - d2 * d2);
                    theta = std::acos(d2 / 0.5 / D[n]);
                    theta_small = std::asin(d1 / D_phi_m);
                    phi[m] += theta_small * D_phi_m_sq + (d2 - dnm) * d1 +
                              0.25 * (PI - theta) * D_sq[n];
                }
                else {
                    d2 = (0.25 * D_sq[n] + dnm_sq - D_phi_m_sq) / 2.0 / dnm;
                    d1 = std::sqrt(0.25 * D_sq[n] - d2 * d2);
                    theta = std::acos(d2 / 0.5 / D[n]);
                    theta_small = std::asin(d1 / D_phi_m);
                    phi[m] += theta_small * D_phi_m_sq - dnm * d1 + 0.25 * theta * D_sq[n];
                }
            }
        }
	}

    // fluid force
    double      f_phi;
    double      two_dl = 2.0 * Dl;
    for (n = 0; n < N_up; n++){
        phi[n] += 0.25 * PI * D_sq[n];
        D_phi_n = 0.5 * D[n] + Dl;
        if (yp[n] < D_phi_n){
            theta = std::acos(yp[n] / D_phi_n);
            phi[n] += (theta * D_phi_n - yp[n] * std::sin(theta)) * D_phi_n;
        }
        f_phi = std::exp(-b * (phi[n] / PI / D_phi_n / D_phi_n - phi_t));
        // dvx = v0 * f_phi - Vx[n];
        // dvy = -Vy[n];
        // dv = std::sqrt(dvx * dvx + dvy * dvy);
        // Fx[n] += (B1 * D[n] + B2 * D_sq[n] * dv) * dvx;
        // Fy[n] += (B1 * D[n] + B2 * D_sq[n] * dv) * dvy;
        Fx[n] += B1 * D[n] * (v0 * f_phi - Vx[n]);
        Fy[n] -= B1 * D[n] * Vy[n];
    }

    // gravity
    for (n = 0; n < N_up; n++)
        Fy[n] -= M[n] * g; // gravity

	return;
}


double Energy_Disk_VL(int N, double Lx, double K, double Kt, double g,
                      vector<double> &D, vector<double> &M,
                      vector<double> &x, vector<double> &y,
                      vector<vector <double>> &Ut, vector<vector <int>> &Touch,
                      vector<double> &Vx, vector<double> &Vy, vector<double> &W,
                      double dt, double mu, vector<double> &Fx, vector<double> &Fy,
                      vector<double> &torque, vector<vector <int>> &VL, int VL_count)
{
    int 		n, m;
    double 		Eval = 0.0;
    vector<vector <int>> Touch_temp(N, vector<int> (N, 0));

	for (n = 0; n < N; n++){
        Fx[n] = 0.0;
        Fy[n] = 0.0;
        torque[n] = 0.0;
    }

	// Disk-disk Force
    double 		Dnm, dnm_sq, dnm, dd;
    double 		dx, dy, dFx, dFy, dFnx, dFny, dFtx, dFty, Fn, Ft, T;
    double      dvx, dvy, dw, vt;
    int 		vl_idx;

	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		n = VL[vl_idx][0];
		m = VL[vl_idx][1];
        Dnm = 0.5 * (D[n] + D[m]);
		dx = x[m] - x[n];
        dx -= std::round(dx / Lx) * Lx;
		if (std::abs(dx) < Dnm){
			dy = y[m] - y[n];
            if (std::abs(dy) < Dnm){
                dnm_sq = dx * dx + dy * dy;
                dnm = std::sqrt(dnm_sq);
                if (dnm < Dnm){
                    Touch_temp[n][m] = 1;
                    if (Touch[n][m] == 0)
                        Ut[n][m] = 0.0;

                    dd = Dnm - dnm;
                    Fn = K * dd;

                    dvx = Vx[n] - Vx[m];
                    dvy = Vy[n] - Vy[m];
                    Ut[n][m] = Ut[n][m] + ((dvx * dy - dvy * dx) / dnm - 0.5 * (W[n] * D[n] + W[m] * D[m])) * dt;
                    Ut[n][m] = sign(Ut[n][m]) * min(mu * Fn / Kt, std::abs(Ut[n][m]));

                    Eval += 0.5 * (K * dd * dd + Kt * Ut[n][m] * Ut[n][m]);

                    Fn /= dnm;
                    Ft = -Kt * Ut[n][m] / dnm;
                    dFnx = Fn * dx;
                    dFny = Fn * dy;
                    dFtx = Ft * dy;
                    dFty = -Ft * dx;
                    dFx = dFtx - dFnx;
                    dFy = dFty - dFny;
                    Fx[n] += dFx;
                    Fy[n] += dFy;
                    Fx[m] -= dFx; // 3rd law
                    Fy[m] -= dFy;

                    T = 0.5 * Kt * Ut[n][m];
                    torque[n] += T * D[n];
                    torque[m] += T * D[m];
                }
            }
		}
	}

    for (n = 0; n < N; n++){
        Dnm = D[n] / 2.0;
        if (y[n] < Dnm) // force from bottom wall
            Fy[n] += K * (Dnm - y[n]);
        Fy[n] -= M[n] * g; // gravity
    }

    for (n = 0; n < N; n++)
        for (m = 0; m < N; m++)
            Touch[n][m] = Touch_temp[n][m];

	return Eval;
}


double Energy_Disk(int N, double Lx, double K, double Kt, double g,
                    vector<double> &D, vector<double> &M,
                    vector<double> &x, vector<double> &y,
                    vector<vector <double>> &Ut, vector<vector <int>> &Touch,
                    vector<double> &Vx, vector<double> &Vy, vector<double> &W,
                    double dt, double mu, vector<double> &Fx, vector<double> &Fy,
                    vector<double> &torque)
{
    int 		n, m;
    double 		Eval = 0.0;
    vector<vector <int>> Touch_temp(N, vector<int> (N, 0));

	for (n = 0; n < N; n++){
        Fx[n] = 0.0;
        Fy[n] = 0.0;
        torque[n] = 0.0;
    }

	// Disk-disk Force
    double 		Dnm, dnm_sq, dnm, dd;
    double 		dx, dy, dFx, dFy, dFnx, dFny, dFtx, dFty, Fn, Ft, T;
    double      dvx, dvy, dw, vt;
    int 		vl_idx;

	for (n = 0; n < N - 1; n++){
        for (m = n + 1; m < N; m++){
            Dnm = 0.5 * (D[n] + D[m]);
            dx = x[m] - x[n];
            dx -= std::round(dx / Lx) * Lx;
            if (std::abs(dx) < Dnm){
                dy = y[m] - y[n];
                if (std::abs(dy) < Dnm){
                    dnm_sq = dx * dx + dy * dy;
                    dnm = std::sqrt(dnm_sq);
                    if (dnm < Dnm){
                        Touch_temp[n][m] = 1;
                        if (Touch[n][m] == 0)
                            Ut[n][m] = 0.0;

                        dd = Dnm - dnm;
                        Fn = K * dd;

                        dvx = Vx[n] - Vx[m];
                        dvy = Vy[n] - Vy[m];
                        Ut[n][m] = Ut[n][m] + ((dvx * dy - dvy * dx) / dnm - 0.5 * (W[n] * D[n] + W[m] * D[m])) * dt;
                        Ut[n][m] = sign(Ut[n][m]) * min(mu * Fn / Kt, std::abs(Ut[n][m]));

                        Eval += 0.5 * (K * dd * dd + Kt * Ut[n][m] * Ut[n][m]);

                        Fn /= dnm;
                        Ft = -Kt * Ut[n][m] / dnm;
                        dFnx = Fn * dx;
                        dFny = Fn * dy;
                        dFtx = Ft * dy;
                        dFty = -Ft * dx;
                        dFx = dFtx - dFnx;
                        dFy = dFty - dFny;
                        Fx[n] += dFx;
                        Fy[n] += dFy;
                        Fx[m] -= dFx; // 3rd law
                        Fy[m] -= dFy;

                        T = 0.5 * Kt * Ut[n][m];
                        torque[n] += T * D[n];
                        torque[m] += T * D[m];
                    }
                }
            }
        }
	}

    for (n = 0; n < N; n++){
        Dnm = D[n] / 2.0;
        if (y[n] < Dnm) // force from bottom wall
            Fy[n] += K * (Dnm - y[n]);
        Fy[n] -= M[n] * g; // gravity
    }

    for (n = 0; n < N; n++)
        for (m = 0; m < N; m++)
            Touch[n][m] = Touch_temp[n][m];

	return Eval;
}


void Force_Disk_VL(int N, int N_up, double Lx, double K, double Kt, double g,
                    double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
                    vector<double> &xp, vector<double> &yp,
                    vector<vector <double>> &Ut, vector<vector <int>> &Touch,
                    vector<double> &Vx, vector<double> &Vy, vector<double> &W,
                    double dt, double mu, double phi_t, double gamma_v,
                    double b, double B1, double B2, double v0,
                    vector<double> &Fx, vector<double> &Fy,  vector<double> &torque,
                    vector<vector <int>> &VL, int VL_count)
{
    int 		n, m;
    vector<vector <int>> Touch_temp(N, vector<int> (N, 0));

	for (n = 0; n < N_up; n++){
        Fx[n] = 0.0;
        Fy[n] = 0.0;
        torque[n] = 0.0;
    }

	// Disk-disk Force
    vector<double>  phi(N_up, 0.0);
    double      D_phi, D_phi_n, D_phi_m, D_phi_n_sq, D_phi_m_sq;
    double      theta, theta_small, d1, d2;
    double      Dnm, dnm_sq, dnm, dd;
    double 		dx, dy, dFx, dFy, dFnx, dFny, dFtx, dFty, dFdx, dFdy, Fn, Ft, Fd, T;
    double      dvx, dvy, dv, dw, vt;
    int 		vl_idx;

    // particle contact & dissipation forces && local packing fraction
	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		n = VL[vl_idx][0];
		m = VL[vl_idx][1];
        Dnm = 0.5 * (D[n] + D[m]);
        D_phi = Dnm + Dl;
        D_phi_n = 0.5 * D[n] + Dl;
        D_phi_m = 0.5 * D[m] + Dl;
        D_phi_n_sq = D_phi_n * D_phi_n;
        D_phi_m_sq = D_phi_m * D_phi_m;
		dx = xp[m] - xp[n];
        dx -= std::round(dx / Lx) * Lx;
        dy = yp[m] - yp[n];
        dnm_sq = dx * dx + dy * dy;
        dnm = std::sqrt(dnm_sq);
        if (dnm < Dnm){ // contact force and local packing fraction
            Touch_temp[n][m] = 1;
            if (Touch[n][m] == 0)
                Ut[n][m] = 0.0;

            dd = Dnm - dnm;
            Fn = K * dd;

            dvx = Vx[n] - Vx[m];
            dvy = Vy[n] - Vy[m];
            Ut[n][m] = Ut[n][m] + ((dvx * dy - dvy * dx) / dnm - 0.5 * (W[n] * D[n] + W[m] * D[m])) * dt;
            Ut[n][m] = sign(Ut[n][m]) * min(mu * Fn / Kt, std::abs(Ut[n][m]));

            Fd = gamma_v / std::sqrt(0.5 * (M[n] + M[m])) * M[n] * M[m] / (M[n] + M[m]) *
                    (dvx * dx + dvy * dy) / dnm_sq; // dissipation force
            Fn /= dnm; // normal force
            Ft = -Kt * Ut[n][m] / dnm; // friction

            dFdx = Fd * dx;
            dFdy = Fd * dy;
            dFnx = Fn * dx;
            dFny = Fn * dy;
            dFtx = Ft * dy;
            dFty = -Ft * dx;
            dFx = dFtx - dFnx - dFdx;
            dFy = dFty - dFny - dFdy;
            Fx[n] += dFx;
            Fy[n] += dFy;
            
            T = 0.5 * Kt * Ut[n][m];
            torque[n] += T * D[n];

            phi[n] += 0.25 * PI * D_sq[m];
            
            if (m < N_up){
                Fx[m] -= dFx; // 3rd law
                Fy[m] -= dFy;
                torque[m] += T * D[m];
                phi[m] += 0.25 * PI * D_sq[n];
            }
        }
        else if (dnm < D_phi){ // local packing fraction only
            if (dnm < D_phi_n - 0.5 * D[m]) // particle m fully inside of the probe circle
                phi[n] += 0.25 * PI * D_sq[m];
            else if (dnm < std::sqrt(D_phi_n_sq - 0.25 * D_sq[m])) { // particle m partially inside of the probe circle
                d2 = (D_phi_n_sq - 0.25 * D_sq[m] - dnm_sq) / 2.0 / dnm;
                d1 = std::sqrt(0.25 * D_sq[m] - d2 * d2);
                theta = std::acos(d2 / 0.5 / D[m]);
                theta_small = std::asin(d1 / D_phi_n);
                phi[n] += theta_small * D_phi_n_sq + (d2 - dnm) * d1 +
                            0.25 * (PI - theta) * D_sq[m];
            }
            else {
                d2 = (0.25 * D_sq[m] + dnm_sq - D_phi_n_sq) / 2.0 / dnm;
                d1 = std::sqrt(0.25 * D_sq[m] - d2 * d2);
                theta = std::acos(d2 / 0.5 / D[m]);
                theta_small = std::asin(d1 / D_phi_n);
                phi[n] += theta_small * D_phi_n_sq - dnm * d1 + 0.25 * theta * D_sq[m];
            }
            if (m < N_up){
                if (dnm < D_phi_m - 0.5 * D[n])  // particle n fully inside of the probe circle
                    phi[m] += 0.25 * PI * D_sq[n];
                else if (dnm < std::sqrt(D_phi_m_sq - 0.25 * D_sq[n])){ // particle n partially inside of the probe circle
                    d2 = (D_phi_m_sq - 0.25 * D_sq[n] - dnm_sq) / 2.0 / dnm;
                    d1 = std::sqrt(0.25 * D_sq[n] - d2 * d2);
                    theta = std::acos(d2 / 0.5 / D[n]);
                    theta_small = std::asin(d1 / D_phi_m);
                    phi[m] += theta_small * D_phi_m_sq + (d2 - dnm) * d1 +
                              0.25 * (PI - theta) * D_sq[n];
                }
                else {
                    d2 = (0.25 * D_sq[n] + dnm_sq - D_phi_m_sq) / 2.0 / dnm;
                    d1 = std::sqrt(0.25 * D_sq[n] - d2 * d2);
                    theta = std::acos(d2 / 0.5 / D[n]);
                    theta_small = std::asin(d1 / D_phi_m);
                    phi[m] += theta_small * D_phi_m_sq - dnm * d1 + 0.25 * theta * D_sq[n];
                }
            }
        }
	}

    // fluid force
    double      f_phi;
    for (n = 0; n < N_up; n++){
        phi[n] += 0.25 * PI * D_sq[n];
        D_phi_n = 0.5 * D[n] + Dl;
        if (yp[n] < D_phi_n){
            theta = std::acos(yp[n] / D_phi_n);
            phi[n] += (theta * D_phi_n - yp[n] * std::sin(theta)) * D_phi_n;
        }
        f_phi = std::exp(-b * (phi[n] / PI / D_phi_n / D_phi_n - phi_t));
        // dvx = v0 * f_phi - Vx[n];
        // dvy = -Vy[n];
        // dv = std::sqrt(dvx * dvx + dvy * dvy);
        // Fx[n] += (B1 * D[n] + B2 * D_sq[n] * dv) * dvx;
        // Fy[n] += (B1 * D[n] + B2 * D_sq[n] * dv) * dvy;
        Fx[n] += B1 * D[n] * (v0 * f_phi - Vx[n]);
        Fy[n] -= B1 * D[n] * Vy[n];
    }

    // gravity
    for (n = 0; n < N_up; n++)
        Fy[n] -= M[n] * g; // gravity

    for (n = 0; n < N; n++)
        for (m = 0; m < N; m++)
            Touch[n][m] = Touch_temp[n][m];

	return;
}


void LocalPhi(int N, int N_up, double Lx,
                double Dl, vector<double> &D, vector<double> &D_sq,
                vector<double> &xp, vector<double> &yp)
{
    int 		n, m;

	// Disk-disk Force
    vector<double>  phi(N_up, 0.0);
    double      D_phi = 1.5 * Dl, D_phi_n, D_phi_m;
    double      theta, theta_small, d1, d2;
    double      Dnm, dnm_sq, dnm, dd;
    double 		dx, dy;

    // particle contact & dissipation forces && local packing fraction
	for (n = 0; n < N_up; n++){
		for (m = n + 1; m < N; m++){
            Dnm = 0.5 * (D[n] + D[m]);
            D_phi_n = 0.5 * D[n] + Dl;
            D_phi_m = 0.5 * D[m] + Dl;
            dx = xp[m] - xp[n];
            dx -= std::round(dx / Lx) * Lx;
            dy = yp[m] - yp[n];
            dnm_sq = dx * dx + dy * dy;
            dnm = std::sqrt(dnm_sq);
            if (dnm < Dnm){ // contact force and local packing fraction
                phi[n] += 0.25 * PI * D_sq[m];
                
                if (m < N_up)
                    phi[m] += 0.25 * PI * D_sq[n];
            }
            else if (dnm < D_phi){ // local packing fraction only
                if (dnm < D_phi_n){
                    if (dnm < D_phi_n - 0.5 * D[m]) // particle m fully inside of the probe circle
                        phi[n] += 0.25 * PI * D_sq[m];
                    else { // particle m partially inside of the probe circle
                        d2 = (D_phi_n * D_phi_n - 0.25 * D_sq[m] - dnm_sq) / 2.0 / dnm;
                        d1 = std::sqrt(0.25 * D_sq[m] - d2 * d2);
                        theta = std::acos(d2 / 0.5 / D[m]);
                        theta_small = std::asin(d1 / D_phi_n);
                        phi[n] += theta_small * D_phi_n * D_phi_n + (d2 - dnm) * d1 +
                                0.25 * (PI - theta) * D_sq[m];
                    }
                }
                if ((dnm < D_phi_m) && (m < N_up)){
                    if (dnm < D_phi_m - 0.5 * D[n])  // particle n fully inside of the probe circle
                        phi[m] += 0.25 * PI * D_sq[n];
                    else { // particle n partially inside of the probe circle
                        d2 = (D_phi_m * D_phi_m - 0.25 * D_sq[n] - dnm_sq) / 2.0 / dnm;
                        d1 = std::sqrt(0.25 * D_sq[n] - d2 * d2);
                        theta = std::acos(d2 / 0.5 / D[n]);
                        theta_small = std::asin(d1 / D_phi_m);
                        phi[m] += theta_small * D_phi_m * D_phi_m + (d2 - dnm) * d1 +
                                0.25 * (PI - theta) * D_sq[n];
                    }
                }
            }
        }
	}

    // fluid force
    double      two_dl = 2.0 * Dl;
    for (n = 0; n < N_up; n++){
        phi[n] += 0.25 * PI * D_sq[n];
        D_phi_n = 0.5 * D[n] + Dl;
        if (yp[n] < D_phi_n){
            theta = std::acos(yp[n] / D_phi_n);
            phi[n] += (theta * D_phi_n - yp[n] * std::sin(theta)) * D_phi_n;
        }
    }

    /*
    FILE*   otptfl;
    otptfl = fopen("./LocalPhi.txt", "w");
    for (n = 0; n < N_up; n++)
        fprintf(otptfl, "%.32e\n", phi[n] / PI / (0.5 * D[n] + Dl) / (0.5 * D[n] + Dl));
    fclose(otptfl);
    */

	return;
}