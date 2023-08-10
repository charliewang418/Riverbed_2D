#ifndef ENERGY_H_
#define ENERGY_H_

using std::vector;

double sign(double x);

double min(double x, double y);

double Energy_Disk_VL(int N, double Lx, double K, double g,
                      vector<double> &D, vector<double> &M,
                      vector<double> &x, vector<double> &y,
                      vector<double> &Fx, vector<double> &Fy,
                      vector<vector <int>> &VL, int VL_count);

void Force_Disk_VL(int N, int N_up, double Lx, double K, double g,
                    double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
                    vector<double> &xp, vector<double> &yp,
                    vector<double> &Vx, vector<double> &Vy,
                    double phi_t, double gamma_v, double b,
                    double B1, double B2, double v0,
                    vector<double> &Fx, vector<double> &Fy,
                    vector<vector <int>> &VL, int VL_count);

double Energy_Disk_VL(int N, double Lx, double K, double Kt, double g,
                      vector<double> &D, vector<double> &M,
                      vector<double> &x, vector<double> &y,
                      vector<vector <double>> &Ut, vector<vector <int>> &Touch,
                      vector<double> &Vx, vector<double> &Vy, vector<double> &W,
                      double dt, double mu, vector<double> &Fx, vector<double> &Fy,
                      vector<double> &torque, vector<vector <int>> &VL, int VL_count);

double Energy_Disk(int N, double Lx, double K, double Kt, double g,
                    vector<double> &D, vector<double> &M,
                    vector<double> &x, vector<double> &y,
                    vector<vector <double>> &Ut, vector<vector <int>> &Touch,
                    vector<double> &Vx, vector<double> &Vy, vector<double> &W,
                    double dt, double mu, vector<double> &Fx, vector<double> &Fy,
                    vector<double> &torque);

void Force_Disk_VL(int N, int N_up, double Lx, double K, double Kt, double g,
                    double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
                    vector<double> &xp, vector<double> &yp,
                    vector<vector <double>> &Ut, vector<vector <int>> &Touch,
                    vector<double> &Vx, vector<double> &Vy, vector<double> &W,
                    double dt, double mu, double phi_t, double gamma_v,
                    double b, double B1, double B2, double v0,
                    vector<double> &Fx, vector<double> &Fy,  vector<double> &torque,
                    vector<vector <int>> &VL, int VL_count);

void LocalPhi(int N, int N_up, double Lx,
                double Dl, vector<double> &D, vector<double> &D_sq,
                vector<double> &xp, vector<double> &yp);

#endif /*MAIN_H_*/
