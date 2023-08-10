#ifndef MD_H_
#define MD_H_

using std::vector;

double Sediment(int N, double Lx, double K, double g,
                vector<double> &D, vector<double> &M,
                vector<double> &x, vector<double> &y, double Fthresh);

void Flow(int N, int N_up, double Lx, double K, double g,
            double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
            vector<double> &xp, vector<double> &yp,
            double phi_t, double gamma_v, double b,
            double B1, double B2, double v0, int Nt, char *rundir);

double Sediment(int N, double Lx, double K, double Kt, double g,
                vector<double> &D, vector<double> &M, vector<double> &I,
                vector<double> &x, vector<double> &y, vector<vector <double>> &Ut,
                vector<vector <int>> &Touch, double mu, double Fthresh);

void Flow(int N, int N_up, double Lx, double K, double Kt, double g,
            double Dl, vector<double> &D, vector<double> &D_sq, vector<double> &M,
            vector<double> &I, vector<double> &xp, vector<double> &yp,
            vector<vector <double>> &Ut, vector<vector <int>> &Touch,
            double mu, double phi_t, double gamma_v, double b,
            double B1, double B2, double v0, int Nt, char *rundir);

#endif /*MAIN_H_*/