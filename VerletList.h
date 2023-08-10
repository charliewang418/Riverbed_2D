#ifndef VERLETLIST_H_
#define VERLETLIST_H_

using std::vector;

int VerletList(int N, double r_cut, double Lx,
			   vector<double> &x, vector<double> &y,
			   vector<double> &x_save, vector<double> &y_save,
			   vector<vector <int>> &VL, int VL_count_old, int first_call);

int VerletList(int N, int N_up, double r_cut, double Lx,
			   vector<double> &xp, vector<double> &yp,
			   vector<double> &xp_save, vector<double> &yp_save,
			   vector<vector <int>> &VL, int VL_count_old, int first_call);

#endif /*MAIN_H_*/
