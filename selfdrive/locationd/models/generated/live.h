/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_932902278235784548);
void inv_err_fun(double *nom_x, double *true_x, double *out_7393045890420779672);
void H_mod_fun(double *state, double *out_4525692026350945785);
void f_fun(double *state, double dt, double *out_4366942196776664177);
void F_fun(double *state, double dt, double *out_8632941670307794939);
void h_3(double *state, double *unused, double *out_452951625100503140);
void H_3(double *state, double *unused, double *out_838983510234956149);
void h_4(double *state, double *unused, double *out_5837194863717773242);
void H_4(double *state, double *unused, double *out_8154432129760471798);
void h_9(double *state, double *unused, double *out_5308196203406902762);
void H_9(double *state, double *unused, double *out_7969806400078223848);
void h_10(double *state, double *unused, double *out_1271034405935205314);
void H_10(double *state, double *unused, double *out_7416623250164249061);
void h_12(double *state, double *unused, double *out_1503866864358181193);
void H_12(double *state, double *unused, double *out_8185756432209870249);
void h_31(double *state, double *unused, double *out_6763410158327123869);
void H_31(double *state, double *unused, double *out_6462488824368921942);
void h_32(double *state, double *unused, double *out_9029394872271426252);
void H_32(double *state, double *unused, double *out_5316710686710810968);
void h_13(double *state, double *unused, double *out_7227059787055038805);
void H_13(double *state, double *unused, double *out_1692690517940362872);
void h_14(double *state, double *unused, double *out_5308196203406902762);
void H_14(double *state, double *unused, double *out_7969806400078223848);
void h_19(double *state, double *unused, double *out_3493869205521582330);
void H_19(double *state, double *unused, double *out_861230834008575461);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);