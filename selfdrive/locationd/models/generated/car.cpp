
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_9071241208105204647) {
   out_9071241208105204647[0] = delta_x[0] + nom_x[0];
   out_9071241208105204647[1] = delta_x[1] + nom_x[1];
   out_9071241208105204647[2] = delta_x[2] + nom_x[2];
   out_9071241208105204647[3] = delta_x[3] + nom_x[3];
   out_9071241208105204647[4] = delta_x[4] + nom_x[4];
   out_9071241208105204647[5] = delta_x[5] + nom_x[5];
   out_9071241208105204647[6] = delta_x[6] + nom_x[6];
   out_9071241208105204647[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7761142807303768284) {
   out_7761142807303768284[0] = -nom_x[0] + true_x[0];
   out_7761142807303768284[1] = -nom_x[1] + true_x[1];
   out_7761142807303768284[2] = -nom_x[2] + true_x[2];
   out_7761142807303768284[3] = -nom_x[3] + true_x[3];
   out_7761142807303768284[4] = -nom_x[4] + true_x[4];
   out_7761142807303768284[5] = -nom_x[5] + true_x[5];
   out_7761142807303768284[6] = -nom_x[6] + true_x[6];
   out_7761142807303768284[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_4663971624270391086) {
   out_4663971624270391086[0] = 1.0;
   out_4663971624270391086[1] = 0.0;
   out_4663971624270391086[2] = 0.0;
   out_4663971624270391086[3] = 0.0;
   out_4663971624270391086[4] = 0.0;
   out_4663971624270391086[5] = 0.0;
   out_4663971624270391086[6] = 0.0;
   out_4663971624270391086[7] = 0.0;
   out_4663971624270391086[8] = 0.0;
   out_4663971624270391086[9] = 1.0;
   out_4663971624270391086[10] = 0.0;
   out_4663971624270391086[11] = 0.0;
   out_4663971624270391086[12] = 0.0;
   out_4663971624270391086[13] = 0.0;
   out_4663971624270391086[14] = 0.0;
   out_4663971624270391086[15] = 0.0;
   out_4663971624270391086[16] = 0.0;
   out_4663971624270391086[17] = 0.0;
   out_4663971624270391086[18] = 1.0;
   out_4663971624270391086[19] = 0.0;
   out_4663971624270391086[20] = 0.0;
   out_4663971624270391086[21] = 0.0;
   out_4663971624270391086[22] = 0.0;
   out_4663971624270391086[23] = 0.0;
   out_4663971624270391086[24] = 0.0;
   out_4663971624270391086[25] = 0.0;
   out_4663971624270391086[26] = 0.0;
   out_4663971624270391086[27] = 1.0;
   out_4663971624270391086[28] = 0.0;
   out_4663971624270391086[29] = 0.0;
   out_4663971624270391086[30] = 0.0;
   out_4663971624270391086[31] = 0.0;
   out_4663971624270391086[32] = 0.0;
   out_4663971624270391086[33] = 0.0;
   out_4663971624270391086[34] = 0.0;
   out_4663971624270391086[35] = 0.0;
   out_4663971624270391086[36] = 1.0;
   out_4663971624270391086[37] = 0.0;
   out_4663971624270391086[38] = 0.0;
   out_4663971624270391086[39] = 0.0;
   out_4663971624270391086[40] = 0.0;
   out_4663971624270391086[41] = 0.0;
   out_4663971624270391086[42] = 0.0;
   out_4663971624270391086[43] = 0.0;
   out_4663971624270391086[44] = 0.0;
   out_4663971624270391086[45] = 1.0;
   out_4663971624270391086[46] = 0.0;
   out_4663971624270391086[47] = 0.0;
   out_4663971624270391086[48] = 0.0;
   out_4663971624270391086[49] = 0.0;
   out_4663971624270391086[50] = 0.0;
   out_4663971624270391086[51] = 0.0;
   out_4663971624270391086[52] = 0.0;
   out_4663971624270391086[53] = 0.0;
   out_4663971624270391086[54] = 1.0;
   out_4663971624270391086[55] = 0.0;
   out_4663971624270391086[56] = 0.0;
   out_4663971624270391086[57] = 0.0;
   out_4663971624270391086[58] = 0.0;
   out_4663971624270391086[59] = 0.0;
   out_4663971624270391086[60] = 0.0;
   out_4663971624270391086[61] = 0.0;
   out_4663971624270391086[62] = 0.0;
   out_4663971624270391086[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6183753959833889470) {
   out_6183753959833889470[0] = state[0];
   out_6183753959833889470[1] = state[1];
   out_6183753959833889470[2] = state[2];
   out_6183753959833889470[3] = state[3];
   out_6183753959833889470[4] = state[4];
   out_6183753959833889470[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6183753959833889470[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6183753959833889470[7] = state[7];
}
void F_fun(double *state, double dt, double *out_3529376992358478119) {
   out_3529376992358478119[0] = 1;
   out_3529376992358478119[1] = 0;
   out_3529376992358478119[2] = 0;
   out_3529376992358478119[3] = 0;
   out_3529376992358478119[4] = 0;
   out_3529376992358478119[5] = 0;
   out_3529376992358478119[6] = 0;
   out_3529376992358478119[7] = 0;
   out_3529376992358478119[8] = 0;
   out_3529376992358478119[9] = 1;
   out_3529376992358478119[10] = 0;
   out_3529376992358478119[11] = 0;
   out_3529376992358478119[12] = 0;
   out_3529376992358478119[13] = 0;
   out_3529376992358478119[14] = 0;
   out_3529376992358478119[15] = 0;
   out_3529376992358478119[16] = 0;
   out_3529376992358478119[17] = 0;
   out_3529376992358478119[18] = 1;
   out_3529376992358478119[19] = 0;
   out_3529376992358478119[20] = 0;
   out_3529376992358478119[21] = 0;
   out_3529376992358478119[22] = 0;
   out_3529376992358478119[23] = 0;
   out_3529376992358478119[24] = 0;
   out_3529376992358478119[25] = 0;
   out_3529376992358478119[26] = 0;
   out_3529376992358478119[27] = 1;
   out_3529376992358478119[28] = 0;
   out_3529376992358478119[29] = 0;
   out_3529376992358478119[30] = 0;
   out_3529376992358478119[31] = 0;
   out_3529376992358478119[32] = 0;
   out_3529376992358478119[33] = 0;
   out_3529376992358478119[34] = 0;
   out_3529376992358478119[35] = 0;
   out_3529376992358478119[36] = 1;
   out_3529376992358478119[37] = 0;
   out_3529376992358478119[38] = 0;
   out_3529376992358478119[39] = 0;
   out_3529376992358478119[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3529376992358478119[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3529376992358478119[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3529376992358478119[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3529376992358478119[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3529376992358478119[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3529376992358478119[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3529376992358478119[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3529376992358478119[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3529376992358478119[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3529376992358478119[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3529376992358478119[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3529376992358478119[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3529376992358478119[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3529376992358478119[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3529376992358478119[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3529376992358478119[56] = 0;
   out_3529376992358478119[57] = 0;
   out_3529376992358478119[58] = 0;
   out_3529376992358478119[59] = 0;
   out_3529376992358478119[60] = 0;
   out_3529376992358478119[61] = 0;
   out_3529376992358478119[62] = 0;
   out_3529376992358478119[63] = 1;
}
void h_25(double *state, double *unused, double *out_2642726364792705572) {
   out_2642726364792705572[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2279510880287765963) {
   out_2279510880287765963[0] = 0;
   out_2279510880287765963[1] = 0;
   out_2279510880287765963[2] = 0;
   out_2279510880287765963[3] = 0;
   out_2279510880287765963[4] = 0;
   out_2279510880287765963[5] = 0;
   out_2279510880287765963[6] = 1;
   out_2279510880287765963[7] = 0;
}
void h_24(double *state, double *unused, double *out_3513645333284058760) {
   out_3513645333284058760[0] = state[4];
   out_3513645333284058760[1] = state[5];
}
void H_24(double *state, double *unused, double *out_412670958672676536) {
   out_412670958672676536[0] = 0;
   out_412670958672676536[1] = 0;
   out_412670958672676536[2] = 0;
   out_412670958672676536[3] = 0;
   out_412670958672676536[4] = 1;
   out_412670958672676536[5] = 0;
   out_412670958672676536[6] = 0;
   out_412670958672676536[7] = 0;
   out_412670958672676536[8] = 0;
   out_412670958672676536[9] = 0;
   out_412670958672676536[10] = 0;
   out_412670958672676536[11] = 0;
   out_412670958672676536[12] = 0;
   out_412670958672676536[13] = 1;
   out_412670958672676536[14] = 0;
   out_412670958672676536[15] = 0;
}
void h_30(double *state, double *unused, double *out_4128108861557645364) {
   out_4128108861557645364[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7334388030661417317) {
   out_7334388030661417317[0] = 0;
   out_7334388030661417317[1] = 0;
   out_7334388030661417317[2] = 0;
   out_7334388030661417317[3] = 0;
   out_7334388030661417317[4] = 1;
   out_7334388030661417317[5] = 0;
   out_7334388030661417317[6] = 0;
   out_7334388030661417317[7] = 0;
}
void h_26(double *state, double *unused, double *out_2960365060481503627) {
   out_2960365060481503627[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3585746652242607108) {
   out_3585746652242607108[0] = 0;
   out_3585746652242607108[1] = 0;
   out_3585746652242607108[2] = 0;
   out_3585746652242607108[3] = 0;
   out_3585746652242607108[4] = 0;
   out_3585746652242607108[5] = 0;
   out_3585746652242607108[6] = 0;
   out_3585746652242607108[7] = 1;
}
void h_27(double *state, double *unused, double *out_952278530057024612) {
   out_952278530057024612[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8621970018498042629) {
   out_8621970018498042629[0] = 0;
   out_8621970018498042629[1] = 0;
   out_8621970018498042629[2] = 0;
   out_8621970018498042629[3] = 1;
   out_8621970018498042629[4] = 0;
   out_8621970018498042629[5] = 0;
   out_8621970018498042629[6] = 0;
   out_8621970018498042629[7] = 0;
}
void h_29(double *state, double *unused, double *out_873564406103200281) {
   out_873564406103200281[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4112157511155392445) {
   out_4112157511155392445[0] = 0;
   out_4112157511155392445[1] = 1;
   out_4112157511155392445[2] = 0;
   out_4112157511155392445[3] = 0;
   out_4112157511155392445[4] = 0;
   out_4112157511155392445[5] = 0;
   out_4112157511155392445[6] = 0;
   out_4112157511155392445[7] = 0;
}
void h_28(double *state, double *unused, double *out_2266907590583632232) {
   out_2266907590583632232[0] = state[5];
   out_2266907590583632232[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6605161519563818330) {
   out_6605161519563818330[0] = 0;
   out_6605161519563818330[1] = 0;
   out_6605161519563818330[2] = 0;
   out_6605161519563818330[3] = 0;
   out_6605161519563818330[4] = 0;
   out_6605161519563818330[5] = 1;
   out_6605161519563818330[6] = 0;
   out_6605161519563818330[7] = 0;
   out_6605161519563818330[8] = 0;
   out_6605161519563818330[9] = 0;
   out_6605161519563818330[10] = 0;
   out_6605161519563818330[11] = 0;
   out_6605161519563818330[12] = 0;
   out_6605161519563818330[13] = 0;
   out_6605161519563818330[14] = 1;
   out_6605161519563818330[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
