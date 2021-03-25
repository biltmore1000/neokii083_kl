/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_9071241208105204647);
void inv_err_fun(double *nom_x, double *true_x, double *out_7761142807303768284);
void H_mod_fun(double *state, double *out_4663971624270391086);
void f_fun(double *state, double dt, double *out_6183753959833889470);
void F_fun(double *state, double dt, double *out_3529376992358478119);
void h_25(double *state, double *unused, double *out_2642726364792705572);
void H_25(double *state, double *unused, double *out_2279510880287765963);
void h_24(double *state, double *unused, double *out_3513645333284058760);
void H_24(double *state, double *unused, double *out_412670958672676536);
void h_30(double *state, double *unused, double *out_4128108861557645364);
void H_30(double *state, double *unused, double *out_7334388030661417317);
void h_26(double *state, double *unused, double *out_2960365060481503627);
void H_26(double *state, double *unused, double *out_3585746652242607108);
void h_27(double *state, double *unused, double *out_952278530057024612);
void H_27(double *state, double *unused, double *out_8621970018498042629);
void h_29(double *state, double *unused, double *out_873564406103200281);
void H_29(double *state, double *unused, double *out_4112157511155392445);
void h_28(double *state, double *unused, double *out_2266907590583632232);
void H_28(double *state, double *unused, double *out_6605161519563818330);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
