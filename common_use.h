#pragma once

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

extern const double pi;

extern const char * DATA_DIR;
//extern const char * PROBLEM_PARAMS_FILE;

// try not to use it
extern const double __DOUBLE_ROUNDING_EPS;

typedef struct _particleParams{
    int particles_amount;
    int im_particles_amount;

    double rho_left;
    double press_left;
    double vel_left;
    double energy_left;

    double rho_right;
    double press_right;
    double vel_right;
    double energy_right;

    bool isGas;
} ParticleParams;

typedef struct _problemParams{
    double T;
    double h;
    double tau;
    double t_stop;
    double K;

    double gamma;
    double membrane;
    double left;
    double right;

    bool haveViscosity;
    double alfa;
    double beta;
    double nu_coef;
} ProblemParams;

double spline_kernel(double x_a, double x_b, ProblemParams params);

double spline_gradient(double x_a, double x_b, ProblemParams params);

double found_next_coordinate(double prev_x, double prev_vel, ProblemParams params);

double found_next_rho(int image_amount, double mass, double prev_x, double * prev_image_x, ProblemParams params);

double found_pressure(double rho, double energy, ProblemParams problemParams);

double found_next_energy(double prev_energy, double mass, double prev_vel, double prev_pressure,
                         double * image_prev_pressure, double prev_rho, double * prev_image_rho,
                         int image_amount, double * prev_image_vel,
                         double prev_x, double * prev_image_x, ProblemParams params);

double found_viscosity(double pres_a, double pres_b, double vel_a, double vel_b, double rho_a, double rho_b,
                       double coord_a, double coord_b, ProblemParams problemParams);

double interpolation_value(double looked_coord, double mass, double * function, double * rho, double * coord,
                           int amount, ProblemParams params);

double interpolation_rho(double looked_coord, double mass, double * coord, int amount, ProblemParams params);

double found_mass(double left_lenght, double left_amount, ParticleParams gasParams);

void fill_initial_coord(double * coord, double * image_coord, int real_left_p_num, int real_right_p_num,
                        int image_left_p_num, int image_right_p_num, ProblemParams problemParams);

void fill_initial_gas_massives(double * grho, double * gvel, double * energy, double * pressure,
                               double * image_grho, double * image_gvel, double * image_energy,
                               double * image_pressure, int gas_real_left_p_num, int gas_real_right_p_num,
                               int gas_image_left_p_num, int gas_image_right_p_num, ParticleParams gasParams);

void fill_initial_dust_massives(double * drho, double * dvel, double * image_drho, double * image_dvel,
                                int dust_real_left_p_num, int dust_real_right_p_num,
                                int dust_image_left_p_num, int dust_image_right_p_num, ParticleParams dustParams);

double found_epsilon(double looked_dcoord, double * dcoord, double dmass, double grho, int damount, ProblemParams params);

bool is_close_to_int(double val, double eps);
