#pragma once

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <vector>

extern const double pi;

extern const char * DATA_DIR;
//extern const char * PROBLEM_PARAMS_FILE;

// try not to use it
extern const double __DOUBLE_ROUNDING_EPS;

typedef unsigned int uint;

typedef struct _particleParams{
    uint particles_amount;
    uint im_particles_amount;

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

double found_next_rho(uint image_amount, double mass, double prev_x, std::vector<double> & prev_image_x,
                      ProblemParams params);

double found_pressure(double rho, double energy, ProblemParams problemParams);

double found_next_energy(double prev_energy, double mass, double prev_vel, double prev_pressure,
                         std::vector<double> & image_prev_pressure, double prev_rho, std::vector<double> & prev_image_rho,
                         uint image_amount, std::vector<double> & prev_image_vel, double prev_x,
                         std::vector<double> & prev_image_x, ProblemParams params);

double found_viscosity(double pres_a, double pres_b, double vel_a, double vel_b, double rho_a, double rho_b,
                       double coord_a, double coord_b, ProblemParams problemParams);

double interpolation_value(double looked_coord, double mass, std::vector<double> & function, std::vector<double> & rho,
                           std::vector<double> & coord, uint amount, ProblemParams params);

double interpolation_rho(double looked_coord, double mass, std::vector<double> & coord, uint amount, ProblemParams params);

double found_mass(double left_lenght, double left_amount, ParticleParams gasParams);

void fill_initial_coord(std::vector<double> & coord, std::vector<double> & image_coord, uint real_left_p_num,
                        uint real_right_p_num, uint image_left_p_num, uint image_right_p_num, ProblemParams problemParams);

void fill_initial_gas_massives(std::vector<double> & grho, std::vector<double> & gvel, std::vector<double> & energy,
                               std::vector<double> & pressure, std::vector<double> & image_grho,
                               std::vector<double> & image_gvel, std::vector<double> & image_energy,
                               std::vector<double> & image_pressure, uint gas_real_left_p_num, uint gas_real_right_p_num,
                               uint gas_image_left_p_num, uint gas_image_right_p_num, ParticleParams gasParams);

void fill_initial_dust_massives(std::vector<double> & drho, std::vector<double> & dvel, std::vector<double> & image_drho,
                                std::vector<double> & image_dvel, uint dust_real_left_p_num, uint dust_real_right_p_num,
                                uint dust_image_left_p_num, uint dust_image_right_p_num, ParticleParams dustParams);

double found_epsilon(double looked_dcoord, std::vector<double> & dcoord, double dmass, double grho, uint damount,
                     ProblemParams params);

bool is_close_to_int(double val, double eps);
