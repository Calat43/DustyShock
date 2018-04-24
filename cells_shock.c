#include "cells_shock.h"

void fill_cells_num(int * cells_num, double * x, double cell_length, int cell_amount, ParticleParams particleParams,
                    ProblemParams problemParams)
{

    for (int i = 0; i < particleParams.particles_amount; ++i)
    {
        for(int j = 1; j <= cell_amount; ++j)
        {
            if((x[i] >= problemParams.left + cell_length * (j - 1)) && (x[i] < problemParams.left + cell_length * j))
            {
                cells_num[i] = j - 1;
            }
        }
    }
}

//в первом столбце - количество частиц в ячейке
void fill_cells(int ** particle_cells, int * cell_num, int cell_amount, ParticleParams params)
{
    int j = 0;

    for(int l = 0; l < cell_amount; ++l)
    {
        j = 1;
        for(int i = 0; i < params.particles_amount; ++i)
        {
            if(cell_num[i] == l)
            {
                particle_cells[l][j] = i;
                ++j;
            }
        }
        particle_cells[l][j] = -1;
        particle_cells[l][0] = j - 1;
    }
}

//int cell_num - номер v*
double vel_asterisk(int cell_num, double * velocity, int ** cells)
{
    double result = 0;
    double neighbors = cells[cell_num][0];
    int neighbor = 0;

    if(neighbors == 0)
    {
        return 0;
    }
    else
    {
        for (int k = 1; k <= neighbors; ++k)
        {
            neighbor = cells[cell_num][k];
            result += velocity[neighbor];
        }
        return result / neighbors;
    }
}

double eps_asterisk(int cell_num, double gmass, double dmass, int ** gas_cells, int ** dust_cells)
{
    double gas_neighbors = (double) gas_cells[cell_num][0];
    double dust_neighbors = (double) dust_cells[cell_num][0];

    if(gas_neighbors == 0)
    {
        return 0;
    }
    else
    {
        return (dmass * dust_neighbors) / (gmass * gas_neighbors);
    }
}

double t_stop_asterisk(int cell_num, double * drho, int ** dust_cells, ProblemParams params)
{
    return vel_asterisk(cell_num, drho, dust_cells) / params.K;
}

double a_p(double grho_a, double gcoord_a, double press_a, double gvel_a, int all_gamount, double gmass,
           double * image_grho, double * image_gcoord, double * image_press, double * image_gvel,
           ProblemParams params)
{
    double result = 0;

    for(int j = 0; j < all_gamount; ++j)
    {
        if (fabs(gcoord_a - image_gcoord[j]) < 2.1 * params.h)
        {
            result += (image_press[j] / pow(image_grho[j], 2) + press_a/ pow(grho_a, 2) +
                       found_viscosity(press_a, image_press[j], gvel_a, image_gvel[j], grho_a, image_grho[j], gcoord_a,
                                       image_gcoord[j], params))
                      * spline_gradient(gcoord_a, image_gcoord[j], params);
            assert(!isnan(image_press[j]));
            assert(!isnan(image_grho[j]));
            assert(!isnan(press_a));
            assert(!isnan(grho_a));
            assert(!isnan(found_viscosity(press_a, image_press[j], gvel_a, image_gvel[j], grho_a, image_grho[j], gcoord_a,
                                          image_gcoord[j], params)));
        }
    }

    return - result * gmass;
}

double a_p_asterisk(int ** gas_cells, int cell_num, double * grho, double * gcoord, double * press, double * gvel,
                    int all_gamount, double gmass, double * image_grho, double * image_gcoord, double * image_press,
                    double * image_gvel, ProblemParams params)
{
    double neighbors = gas_cells[cell_num][0];
    int neighbor = 0;
    double result = 0;

    if(neighbors == 0)
    {
        return 0;
    }
    else
    {
        for (int k = 1; k <= neighbors; ++k)
        {
            neighbor = gas_cells[cell_num][k];
            result += a_p(grho[neighbor], gcoord[neighbor], press[neighbor], gvel[neighbor], all_gamount, gmass,
                          image_grho, image_gcoord, image_press, image_gvel, params);
            assert(!isnan(grho[neighbor]));
            assert(!isnan(gcoord[neighbor]));
            assert(!isnan(press[neighbor]));
            assert(!isnan(result));
        }
        assert(neighbors != 0);
        return result / neighbors;
    }
}

static void compute_x(double * x, double * gvel_astr, double * dvel_astr, int cell_amount)
{
    for(int k = 0; k < cell_amount; ++k)
    {
        x[k] = gvel_astr[k] - dvel_astr[k];
    }
}

static void compute_y(double * y, double * gvel_astr, double * dvel_astr, double * eps_astr, int cell_amount)
{
    for(int k = 0; k < cell_amount; ++k)
    {
        y[k] = gvel_astr[k] + eps_astr[k] * dvel_astr[k];
    }
}

static double found_next_x(double prev_x, double a_p_astr, double t_stop_astr, double eps_astr, ProblemParams params)
{
    double denom = 0;

    if(t_stop_astr == 0)
    {
        denom = 1;
    }
    else
    {
        denom = 1 + (eps_astr - 1) * params.tau / t_stop_astr;
    }

    return (params.tau * a_p_astr + prev_x) / denom;
}

static double found_next_y(double prev_y, double a_p_astr, ProblemParams params)
{
    return params.tau * a_p_astr + prev_y;
}

double found_next_gvel_astr(double x, double y, double eps_astr)
{
    return (y + eps_astr * x) / (1 + eps_astr);
}

double found_next_dvel_astr(double x, double y, double eps_astr)
{
    return (y - x) / (1 + eps_astr);
}

double next_gvelocity(double eps_astr, double t_stop_astr, double next_dvel_astr, double grho_a,
                 double gcoord_a, double press_a, double gvel_a, int all_gamount, double gmass, double * image_grho,
                 double * image_gcoord, double * image_press, double * image_gvel, ProblemParams params)
{
    double result = 0;
    double denom = 0;
    double tau = params.tau;
    double ap = a_p(grho_a, gcoord_a, press_a, gvel_a, all_gamount, gmass, image_grho, image_gcoord, image_press,
                    image_gvel, params);

    if(t_stop_astr == 0 || next_dvel_astr == 0)
    {
        result = gvel_a + tau * ap;
    }
    else
    {
        denom = 1. + tau * eps_astr / t_stop_astr;
        result = (gvel_a + tau * ap + tau * eps_astr / t_stop_astr * next_dvel_astr) / denom;
    }

    return result;
}

double next_dvelocity(double prev_dvel, double t_stop_astr, double next_gvel_astr, ProblemParams params)
{
    double result = 0;
    double denom = 0;
    double tau = params.tau;

    if(t_stop_astr == 0 || next_gvel_astr == 0)
    {
        result = prev_dvel;
    }
    else
    {
        denom = 1 + tau / t_stop_astr;
        result = (prev_dvel + tau * next_gvel_astr / t_stop_astr) / denom;
    }

    return result;
}

static double found_momentum(double gmass, double dmass, double * gvel, double * dvel, int re_gamount, int re_damount)
{
	double gasResult = 0;
	double dustResult = 0;
	for(int i = 0; i < re_gamount; ++i)
	{
		gasResult += gvel[i];
	}
	for(int i = 0; i < re_damount; ++i)
	{
		dustResult += dvel[i];
	}
	
	return gmass * gasResult + dmass * dustResult;
}

static double found_gvel_drag(double * epsilon_astr, double * t_stop_astr, double * gvel, double * dvel_astr,
                       int * gas_cells_num, int re_gamount, ProblemParams params)
{
	int cell_num = 0;
	double result = 0;
	for(int i = 0; i < re_gamount; ++i)
	{
		cell_num = gas_cells_num[i];
		result += epsilon_astr[cell_num] * (gvel[i] - dvel_astr[cell_num]) / t_stop_astr[cell_num];
	}
	return result * params.K;
}

static double found_dvel_drag(double * t_stop_astr, double * dvel, double * gvel_astr, int * dust_cells_num, int re_damount,
                       ProblemParams params)
{
	int cell_num = 0;
	double result = 0;
	for(int i = 0; i < re_damount; ++i)
	{
		cell_num = dust_cells_num[i];
		result += (dvel[i] - gvel_astr[cell_num]) / t_stop_astr[cell_num];
	}
	return result * params.K;
}

static double found_gvel_cell_drag(int cell_num, double * epsilon_astr, double * t_stop_astr, double * gvel,
                                   double * dvel_astr, int ** gas_cells, ProblemParams params)
{
    double result = 0;
    double neighbors = gas_cells[cell_num][0];
    int neighbor = 0;

    for(int k = 1; k <= neighbors; ++k)
    {
        neighbor = gas_cells[cell_num][k];
        result += epsilon_astr[cell_num] * (gvel[neighbor] - dvel_astr[cell_num]) / t_stop_astr[cell_num];
    }
    return result * params.K;
}

static double found_dvel_cell_drag(int cell_num, double * t_stop_astr, double * dvel, double * gvel_astr,
                                   int ** dust_cells, ProblemParams params)
{
    double result = 0;
    double neighbors = dust_cells[cell_num][0];
    int neighbor = 0;

    for(int k = 1; k <= neighbors; ++k)
    {
        neighbor = dust_cells[cell_num][k];
        result += (dvel[neighbor] - gvel_astr[cell_num]) / t_stop_astr[cell_num];
    }
    return result * params.K;
}

void cells_shock(double cell_length, ParticleParams gasParams, ParticleParams dustParams, ProblemParams problemParams)
{
    //Блок для газа. BEGIN
    int re_gamount = gasParams.particles_amount;
    int all_gamount = gasParams.particles_amount + gasParams.im_particles_amount;

    double gas_l2r = gasParams.rho_left / gasParams.rho_right;

    // считаем количества реальных частиц
    int gas_real_right_p_num = (int) round((double) gasParams.particles_amount / (gas_l2r + 1));
    int gas_real_left_p_num = gasParams.particles_amount - gas_real_right_p_num;

    // аналогично для виртуальных
    int gas_image_right_p_num = (int) round((double) gasParams.im_particles_amount / (gas_l2r + 1));
    int gas_image_left_p_num = gasParams.im_particles_amount - gas_image_right_p_num;

    double prev_gcoord[re_gamount];
    double image_prev_gcoord[all_gamount];
    double next_gcoord[re_gamount];
    double image_next_gcoord[all_gamount];

    fill_initial_coord(prev_gcoord, image_prev_gcoord, gas_real_left_p_num, gas_real_right_p_num, gas_image_left_p_num,
                       gas_image_right_p_num, problemParams);

    double gas_image_left_lenght = image_prev_gcoord[gas_image_left_p_num + gas_real_left_p_num] - image_prev_gcoord[0];
    double gmass = found_mass(gas_image_left_lenght, gas_image_left_p_num + gas_real_left_p_num, gasParams);

    double prev_grho[re_gamount];
    double next_grho[re_gamount];
    double prev_gvel[re_gamount];
    double next_gvel[re_gamount];
    double prev_energy[re_gamount];
    double next_energy[re_gamount];
    double prev_pressure[re_gamount];
    double next_pressure[re_gamount];

    double image_prev_grho[all_gamount];
    double image_next_grho[all_gamount];
    double image_prev_gvel[all_gamount];
    double image_next_gvel[all_gamount];
    double image_prev_energy[all_gamount];
    double image_next_energy[all_gamount];
    double image_prev_pressure[all_gamount];
    double image_next_pressure[all_gamount];

    fill_initial_gas_massives(prev_grho, prev_gvel, prev_energy, prev_pressure, image_prev_grho, image_prev_gvel,
                              image_prev_energy, image_prev_pressure, gas_real_left_p_num, gas_real_right_p_num,
                              gas_image_left_p_num, gas_image_right_p_num, gasParams);

    for (int i = 0; i < re_gamount; ++i)
    {
        prev_grho[i] = found_next_rho(all_gamount, gmass, prev_gcoord[i], image_prev_gcoord, problemParams);
    }

    for (int j = 0; j < all_gamount; ++j)
    {
        image_prev_grho[j] = found_next_rho(all_gamount, gmass, image_prev_gcoord[j], image_prev_gcoord, problemParams);
    }
    for (int j = 0; j < gas_image_left_p_num; ++j)
    {
        image_prev_grho[j] = prev_grho[0];
    }

    for(int i = 0; i < re_gamount; ++i)
    {
        prev_pressure[i] = found_pressure(prev_grho[i], prev_energy[i], problemParams);
    }
    for(int j = 0; j < all_gamount; ++j)
    {
        image_prev_pressure[j] = found_pressure(image_prev_grho[j], image_prev_energy[j], problemParams);
    }

    for(int j = 0; j < all_gamount; ++j)
    {
        image_next_gcoord[j] = image_prev_gcoord[j];
        image_next_gvel[j] = image_prev_gvel[j];
        image_next_grho[j] = image_prev_grho[j];
        image_next_energy[j] = image_prev_energy[j];
        image_next_pressure[j] = image_prev_pressure[j];
    }
    //Блок для газа. END

    //Блок для пыли. BEGIN
    int re_damount = dustParams.particles_amount;
    int all_damount = dustParams.particles_amount + dustParams.im_particles_amount;

    double dust_l2r = dustParams.rho_left / dustParams.rho_right;

    // считаем количества реальных частиц
    int dust_real_right_p_num = (int) round((double) dustParams.particles_amount / (dust_l2r + 1));
    int dust_real_left_p_num = dustParams.particles_amount - dust_real_right_p_num;

    // аналогично для виртуальных
    int dust_image_right_p_num = (int) round((double) dustParams.im_particles_amount / (dust_l2r + 1));
    int dust_image_left_p_num = dustParams.im_particles_amount - dust_image_right_p_num;

    double prev_dcoord[re_damount];
    double image_prev_dcoord[all_damount];
    double next_dcoord[re_damount];
    double image_next_dcoord[all_damount];

    fill_initial_coord(prev_dcoord, image_prev_dcoord, dust_real_left_p_num, dust_real_right_p_num, dust_image_left_p_num,
                       dust_image_right_p_num, problemParams);

    double dust_image_left_lenght = image_prev_dcoord[dust_image_left_p_num + dust_real_left_p_num] - image_prev_dcoord[0];
    double dmass = found_mass(dust_image_left_lenght, dust_image_left_p_num + dust_real_left_p_num, dustParams);

    double prev_drho[re_damount];
    double next_drho[re_damount];
    double prev_dvel[re_damount];
    double next_dvel[re_damount];

    double image_prev_drho[all_damount];
    double image_next_drho[all_damount];
    double image_prev_dvel[all_damount];
    double image_next_dvel[all_damount];

    fill_initial_dust_massives(prev_drho, prev_dvel, image_prev_drho, image_prev_dvel,
                               dust_real_left_p_num, dust_real_right_p_num,
                               dust_image_left_p_num, dust_image_right_p_num, dustParams);

    for(int i = 0; i < re_damount; ++i)
    {
        prev_drho[i] = found_next_rho(all_damount, dmass, prev_dcoord[i], image_prev_dcoord, problemParams);
    }
    for(int j = 0; j < all_damount; ++j)
    {
        image_prev_drho[j] = found_next_rho(all_damount, dmass, image_prev_dcoord[j], image_prev_dcoord, problemParams);
    }
    for (int j = 0; j < dust_image_left_p_num; ++j)
    {
        image_prev_drho[j] = prev_drho[0];
    }

    for(int j = 0; j < all_damount; ++j)
    {
        image_next_dcoord[j] = image_prev_dcoord[j];
        image_next_dvel[j] = image_prev_dvel[j];
        image_next_drho[j] = image_prev_drho[j];
    }
    //Блок для пыли. END

    int cell_num = 0;

    int cell_amount = (int)((problemParams.right - problemParams.left) / cell_length);
    //re_damount;

    double a_p_astr[cell_amount];
    double eps_astr[cell_amount];
    double t_stop_astr[cell_amount];

    int gas_cells_num[re_gamount];
    //int _gas_cells[cell_amount][re_gamount];
    //int * gas_cells[cell_amount];

    int ** gas_cells = (int**) malloc(cell_amount * sizeof(int*));
    for (int i = 0; i < cell_amount; ++i)
    {
        gas_cells[i] = (int*) malloc(re_gamount * sizeof(int));
    }
    /*
    int ** gas_cells = (int**) malloc(re_gamount * sizeof(int*));
    for (int i = 0; i < re_gamount; ++i)
    {
        gas_cells[i] = (int *) malloc((re_gamount + 1) * sizeof(int));
    }
     */

    int dust_cells_num[re_damount];
    //int _dust_cells[cell_amount][re_damount];
    //int * dust_cells[cell_amount];

    int ** dust_cells = (int**) malloc(cell_amount * sizeof(int*));
    for (int i = 0; i < cell_amount; ++i)
    {
        dust_cells[i] = (int*) malloc(re_damount * sizeof(int));
    }
/*
    int ** dust_cells = (int**) malloc(re_damount * sizeof(int*));
    for (int i = 0; i < re_damount; ++i)
    {
        dust_cells[i] = (int*) malloc((re_damount + 1) * sizeof(int));
    }
*/
    double prev_gvel_astr[cell_amount];
    double next_gvel_astr[cell_amount];

    double prev_dvel_astr[cell_amount];
    double next_dvel_astr[cell_amount];

    double prev_x[cell_amount];
    double next_x[cell_amount];

    double prev_y[cell_amount];
    double next_y[cell_amount];

    for(int k = 0; k < cell_amount; ++k)
    {
        prev_gvel_astr[k] = NAN;
        next_gvel_astr[k] = NAN;
        prev_dvel_astr[k] = NAN;
        next_dvel_astr[k] = NAN;

        prev_x[k] = NAN;
        next_x[k] = NAN;
        prev_y[k] = NAN;
        next_y[k] = NAN;
    }

    double gasDrag = 0;
    double dustDrag = 0;
    double momentum = 0;

    char filename[512];

    sprintf(filename, "%s/im_cellsShock_gas_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * gas0_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(gas0_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord[i], image_prev_grho[i], image_prev_gvel[i],
                image_prev_energy[i], image_prev_pressure[i]);
    }
    fclose(gas0_out);

    sprintf(filename, "%s/im_cellsShock_dust_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * dust0_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(dust0_out, "%lf %lf %lf\n", image_prev_dcoord[i], image_prev_drho[i], image_prev_dvel[i]);
    }
    fclose(dust0_out);

    //sprintf(filename, "%s/cell_momentum_tau%lg_K%lg.dat", DATA_DIR, problemParams.tau, problemParams.K);
    //FILE * moment = fopen(filename, "w");

    /*
    for(int i = 0; i < re_gamount; ++i)
    {
        gas_cells[i][0] = 1;
        gas_cells_num[i] = i;
        gas_cells[i][1] = i;
    }

    for(int i = 0; i < re_damount; ++i)
    {
        dust_cells[i][0] = 1;
        dust_cells_num[i] = i;
        dust_cells[i][1] = i;
    }
     */

    for(int frameId = 0; frameId < floor(problemParams.T / problemParams.tau); ++frameId)
    {
        printf("%d\n", frameId);

        fill_cells_num(gas_cells_num, prev_gcoord, cell_length, cell_amount, gasParams, problemParams);
        fill_cells(gas_cells, gas_cells_num, cell_amount, gasParams);

        fill_cells_num(dust_cells_num, prev_dcoord, cell_length, cell_amount, dustParams, problemParams);
        fill_cells(dust_cells, dust_cells_num, cell_amount, dustParams);

        /*
        for(int k = 0; k < cell_amount; ++k)
        {
            assert(gas_cells[k][0] != 0);
            assert(dust_cells[k][0] != 0);
        }
         */

        for(int k = 0; k < cell_amount; ++k)
        {
            a_p_astr[k] = a_p_asterisk(gas_cells, k, prev_grho, prev_gcoord, prev_pressure, prev_gvel, all_gamount,
                                       gmass, image_prev_grho, image_prev_gcoord, image_prev_pressure, image_prev_gvel,
                                       problemParams);
            eps_astr[k] = eps_asterisk(k, gmass, dmass, gas_cells, dust_cells);

            t_stop_astr[k] = t_stop_asterisk(k, prev_drho, dust_cells, problemParams);

            assert(!isnan(a_p_astr[k]));
            assert(!isnan(eps_astr[k]));
            assert(!isnan(t_stop_astr[k]));
        }

        for(int k = 0; k < cell_amount; ++k)
        {
            prev_gvel_astr[k] = vel_asterisk(k, prev_gvel, gas_cells);
            prev_dvel_astr[k] = vel_asterisk(k, prev_dvel, dust_cells);
            assert(!isnan(prev_gvel_astr[k]));
            assert(!isnan(prev_dvel_astr[k]));
        }

        compute_x(prev_x, prev_gvel_astr, prev_dvel_astr, cell_amount);
        compute_y(prev_y, prev_gvel_astr, prev_dvel_astr, eps_astr, cell_amount);

        for(int k = 0; k < cell_amount; ++k)
        {
            next_x[k] = found_next_x(prev_x[k], a_p_astr[k], t_stop_astr[k], eps_astr[k], problemParams);
            next_y[k] = found_next_y(prev_y[k], a_p_astr[k], problemParams);
        }

        for(int k = 0; k < cell_amount; ++k)
        {
            next_dvel_astr[k] = found_next_dvel_astr(next_x[k], next_y[k], eps_astr[k]);
            next_gvel_astr[k] = found_next_gvel_astr(next_x[k], next_y[k], eps_astr[k]);
        }

        for(int i = 0; i < re_gamount; ++i)
        {
            cell_num = gas_cells_num[i];
            next_gvel[i] = next_gvelocity(eps_astr[cell_num], t_stop_astr[cell_num], next_dvel_astr[cell_num],
                                          prev_grho[i], prev_gcoord[i], prev_pressure[i], prev_gvel[i], all_gamount,
                                          gmass, image_prev_grho, image_prev_gcoord, image_prev_pressure,
                                          image_prev_gvel, problemParams);
            assert(!isnan(next_gvel[i]));
        }

        for(int i = 0; i < re_damount; ++i)
        {
            cell_num = dust_cells_num[i];
            next_dvel[i] = next_dvelocity(prev_dvel[i], t_stop_astr[cell_num], next_gvel_astr[cell_num], problemParams);
            assert(!isnan(next_dvel[i]));
        }

        for(int i = 0; i < re_gamount; ++i)
        {
            if(frameId == 186)
            {
                double x = 8;
            }
            next_gcoord[i] = found_next_coordinate(prev_gcoord[i], prev_gvel[i], problemParams);
            next_energy[i] = found_next_energy(prev_energy[i], gmass, prev_gvel[i], prev_pressure[i],
                                               image_prev_pressure, prev_grho[i], image_prev_grho, all_gamount,
                                               image_prev_gvel, prev_gcoord[i], image_prev_gcoord, problemParams);
            assert(!isnan(next_gcoord[i]));
            assert(!isnan(next_energy[i]));
        }
        for(int i = 0; i < re_damount; ++i)
        {
            next_dcoord[i] = found_next_coordinate(prev_dcoord[i], prev_dvel[i], problemParams);
            assert(!isnan(next_dcoord[i]));
        }

        for(int j = gas_image_left_p_num; j < gas_image_left_p_num + re_gamount; ++j)
        {
            image_next_gcoord[j] = next_gcoord[j - gas_image_left_p_num];
            image_next_gvel[j] = next_gvel[j - gas_image_left_p_num];
            image_next_energy[j] = next_energy[j - gas_image_left_p_num];

            assert(!isnan(image_next_gcoord[j]));
            assert(!isnan(image_next_gvel[j]));
            assert(!isnan(image_next_energy[j]));
        }
        for(int j = dust_image_left_p_num; j < dust_image_left_p_num + re_damount; ++j)
        {
            image_next_dcoord[j] = next_dcoord[j - dust_image_left_p_num];
            image_next_dvel[j] = next_dvel[j - dust_image_left_p_num];

            assert(!isnan(image_next_dcoord[j]));
            assert(!isnan(image_next_dvel[j]));
        }

        for(int i = 0; i < re_gamount; ++i)
        {
            next_grho[i] = found_next_rho(all_gamount, gmass, next_gcoord[i], image_next_gcoord, problemParams);
            assert(!isnan(next_grho[i]));
        }
        for(int j = gas_image_left_p_num; j < gas_image_left_p_num + re_gamount; ++j)
        {
            image_next_grho[j] = next_grho[j - gas_image_left_p_num];
            assert(!isnan(image_next_grho[j]));
        }

        for(int i = 0; i < re_damount; ++i)
        {
            next_drho[i] = found_next_rho(all_damount, dmass, next_dcoord[i], image_next_dcoord, problemParams);
            assert(!isnan(next_drho[i]));
        }
        for(int j = dust_image_left_p_num; j < dust_image_left_p_num + re_damount; ++j)
        {
            image_next_drho[j] = next_drho[j - dust_image_left_p_num];
            assert(!isnan(image_next_drho[j]));
        }

        for(int i = 0; i < re_gamount; ++i)
        {
            next_pressure[i] = found_pressure(next_grho[i], next_energy[i], problemParams);
            assert(!isnan(next_pressure[i]));
        }
        for(int j = gas_image_left_p_num; j < gas_image_left_p_num + re_gamount; ++j)
        {
            image_next_pressure[j] = next_pressure[j - gas_image_left_p_num];
            assert(!isnan(image_next_pressure[j]));
        }

        //momentum = found_momentum(gmass, dmass, prev_gvel, prev_dvel, re_gamount, re_damount);
        //gasDrag = found_gvel_drag(eps_astr, t_stop_astr, prev_gvel, prev_dvel_astr, gas_cells_num, re_gamount, problemParams);
        //dustDrag = found_dvel_drag(t_stop_astr, prev_dvel, prev_gvel_astr, dust_cells_num, re_damount, problemParams);
        //fprintf(moment, "%lf |", momentum);


        /*
        for(int k = 0; k < cell_amount; ++k)
        {
            gasDrag = found_gvel_cell_drag(k, eps_astr, t_stop_astr, prev_gvel, prev_dvel_astr, gas_cells, problemParams);
            dustDrag = found_dvel_cell_drag(k, t_stop_astr, prev_dvel, prev_gvel_astr, dust_cells, problemParams);
            fprintf(moment, "cell_num_%d: %0.10lf %0.10lf %lf \n", k, gasDrag, dustDrag, (gasDrag + dustDrag) / gasDrag * 100.);
        }

        fprintf(moment, "-----\n");
         */

		//fprintf(moment, "%lf %lf %lf %lf\n", momentum, gasDrag, dustDrag, (gasDrag + dustDrag) / gasDrag * 100);

        for(int i = 0; i < re_gamount; ++i)
        {
            prev_gcoord[i] = next_gcoord[i];
            prev_gvel[i] = next_gvel[i];
            prev_grho[i] = next_grho[i];
            prev_energy[i] = next_energy[i];
            prev_pressure[i] = next_pressure[i];
        }

        for(int j = 0; j < all_gamount; ++j)
        {
            image_prev_gcoord[j] = image_next_gcoord[j];
            image_prev_gvel[j] = image_next_gvel[j];
            image_prev_grho[j] = image_next_grho[j];
            image_prev_energy[j] = image_next_energy[j];
            image_prev_pressure[j] = image_next_pressure[j];
        }

        for(int i = 0; i < re_damount; ++i)
        {
            prev_dcoord[i] = next_dcoord[i];
            prev_dvel[i] = next_dvel[i];
            prev_drho[i] = next_drho[i];
        }
        for(int j = 0; j < all_damount; ++j)
        {
            image_prev_dcoord[j] = image_next_dcoord[j];
            image_prev_dvel[j] = image_next_dvel[j];
            image_prev_drho[j] = image_next_drho[j];
        }

        /*
        if(frameId == 211)
        {
            sprintf(filename, "%s/im_cellsShock_gas_T0.023_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
                    problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
                    problemParams.nu_coef, problemParams.K);
            FILE * gas1_out = fopen(filename, "w");
            for (int i = 0; i < all_gamount; i++) {
                fprintf(gas1_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord[i], image_prev_grho[i], image_prev_gvel[i],
                        image_prev_energy[i], image_prev_pressure[i]);
            }
            fclose(gas1_out);

            sprintf(filename, "%s/im_cellsShock_dust_T0.023_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
                    problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
                    problemParams.nu_coef, problemParams.K);
            FILE * dust1_out = fopen(filename, "w");
            for (int i = 0; i < all_gamount; i++) {
                fprintf(dust1_out, "%lf %lf %lf\n", image_prev_dcoord[i], image_prev_drho[i], image_prev_dvel[i]);
            }
            fclose(dust1_out);
        }
         */
    }

    sprintf(filename, "%s/im_cellsShock_gas_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * gas_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(gas_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord[i], image_prev_grho[i], image_prev_gvel[i],
                image_prev_energy[i], image_prev_pressure[i]);
    }
    fclose(gas_out);

    sprintf(filename, "%s/im_cellsShock_dust_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * dust_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(dust_out, "%lf %lf %lf\n", image_prev_dcoord[i], image_prev_drho[i], image_prev_dvel[i]);
    }
    fclose(dust_out);

    //fclose(moment);
}