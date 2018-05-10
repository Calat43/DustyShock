#include "cells_shock.hpp"

void fill_cells_num(std::vector<int> & cells_num, std::vector<double> & x, double cell_length, uint cell_amount,
                    ParticleParams particleParams, ProblemParams problemParams)
{

    for (uint i = 0; i < particleParams.particles_amount; ++i)
    {
        for(uint j = 1; j <= cell_amount; ++j)
        {
            if((x.at(i) >= problemParams.left + cell_length * (j - 1)) && (x.at(i) < problemParams.left + cell_length * j))
            {
                cells_num.at(i) = j - 1;
            }
        }
    }
}

//в первом столбце - количество частиц в ячейке
void fill_cells(std::vector<std::vector<int>> & particle_cells, std::vector<int> & cell_num, uint cell_amount,
                ParticleParams params)
{
    int j = 0;

    for(uint l = 0; l < cell_amount; ++l)
    {
        j = 1;
        for(uint i = 0; i < params.particles_amount; ++i)
        {
            if(cell_num.at(i) == (int)l)
            {
                particle_cells.at(l).at(j) = i;
                ++j;
            }
        }
        particle_cells.at(l).at(j) = -1;
        particle_cells.at(l).at(0) = j - 1;
    }
}

//int cell_num - номер v*
double vel_asterisk(uint cell_num, std::vector<double> & velocity, std::vector<std::vector<int>> & cells)
{
    double result = 0;
    double neighbors = cells.at(cell_num).at(0);
    int neighbor = 0;

    if(neighbors == 0)
    {
        return 0;
    }
    else
    {
        for (int k = 1; k <= neighbors; ++k)
        {
            neighbor = cells.at(cell_num).at(k);
            result += velocity.at(neighbor);
        }
        return result / neighbors;
    }
}

double eps_asterisk(uint cell_num, double gmass, double dmass, std::vector<std::vector<int>> & gas_cells,
                    std::vector<std::vector<int>> & dust_cells)
{
    double gas_neighbors = (double) gas_cells.at(cell_num).at(0);
    double dust_neighbors = (double) dust_cells.at(cell_num).at(0);

    if(gas_neighbors == 0)
    {
        return 0;
    }
    else
    {
        return (dmass * dust_neighbors) / (gmass * gas_neighbors);
    }
}

double t_stop_asterisk(uint cell_num, std::vector<double> & drho, std::vector<std::vector<int>> & dust_cells,
                       ProblemParams params)
{
    return vel_asterisk(cell_num, drho, dust_cells) / params.K;
}

double a_p(double grho_a, double gcoord_a, double press_a, double gvel_a, uint all_gamount, double gmass,
           std::vector<double> & image_grho, std::vector<double> & image_gcoord, std::vector<double> & image_press,
           std::vector<double> & image_gvel, ProblemParams params)
{
    double result = 0;

    for(uint j = 0; j < all_gamount; ++j)
    {
        if (fabs(gcoord_a - image_gcoord.at(j)) < 2.1 * params.h)
        {
            result += (image_press.at(j) / pow(image_grho.at(j), 2) + press_a/ pow(grho_a, 2) +
                       found_viscosity(press_a, image_press.at(j), gvel_a, image_gvel.at(j), grho_a, image_grho.at(j), gcoord_a,
                                       image_gcoord.at(j), params))
                      * spline_gradient(gcoord_a, image_gcoord.at(j), params);
            assert(!isnan(image_press.at(j)));
            assert(!isnan(image_grho.at(j)));
            assert(!isnan(press_a));
            assert(!isnan(grho_a));
            assert(!isnan(found_viscosity(press_a, image_press.at(j), gvel_a, image_gvel.at(j), grho_a, image_grho.at(j), gcoord_a,
                                          image_gcoord.at(j), params)));
        }
    }

    return - result * gmass;
}

double a_p_asterisk(std::vector<std::vector<int>> & gas_cells, uint cell_num, std::vector<double> & grho,
                    std::vector<double> & gcoord, std::vector<double> & press, std::vector<double> & gvel,
                    uint all_gamount, double gmass, std::vector<double> & image_grho, std::vector<double> & image_gcoord,
                    std::vector<double> & image_press, std::vector<double> & image_gvel, ProblemParams params)
{
    double neighbors = gas_cells.at(cell_num).at(0);
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
            neighbor = gas_cells.at(cell_num).at(k);
            result += a_p(grho.at(neighbor), gcoord.at(neighbor), press.at(neighbor), gvel.at(neighbor), all_gamount, gmass,
                          image_grho, image_gcoord, image_press, image_gvel, params);
            assert(!isnan(grho.at(neighbor)));
            assert(!isnan(gcoord.at(neighbor)));
            assert(!isnan(press.at(neighbor)));
            assert(!isnan(result));
        }
        assert(neighbors != 0);
        return result / neighbors;
    }
}

static void compute_x(std::vector<double> & x, std::vector<double> & gvel_astr, std::vector<double> & dvel_astr,
                      uint cell_amount)
{
    for(uint k = 0; k < cell_amount; ++k)
    {
        x.at(k) = gvel_astr.at(k) - dvel_astr.at(k);
    }
}

static void compute_y(std::vector<double> & y, std::vector<double> & gvel_astr, std::vector<double> & dvel_astr,
                      std::vector<double> & eps_astr, uint cell_amount)
{
    for(uint k = 0; k < cell_amount; ++k)
    {
        y.at(k) = gvel_astr.at(k) + eps_astr.at(k) * dvel_astr.at(k);
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
        denom = 1 + (eps_astr + 1) * params.tau / t_stop_astr;
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

double found_next_dvel_astr(double x, double y, double eps_astr, int is_drag, double prev_dvel_astr)
{
    double result;
    if(is_drag == 0)
    {
        result = prev_dvel_astr;
    }
    else
    {
        result = (y - x) / (1 + eps_astr);
    }
    return result;
}

double next_gvelocity(double eps_astr, double t_stop_astr, double next_dvel_astr, double grho_a,
                      double gcoord_a, double press_a, double gvel_a, uint all_gamount, double gmass,
                      std::vector<double> & image_grho, std::vector<double> & image_gcoord,
                      std::vector<double> & image_press, std::vector<double> & image_gvel, ProblemParams params)
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

double next_dvelocity(double prev_dvel, double t_stop_astr, double next_gvel_astr, double is_drag,
                      ProblemParams params)
{
    double result = 0;
    double denom = 0;
    double tau = params.tau;

    if(is_drag == 0)
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
/*
static double found_momentum(double gmass, double dmass, std::vector<double> & gvel, std::vector<double> & dvel,
                             uint re_gamount, uint re_damount)
{
	double gasResult = 0;
	double dustResult = 0;
	for(uint i = 0; i < re_gamount; ++i)
	{
		gasResult += gvel.at(i);
	}
	for(uint i = 0; i < re_damount; ++i)
	{
		dustResult += dvel.at(i);
	}
	
	return gmass * gasResult + dmass * dustResult;
}

static double found_gvel_drag(std::vector<double> & epsilon_astr, std::vector<double> & t_stop_astr,
                              std::vector<double> & gvel, std::vector<double> & dvel_astr,
                              std::vector<int> & gas_cells_num, uint re_gamount, ProblemParams params)
{
	int cell_num = 0;
	double result = 0;
	for(uint i = 0; i < re_gamount; ++i)
	{
		cell_num = gas_cells_num.at(i);
		result += epsilon_astr.at(cell_num) * (gvel.at(i) - dvel_astr.at(cell_num)) / t_stop_astr.at(cell_num);
	}
	return result * params.K;
}
 */
/*

static double found_dvel_drag(std::vector<double> & t_stop_astr, std::vector<double> & dvel,
                              std::vector<double> & gvel_astr, std::vector<int> & dust_cells_num, uint re_damount,
                              ProblemParams params)
{
	int cell_num = 0;
	double result = 0;
	for(uint i = 0; i < re_damount; ++i)
	{
		cell_num = dust_cells_num.at(i);
		result += (dvel.at(i) - gvel_astr.at(cell_num)) / t_stop_astr.at(cell_num);
	}
	return result * params.K;
}

static double found_gvel_cell_drag(int cell_num, std::vector<double> & epsilon_astr, std::vector<double> & t_stop_astr,
                                   std::vector<double> & gvel, std::vector<double> & dvel_astr, std::vector<std::vector<int>> & gas_cells,
                                   ProblemParams params)
{
    double result = 0;
    double neighbors = gas_cells.at(cell_num).at(0);
    int neighbor = 0;

    for(int k = 1; k <= neighbors; ++k)
    {
        neighbor = gas_cells.at(cell_num).at(k);
        result += epsilon_astr.at(cell_num) * (gvel.at(neighbor) - dvel_astr.at(cell_num)) / t_stop_astr.at(cell_num);
    }
    return result * params.K;
}

static double found_dvel_cell_drag(int cell_num, std::vector<double> & t_stop_astr, std::vector<double> & dvel,
                                   std::vector<double> & gvel_astr, std::vector<std::vector<int>> & dust_cells, ProblemParams params)
{
    double result = 0;
    double neighbors = dust_cells.at(cell_num).at(0);
    int neighbor = 0;

    for(int k = 1; k <= neighbors; ++k)
    {
        neighbor = dust_cells.at(cell_num).at(k);
        result += (dvel.at(neighbor) - gvel_astr.at(cell_num)) / t_stop_astr.at(cell_num);
    }
    return result * params.K;
}
*/

void cells_shock(double cell_length, ParticleParams gasParams, ParticleParams dustParams, ProblemParams problemParams)
{
    //Блок для газа. BEGIN
    uint re_gamount = gasParams.particles_amount;
    uint all_gamount = gasParams.particles_amount + gasParams.im_particles_amount;

    double gas_l2r = gasParams.rho_left / gasParams.rho_right;

    // считаем количества реальных частиц
    uint gas_real_right_p_num = (uint) round((double) gasParams.particles_amount / (gas_l2r + 1));
    uint gas_real_left_p_num = gasParams.particles_amount - gas_real_right_p_num;

    // аналогично для виртуальных
    uint gas_image_right_p_num = (uint) round((double) gasParams.im_particles_amount / (gas_l2r + 1));
    uint gas_image_left_p_num = gasParams.im_particles_amount - gas_image_right_p_num;

    std::vector<double> prev_gcoord(re_gamount);
    std::vector<double> image_prev_gcoord(all_gamount);
    std::vector<double> next_gcoord(re_gamount);
    std::vector<double> image_next_gcoord(all_gamount);

    fill_initial_coord(prev_gcoord, image_prev_gcoord, gas_real_left_p_num, gas_real_right_p_num, gas_image_left_p_num,
                       gas_image_right_p_num, problemParams);

    double gas_image_left_lenght = image_prev_gcoord.at(gas_image_left_p_num + gas_real_left_p_num) - image_prev_gcoord.at(0);
    double gmass = found_mass(gas_image_left_lenght, gas_image_left_p_num + gas_real_left_p_num, gasParams);

    std::vector<double> prev_grho(re_gamount);
    std::vector<double> next_grho(re_gamount);
    std::vector<double> prev_gvel(re_gamount);
    std::vector<double> next_gvel(re_gamount);
    std::vector<double> prev_energy(re_gamount);
    std::vector<double> next_energy(re_gamount);
    std::vector<double> prev_pressure(re_gamount);
    std::vector<double> next_pressure(re_gamount);

    std::vector<double> image_prev_grho(all_gamount);
    std::vector<double> image_next_grho(all_gamount);
    std::vector<double> image_prev_gvel(all_gamount);
    std::vector<double> image_next_gvel(all_gamount);
    std::vector<double> image_prev_energy(all_gamount);
    std::vector<double> image_next_energy(all_gamount);
    std::vector<double> image_prev_pressure(all_gamount);
    std::vector<double> image_next_pressure(all_gamount);

    std::vector<double> check_gvel(re_gamount);
    std::vector<double> check_image_gvel(all_gamount);

    fill_initial_gas_massives(prev_grho, prev_gvel, prev_energy, prev_pressure, image_prev_grho, image_prev_gvel,
                              image_prev_energy, image_prev_pressure, gas_real_left_p_num, gas_real_right_p_num,
                              gas_image_left_p_num, gas_image_right_p_num, gasParams);

    for (uint i = 0; i < re_gamount; ++i)
    {
        prev_grho.at(i) = found_next_rho(all_gamount, gmass, prev_gcoord.at(i), image_prev_gcoord, problemParams);
    }

    for (uint j = 0; j < all_gamount; ++j)
    {
        image_prev_grho.at(j) = found_next_rho(all_gamount, gmass, image_prev_gcoord.at(j), image_prev_gcoord, problemParams);
    }
    for (uint j = 0; j < gas_image_left_p_num; ++j)
    {
        image_prev_grho.at(j) = prev_grho.at(0);
    }

    for(uint i = 0; i < re_gamount; ++i)
    {
        prev_pressure.at(i) = found_pressure(prev_grho.at(i), prev_energy.at(i), problemParams);
    }
    for(uint j = 0; j < all_gamount; ++j)
    {
        image_prev_pressure.at(j) = found_pressure(image_prev_grho.at(j), image_prev_energy.at(j), problemParams);
    }

    for(uint j = 0; j < all_gamount; ++j)
    {
        image_next_gcoord.at(j) = image_prev_gcoord.at(j);
        image_next_gvel.at(j) = image_prev_gvel.at(j);
        image_next_grho.at(j) = image_prev_grho.at(j);
        image_next_energy.at(j) = image_prev_energy.at(j);
        image_next_pressure.at(j) = image_prev_pressure.at(j);
    }

    //Блок для газа. END

    //Блок для пыли. BEGIN
    uint re_damount = dustParams.particles_amount;
    uint all_damount = dustParams.particles_amount + dustParams.im_particles_amount;

    double dust_l2r = dustParams.rho_left / dustParams.rho_right;

    // считаем количества реальных частиц
    uint dust_real_right_p_num = (uint) round((double) dustParams.particles_amount / (dust_l2r + 1));
    uint dust_real_left_p_num = dustParams.particles_amount - dust_real_right_p_num;

    // аналогично для виртуальных
    uint dust_image_right_p_num = (uint) round((double) dustParams.im_particles_amount / (dust_l2r + 1));
    uint dust_image_left_p_num = dustParams.im_particles_amount - dust_image_right_p_num;

    std::vector<double> prev_dcoord(re_damount);
    std::vector<double> image_prev_dcoord(all_damount);
    std::vector<double> next_dcoord(re_damount);
    std::vector<double> image_next_dcoord(all_damount);

    fill_initial_coord(prev_dcoord, image_prev_dcoord, dust_real_left_p_num, dust_real_right_p_num, dust_image_left_p_num,
                       dust_image_right_p_num, problemParams);

    double dust_image_left_lenght = image_prev_dcoord.at(dust_image_left_p_num + dust_real_left_p_num) - image_prev_dcoord.at(0);
    double dmass = found_mass(dust_image_left_lenght, dust_image_left_p_num + dust_real_left_p_num, dustParams);

    std::vector<double> prev_drho(re_damount);
    std::vector<double> next_drho(re_damount);
    std::vector<double> prev_dvel(re_damount);
    std::vector<double> next_dvel(re_damount);

    std::vector<double> image_prev_drho(all_damount);
    std::vector<double> image_next_drho(all_damount);
    std::vector<double> image_prev_dvel(all_damount);
    std::vector<double> image_next_dvel(all_damount);

    std::vector<double> check_dvel(re_damount);
    std::vector<double> check_image_dvel(all_damount);

    fill_initial_dust_massives(prev_drho, prev_dvel, image_prev_drho, image_prev_dvel,
                               dust_real_left_p_num, dust_real_right_p_num,
                               dust_image_left_p_num, dust_image_right_p_num, dustParams);

    for(uint i = 0; i < re_damount; ++i)
    {
        prev_drho.at(i) = found_next_rho(all_damount, dmass, prev_dcoord.at(i), image_prev_dcoord, problemParams);
    }
    for(uint j = 0; j < all_damount; ++j)
    {
        image_prev_drho.at(j) = found_next_rho(all_damount, dmass, image_prev_dcoord.at(j), image_prev_dcoord, problemParams);
    }
    for (uint j = 0; j < dust_image_left_p_num; ++j)
    {
        image_prev_drho.at(j) = prev_drho.at(0);
    }

    for(uint j = 0; j < all_damount; ++j)
    {
        image_next_dcoord.at(j) = image_prev_dcoord.at(j);
        image_next_dvel.at(j) = image_prev_dvel.at(j);
        image_next_drho.at(j) = image_prev_drho.at(j);
    }
    //Блок для пыли. END

    int cell_num = 0;

    uint cell_amount = (uint)((problemParams.right - problemParams.left) / cell_length);
    //re_damount;

    std::vector<double> a_p_astr(cell_amount);
    std::vector<double> eps_astr(cell_amount);
    std::vector<double> t_stop_astr(cell_amount);

    std::vector<int> gas_cells_num(re_gamount);

    std::vector<std::vector<int>> gas_cells(cell_amount, std::vector<int>(re_gamount));

    /*
    int ** gas_cells = (int**) malloc(cell_amount * sizeof(int*));
    for (int i = 0; i < cell_amount; ++i)
    {
        gas_cells.at(i) = (int*) malloc(re_gamount * sizeof(int));
    }
     */

    std::vector<int> dust_cells_num(re_damount);

    std::vector<std::vector<int>> dust_cells(cell_amount, std::vector<int>(re_damount));

    /*
    int ** dust_cells = (int**) malloc(cell_amount * sizeof(int*));
    for (int i = 0; i < cell_amount; ++i)
    {
        dust_cells.at(i) = (int*) malloc(re_damount * sizeof(int));
    }
    */
    std::vector<double> prev_gvel_astr(cell_amount);
    std::vector<double> next_gvel_astr(cell_amount);

    std::vector<double> prev_dvel_astr(cell_amount);
    std::vector<double> next_dvel_astr(cell_amount);

    std::vector<double> prev_x(cell_amount);
    std::vector<double> next_x(cell_amount);

    std::vector<double> prev_y(cell_amount);
    std::vector<double> next_y(cell_amount);

    std::vector<int> is_drag(cell_amount);

    for(uint k = 0; k < cell_amount; ++k)
    {
        prev_gvel_astr.at(k) = NAN;
        next_gvel_astr.at(k) = NAN;
        prev_dvel_astr.at(k) = NAN;
        next_dvel_astr.at(k) = NAN;

        prev_x.at(k) = NAN;
        next_x.at(k) = NAN;
        prev_y.at(k) = NAN;
        next_y.at(k) = NAN;
    }

    //double gasDrag = 0;
    //double dustDrag = 0;
    //double momentum = 0;

    char filename[512];

    sprintf(filename, "%s/im_cellsShock_gas_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg_d2g%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K, problemParams.d2g);
    FILE * gas0_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(gas0_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord.at(i), image_prev_grho.at(i), image_prev_gvel.at(i),
                image_prev_energy.at(i), image_prev_pressure.at(i));
    }
    fclose(gas0_out);

    sprintf(filename, "%s/im_cellsShock_dust_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg_d2g%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K, problemParams.d2g);
    FILE * dust0_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(dust0_out, "%lf %lf %lf\n", image_prev_dcoord.at(i), image_prev_drho.at(i), image_prev_dvel.at(i));
    }
    fclose(dust0_out);

    sprintf(filename, "%s/gcells_tau%lg_K%lg.dat", DATA_DIR, problemParams.tau, problemParams.K);
    FILE * gcells = fopen(filename, "w");
    sprintf(filename, "%s/dcells_tau%lg_K%lg.dat", DATA_DIR, problemParams.tau, problemParams.K);
    FILE * dcells = fopen(filename, "w");

    /*
    for(int i = 0; i < re_gamount; ++i)
    {
        gas_cells.at(i).at(0) = 1;
        gas_cells_num.at(i) = i;
        gas_cells.at(i).at(1) = i;
    }

    for(int i = 0; i < re_damount; ++i)
    {
        dust_cells.at(i).at(0) = 1;
        dust_cells_num.at(i) = i;
        dust_cells.at(i).at(1) = i;
    }
    */

    for(int frameId = 0; frameId < floor(problemParams.T / problemParams.tau); ++frameId)
    {
        printf("%d\n", frameId);

        fill_cells_num(gas_cells_num, prev_gcoord, cell_length, cell_amount, gasParams, problemParams);
        fill_cells(gas_cells, gas_cells_num, cell_amount, gasParams);

        fill_cells_num(dust_cells_num, prev_dcoord, cell_length, cell_amount, dustParams, problemParams);
        fill_cells(dust_cells, dust_cells_num, cell_amount, dustParams);

        for(uint k = 0; k < cell_amount; ++k)
        {
            if (gas_cells[k][0] == 0 || dust_cells[k][0] == 0)
            {
                is_drag[k] = 0;
            }
            else
            {
                is_drag[k] = 1;
            }
        }

        for(uint k = 0; k < cell_amount; ++k)
        {
            a_p_astr.at(k) = a_p_asterisk(gas_cells, k, prev_grho, prev_gcoord, prev_pressure, prev_gvel, all_gamount,
                                       gmass, image_prev_grho, image_prev_gcoord, image_prev_pressure, image_prev_gvel,
                                       problemParams);
            eps_astr.at(k) = eps_asterisk(k, gmass, dmass, gas_cells, dust_cells);

            t_stop_astr.at(k) = t_stop_asterisk(k, prev_drho, dust_cells, problemParams);

            assert(!isnan(a_p_astr.at(k)));
            assert(!isnan(eps_astr.at(k)));
            assert(!isnan(t_stop_astr.at(k)));
        }

        for(uint k = 0; k < cell_amount; ++k)
        {
            prev_gvel_astr.at(k) = vel_asterisk(k, prev_gvel, gas_cells);
            prev_dvel_astr.at(k) = vel_asterisk(k, prev_dvel, dust_cells);
            assert(!isnan(prev_gvel_astr.at(k)));
            assert(!isnan(prev_dvel_astr.at(k)));
        }

        compute_x(prev_x, prev_gvel_astr, prev_dvel_astr, cell_amount);
        compute_y(prev_y, prev_gvel_astr, prev_dvel_astr, eps_astr, cell_amount);

        for(uint k = 0; k < cell_amount; ++k)
        {
            next_x.at(k) = found_next_x(prev_x.at(k), a_p_astr.at(k), t_stop_astr.at(k), eps_astr.at(k), problemParams);
            next_y.at(k) = found_next_y(prev_y.at(k), a_p_astr.at(k), problemParams);
        }

        for(uint k = 0; k < cell_amount; ++k)
        {
            next_dvel_astr.at(k) = found_next_dvel_astr(next_x.at(k), next_y.at(k), eps_astr.at(k), is_drag.at(k),
                                                        prev_dvel_astr.at(k));
            next_gvel_astr.at(k) = found_next_gvel_astr(next_x.at(k), next_y.at(k), eps_astr.at(k));
        }

        for(uint i = 0; i < re_gamount; ++i)
        {
            cell_num = gas_cells_num.at(i);
            next_gvel.at(i) = next_gvelocity(eps_astr.at(cell_num), t_stop_astr.at(cell_num), next_dvel_astr.at(cell_num),
                                             prev_grho.at(i), prev_gcoord.at(i), prev_pressure.at(i), prev_gvel.at(i), all_gamount,
                                             gmass, image_prev_grho, image_prev_gcoord, image_prev_pressure,
                                             image_prev_gvel, problemParams);
            assert(!isnan(next_gvel.at(i)));
        }

        for(uint i = 0; i < re_damount; ++i)
        {
            cell_num = dust_cells_num.at(i);
            next_dvel.at(i) = next_dvelocity(prev_dvel.at(i), t_stop_astr.at(cell_num), next_gvel_astr.at(cell_num),
                                             is_drag.at(cell_num), problemParams);
            assert(!isnan(next_dvel.at(i)));
        }

        for(uint i = 0; i < re_gamount; ++i)
        {
            next_gcoord.at(i) = found_next_coordinate(prev_gcoord.at(i), prev_gvel.at(i), problemParams);
            next_energy.at(i) = found_next_energy(prev_energy.at(i), gmass, prev_gvel.at(i), prev_pressure.at(i),
                                               image_prev_pressure, prev_grho.at(i), image_prev_grho, all_gamount,
                                               image_prev_gvel, prev_gcoord.at(i), image_prev_gcoord, problemParams);
            assert(!isnan(next_gcoord.at(i)));
            assert(!isnan(next_energy.at(i)));
        }
        for(uint i = 0; i < re_damount; ++i)
        {
            next_dcoord.at(i) = found_next_coordinate(prev_dcoord.at(i), prev_dvel.at(i), problemParams);
            assert(!isnan(next_dcoord.at(i)));
        }

        for(uint j = gas_image_left_p_num; j < gas_image_left_p_num + re_gamount; ++j)
        {
            image_next_gcoord.at(j) = next_gcoord.at(j - gas_image_left_p_num);
            image_next_gvel.at(j) = next_gvel.at(j - gas_image_left_p_num);
            image_next_energy.at(j) = next_energy.at(j - gas_image_left_p_num);

            assert(!isnan(image_next_gcoord.at(j)));
            assert(!isnan(image_next_gvel.at(j)));
            assert(!isnan(image_next_energy.at(j)));
        }

        for(uint j = dust_image_left_p_num; j < dust_image_left_p_num + re_damount; ++j)
        {
            image_next_dcoord.at(j) = next_dcoord.at(j - dust_image_left_p_num);
            image_next_dvel.at(j) = next_dvel.at(j - dust_image_left_p_num);

            assert(!isnan(image_next_dcoord.at(j)));
            assert(!isnan(image_next_dvel.at(j)));
        }

        for(uint i = 0; i < re_gamount; ++i)
        {
            next_grho.at(i) = found_next_rho(all_gamount, gmass, next_gcoord.at(i), image_next_gcoord, problemParams);
            assert(!isnan(next_grho.at(i)));
        }
        for(uint j = gas_image_left_p_num; j < gas_image_left_p_num + re_gamount; ++j)
        {
            image_next_grho.at(j) = next_grho.at(j - gas_image_left_p_num);
            assert(!isnan(image_next_grho.at(j)));
        }

        for(uint i = 0; i < re_damount; ++i)
        {
            next_drho.at(i) = found_next_rho(all_damount, dmass, next_dcoord.at(i), image_next_dcoord, problemParams);
            assert(!isnan(next_drho.at(i)));
        }
        for(uint j = dust_image_left_p_num; j < dust_image_left_p_num + re_damount; ++j)
        {
            image_next_drho.at(j) = next_drho.at(j - dust_image_left_p_num);
            assert(!isnan(image_next_drho.at(j)));
        }

        for(uint i = 0; i < re_gamount; ++i)
        {
            next_pressure.at(i) = found_pressure(next_grho.at(i), next_energy.at(i), problemParams);
            assert(!isnan(next_pressure.at(i)));
        }
        for(uint j = gas_image_left_p_num; j < gas_image_left_p_num + re_gamount; ++j)
        {
            image_next_pressure.at(j) = next_pressure.at(j - gas_image_left_p_num);
            assert(!isnan(image_next_pressure.at(j)));
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

        for(uint i = 0; i < re_gamount; ++i)
        {
            prev_gcoord.at(i) = next_gcoord.at(i);
            prev_gvel.at(i) = next_gvel.at(i);
            prev_grho.at(i) = next_grho.at(i);
            prev_energy.at(i) = next_energy.at(i);
            prev_pressure.at(i) = next_pressure.at(i);
        }

        for(uint j = 0; j < all_gamount; ++j)
        {
            image_prev_gcoord.at(j) = image_next_gcoord.at(j);
            image_prev_gvel.at(j) = image_next_gvel.at(j);
            image_prev_grho.at(j) = image_next_grho.at(j);
            image_prev_energy.at(j) = image_next_energy.at(j);
            image_prev_pressure.at(j) = image_next_pressure.at(j);
        }

        for(uint i = 0; i < re_damount; ++i)
        {
            prev_dcoord.at(i) = next_dcoord.at(i);
            prev_dvel.at(i) = next_dvel.at(i);
            prev_drho.at(i) = next_drho.at(i);
        }
        for(uint j = 0; j < all_damount; ++j)
        {
            image_prev_dcoord.at(j) = image_next_dcoord.at(j);
            image_prev_dvel.at(j) = image_next_dvel.at(j);
            image_prev_drho.at(j) = image_next_drho.at(j);
        }

        /*
        if(frameId == 250 || frameId == 300)
        {
            sprintf(filename, "%s/im_cellsShock_gas_%dframeId_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
                    frameId, problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta,
                    gasParams.particles_amount, problemParams.nu_coef, problemParams.K);
            FILE * gas1_out = fopen(filename, "w");
            for (uint i = 0; i < all_gamount; i++) {
                fprintf(gas1_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord.at(i), image_prev_grho.at(i), image_prev_gvel.at(i),
                        image_prev_energy.at(i), image_prev_pressure.at(i));
            }
            fclose(gas1_out);

            sprintf(filename, "%s/im_cellsShock_dust_%dframeId_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
                    frameId, problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta,
                    dustParams.particles_amount, problemParams.nu_coef, problemParams.K);
            FILE * dust1_out = fopen(filename, "w");
            for (uint i = 0; i < all_gamount; i++) {
                fprintf(dust1_out, "%lf %lf %lf\n", image_prev_dcoord.at(i), image_prev_drho.at(i), image_prev_dvel.at(i));
            }
            fclose(dust1_out);

            for(uint i = 0; i < cell_amount; ++i)
            {
                fprintf(dcells, "%d ", is_drag.at(i));
            }
        }
        */

    }

    sprintf(filename, "%s/im_cellsShock_gas_cell%lg_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg_d2g%lg.dat",
            DATA_DIR, cell_length, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K, problemParams.d2g);
    FILE * gas_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(gas_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord.at(i), image_prev_grho.at(i), image_prev_gvel.at(i),
                image_prev_energy.at(i), image_prev_pressure.at(i));
    }
    fclose(gas_out);

    sprintf(filename, "%s/im_cellsShock_dust_cell%lg_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg_d2g%lg.dat",
            DATA_DIR, cell_length, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K, problemParams.d2g);
    FILE * dust_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(dust_out, "%lf %lf %lf\n", image_prev_dcoord.at(i), image_prev_drho.at(i), image_prev_dvel.at(i));
    }
    fclose(dust_out);

    //fclose(moment);
}
