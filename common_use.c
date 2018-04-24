#include "common_use.h"

const double pi = 3.14159265358;

const char * DATA_DIR = "/home/calat/CLionProjects/DustyShock";

//const char * PROBLEM_PARAMS_FILE = "/home/calat/CLionProjects/DustyShock/...";

const double __DOUBLE_ROUNDING_EPS = 0.0000001;

double spline_kernel(double x_a, double x_b, ProblemParams params)
{
    double h = params.h;
    double r = fabs(x_a - x_b);
    double q = r / h;
    double result = 0;

    if (r / h >= 0 && r / h <= 1)
    {
        result = 1 - 3. / 2 * pow(q, 2) + 3. / 4 * pow(q, 3);
        return 2./ 3. / h * result;
    }
    if (r / h > 1 && r / h <= 2)
    {
        result = 1. / 4 * pow((2. - q), 3);
        return 2./ 3. / h * result;
    }
    if(r / h > 2)
    {
        return 0;
    }
    assert(false);
}

double spline_gradient(double x_a, double x_b, ProblemParams params)
{
    double h = params.h;
    double r = fabs(x_a - x_b);
    double q = r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }
    if (x_a > x_b)
    {
        return 2./3. / h / h * result;
    }
    if (x_a == x_b)
    {
        return 0;
    }
    if (x_a < x_b)
    {
        return - 2./3. / h / h * result;
    }

    assert(false);
}

double found_next_coordinate(double prev_x, double prev_vel, ProblemParams params)
{
    return prev_x + params.tau * prev_vel;
}

double found_next_rho(int image_amount, double mass, double prev_x, double * prev_image_x, ProblemParams params)
{
    double rho = 0;
    for(int j = 0; j < image_amount; ++j)
    {
        if (fabs(prev_x - prev_image_x[j]) < 2.1 * params.h)
        {
            rho += spline_kernel(prev_x, prev_image_x[j], params);
        }
    }
    assert(rho > 0);
    return mass * rho;
}

double found_pressure(double rho, double energy, ProblemParams problemParams)
{
    return rho * energy * (problemParams.gamma - 1);
}

//for gas
double found_next_energy(double prev_energy, double mass, double prev_vel, double prev_pressure,
                         double * image_prev_pressure, double prev_rho, double * prev_image_rho,
                         int image_amount, double * prev_image_vel,
                         double prev_x, double * prev_image_x, ProblemParams params)
{
    double pres_member = 0;
    double visc_member = 0;
    double result = 0;

    for(int j = 0; j < image_amount; ++j)
    {
        if (fabs(prev_x - prev_image_x[j]) < 2.1 * params.h)
        {
            pres_member += (prev_vel - prev_image_vel[j]) * spline_gradient(prev_x, prev_image_x[j], params);
            visc_member += (prev_vel - prev_image_vel[j]) * spline_gradient(prev_x, prev_image_x[j], params) *
                           found_viscosity(prev_pressure, image_prev_pressure[j], prev_vel, prev_image_vel[j],
                                           prev_rho, prev_image_rho[j], prev_x, prev_image_x[j], params);
        }
    }

    result = prev_pressure /  pow(prev_rho, 2) * mass * params.tau * pres_member +
             params.tau * mass / 2. * visc_member + prev_energy;
    assert(visc_member >= 0);
    assert(result >= 0);

    return result;
}

bool is_close_to_int(double val, double eps) {
    return (floor( val + 0.5 ) - eps <= val) && (val <= floor( val + 0.5 ) + eps);
}

int find_nearest_int_that_divides(int nearest_to, int divides_by) {
    return divides_by * (int) round((double) nearest_to / (double) divides_by);
}

double found_sound_speed(double pressure, double rho, double gamma)
{
    assert(!isnan(sqrt(gamma * pressure / rho)));
    return sqrt(gamma * pressure / rho);
}

double found_mu(double vel_ab, double coord_ab, ProblemParams problemParams)
{
    return (problemParams.h * vel_ab * coord_ab) /
            (pow(coord_ab, 2) + pow(problemParams.nu_coef * problemParams.h, 2));
}

double found_viscosity(double pres_a, double pres_b, double vel_a, double vel_b, double rho_a, double rho_b,
                       double coord_a, double coord_b, ProblemParams problemParams)
{
    if(problemParams.haveViscosity == false)
    {
        return 0;
    }
    else
    {
        double vel_ab = vel_a - vel_b;
        double coord_ab = coord_a - coord_b;

        if(vel_ab * coord_ab >= 0)
        {
            return 0;
        }
        if(vel_ab * coord_ab < 0)
        {
            double c_a = found_sound_speed(pres_a, rho_a, problemParams.gamma);
            double c_b = found_sound_speed(pres_b, rho_b, problemParams.gamma);

            double rho_ab = 1. / 2 * (rho_a + rho_b);
            double c_ab = 1. / 2 * (c_a + c_b);

            double mu_ab = found_mu(vel_ab, coord_ab, problemParams);

            assert(!isnan(c_a));
            assert(!isnan(c_b));
            assert(!isnan(mu_ab));

            return (- problemParams.alfa * c_ab * mu_ab + problemParams.beta * pow(mu_ab, 2)) / rho_ab;
        }
        assert(false);
    }
}

double interpolation_value(double looked_coord, double mass, double * function, double * rho, double * coord,
                           int amount, ProblemParams params)
{
    double result = 0;

    for(int i = 0; i < amount; ++i)
    {
        if (fabs(looked_coord - coord[i]) < 2.1 * params.h)
        {
            result += function[i] / rho[i] * spline_kernel(looked_coord, coord[i], params);
        }
    }
    return mass * result;
}

double interpolation_rho(double looked_coord, double mass, double * coord, int amount, ProblemParams params)
{
    double result = 0;

    for(int i = 0; i < amount; ++i)
    {
        if(fabs(looked_coord - coord[i]) < 2.1 * params.h)
        {
            result += spline_kernel(looked_coord, coord[i], params);
        }
    }
    return mass * result;
}

double found_epsilon(double looked_dcoord, double * dcoord, double dmass, double grho, int damount, ProblemParams params)
{
    return interpolation_rho(looked_dcoord, dmass, dcoord, damount, params) / grho;
}

double found_mass(double left_lenght, double left_amount, ParticleParams params)
{
    return left_lenght * params.rho_left / left_amount;
}

void fill_initial_coord(double * coord, double * image_coord, int real_left_p_num, int real_right_p_num,
                        int image_left_p_num, int image_right_p_num, ProblemParams problemParams)
{
    // шаг = длина / кол-во
    double step_left = (problemParams.membrane - problemParams.left) / real_left_p_num;
    double step_right = (problemParams.right - problemParams.membrane) / real_right_p_num;

    // заполняем левые реальные в оба массива
    coord[0] = problemParams.left + (step_left / 2);
    image_coord[image_left_p_num] = coord[0];
    for (int i = 1; i < real_left_p_num; i++) {
        coord[i] = coord[i - 1] + step_left;
        image_coord[image_left_p_num + i] = coord[i];
    }
    // досчитываем отдельно левые мнимые в image массив
    for (int i = image_left_p_num - 1; i >= 0; i--) {
        image_coord[i] = image_coord[i + 1] - step_left;
    }

    // пишем правые реальные в оба масива сразу
    coord[real_left_p_num] = problemParams.membrane + (step_right / 2);
    image_coord[image_left_p_num + real_left_p_num] = coord[real_left_p_num];
    for (int i = 1; i < real_right_p_num; i++) {
        coord[real_left_p_num + i] = coord[real_left_p_num + i - 1] + step_right;
        image_coord[image_left_p_num + real_left_p_num + i] = coord[real_left_p_num + i];
    }
    // отдельно правые мнимые в большой массив
    for (int i = 0; i < image_right_p_num; i++) {
        int r_image_id = image_left_p_num + real_left_p_num + real_right_p_num + i;
        image_coord[r_image_id] = image_coord[r_image_id - 1] + step_right;
    }
}

void fill_initial_gas_massives(double * grho, double * gvel, double * energy, double * pressure,
                               double * image_grho, double * image_gvel, double * image_energy,
                               double * image_pressure, int gas_real_left_p_num, int gas_real_right_p_num,
                               int gas_image_left_p_num, int gas_image_right_p_num, ParticleParams gasParams)
{
    int all_left = gas_real_left_p_num + gas_image_left_p_num;
    int all_right = gas_real_right_p_num + gas_image_right_p_num;

    for (int j = 0; j < all_left; ++j)
    {
        image_grho[j] = gasParams.rho_left;
        image_gvel[j] = gasParams.vel_left;
        image_energy[j] = gasParams.energy_left;
        image_pressure[j] = gasParams.press_left;
    }
    for (int j = all_left; j < all_left + all_right; ++j)
    {
        image_grho[j] = gasParams.rho_right;
        image_gvel[j] = gasParams.vel_right;
        image_energy[j] = gasParams.energy_right;
        image_pressure[j] = gasParams.press_right;
    }

    for (int i = 0; i < gas_real_left_p_num; ++i)
    {
        grho[i] = gasParams.rho_left;
        gvel[i] = gasParams.vel_left;
        energy[i] = gasParams.energy_left;
        pressure[i] = gasParams.press_left;
    }
    for(int i = gas_real_left_p_num; i < gas_real_left_p_num + gas_real_right_p_num; ++i)
    {
        grho[i] = gasParams.rho_right;
        gvel[i] = gasParams.vel_right;
        energy[i] = gasParams.energy_right;
        pressure[i] = gasParams.press_right;
    }
}

void fill_initial_dust_massives(double * drho, double * dvel, double * image_drho, double * image_dvel,
                                int dust_real_left_p_num, int dust_real_right_p_num,
                                int dust_image_left_p_num, int dust_image_right_p_num, ParticleParams dustParams)
{
    int all_left = dust_real_left_p_num + dust_image_left_p_num;
    int all_right = dust_real_right_p_num + dust_image_right_p_num;

    for (int j = 0; j < all_left; ++j)
    {
        image_drho[j] = dustParams.rho_left;
        image_dvel[j] = dustParams.vel_left;
    }
    for (int j = all_left; j < all_left + all_right; ++j)
    {
        image_drho[j] = dustParams.rho_right;
        image_dvel[j] = dustParams.vel_right;
    }

    for (int i = 0; i < dust_real_left_p_num; ++i)
    {
        drho[i] = dustParams.rho_left;
        dvel[i] = dustParams.vel_left;
    }
    for(int i = dust_real_left_p_num; i < dust_real_left_p_num + dust_real_right_p_num; ++i)
    {
        drho[i] = dustParams.rho_right;
        dvel[i] = dustParams.vel_right;
    }
}