#include "explicit_shock.h"

static double found_next_gas_vel(double prev_gvel, double prev_grho, double prev_gpress, double prev_gcoord,
                                 double gmass, double * image_prev_gvel, double * image_prev_grho,
                                 double * image_prev_gpress, double * image_prev_gcoord,  int image_amount,
                                 double prev_dvel_in_gas, double epsilon, ProblemParams params)
{
    double result = 0;
    for(int j = 0; j < image_amount; ++j)
    {
        if (fabs(prev_gcoord - image_prev_gcoord[j]) < 2.1 * params.h)
        {
            result += (image_prev_gpress[j] / pow(image_prev_grho[j], 2) + prev_gpress/ pow(prev_grho, 2) +
                       found_viscosity(prev_gpress, image_prev_gpress[j], prev_gvel, image_prev_gvel[j],
                                       prev_grho, image_prev_grho[j], prev_gcoord, image_prev_gcoord[j], params))
                      * spline_gradient(prev_gcoord, image_prev_gcoord[j], params);
        }
    }

    return - params.tau * gmass * result - params.tau * params.K * (prev_gvel - prev_dvel_in_gas) / prev_grho +
            prev_gvel;
}

static double found_next_dust_vel(double prev_dvel, double prev_drho, double prev_gvel_in_dust, ProblemParams params)
{
    return params.tau * params.K * (prev_gvel_in_dust - prev_dvel) / prev_drho + prev_dvel;
}

double near_velocity(double r, double * coordinate, double * velocity, int amount)
{
    int nearly = 0;

    for (int i = 0; i < amount; ++i)
    {
        double new_fabs = fabs(r - coordinate[i]);
        double old_fabs = fabs(r - coordinate[nearly]);
        if (new_fabs < old_fabs)
        {
            nearly = i;
        }
    }

    return velocity[nearly];
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

static double found_gvel_drag(double * grho, double * gvel, double * dvel_in_gas, double K, int re_gamount)
{
    double result = 0;

    for(int i = 0; i < re_gamount; ++i)
    {
        result += (gvel[i] - dvel_in_gas[i]) / grho[i];
    }
    return result * K;
}

static double found_dvel_drag(double * drho, double * dvel, double * gvel_in_dust, double K, int re_damount)
{
    double result = 0;

    for(int i = 0; i < re_damount; ++i)
    {
        result += (dvel[i] - gvel_in_dust[i]) / drho[i];
    }
    return result * K;
}

void explicit_shock(ParticleParams gasParams, ParticleParams dustParams, ProblemParams problemParams)
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

    //а теперь смотрим, когда количество частиц газа = количеству частиц пыли
    double epsilon[re_gamount];
    double prev_gvel_in_dust[re_damount];
    double prev_dvel_in_gas[re_gamount];

    double gasDrag = 0;
    double dustDrag = 0;
    double momentum = 0;

    for(int i = 0; i < re_gamount; ++i)
    {
        epsilon[i] = found_epsilon(prev_dcoord[i], prev_dcoord, dmass, prev_grho[i], re_damount, problemParams);
    }

    char filename[512];

    sprintf(filename, "%s/im_explShock_gas_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * gas0_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(gas0_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord[i], image_prev_grho[i], image_prev_gvel[i],
                image_prev_energy[i], image_prev_pressure[i]);
    }
    fclose(gas0_out);

    sprintf(filename, "%s/im_explShock_dust_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * dust0_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(dust0_out, "%lf %lf %lf\n", image_prev_dcoord[i], image_prev_drho[i], image_prev_dvel[i]);
    }
    fclose(dust0_out);

    //sprintf(filename, "%s/near_momentum_tau%lg.dat", DATA_DIR, problemParams.tau);
    //FILE * moment = fopen(filename, "w");

    for(int frameId = 0; frameId < floor(problemParams.T / problemParams.tau); ++frameId)
    {
        printf("%d\n", frameId);

        for(int i = 0; i < re_damount; ++i)
        {
            prev_gvel_in_dust[i] = //near_velocity(prev_dcoord[i], image_prev_gcoord, image_prev_gvel, all_gamount);
                    interpolation_value(prev_dcoord[i], gmass, image_prev_gvel, image_prev_grho,
                    image_prev_gcoord, all_gamount, problemParams);
        }
        for(int i = 0; i < re_gamount; ++i)
        {
            prev_dvel_in_gas[i] = //near_velocity(prev_gcoord[i], image_prev_dcoord, image_prev_dvel, all_damount);
                    interpolation_value(prev_gcoord[i], dmass, image_prev_dvel, image_prev_drho,
                    image_prev_dcoord, all_damount, problemParams);
        }

        for(int i = 0; i < re_gamount; ++i)
        {
            next_gcoord[i] = found_next_coordinate(prev_gcoord[i], prev_gvel[i], problemParams);
            next_gvel[i] = found_next_gas_vel(prev_gvel[i], prev_grho[i], prev_pressure[i], prev_gcoord[i], gmass,
                                              image_prev_gvel, image_prev_grho, image_prev_pressure, image_prev_gcoord,
                                              all_gamount, prev_dvel_in_gas[i], epsilon[i], problemParams);
            next_energy[i] = found_next_energy(prev_energy[i], gmass, prev_gvel[i], prev_pressure[i],
                                               image_prev_pressure, prev_grho[i], image_prev_grho, all_gamount,
                                               image_prev_gvel, prev_gcoord[i], image_prev_gcoord, problemParams);
            assert(!isnan(next_gcoord[i]));
            assert(!isnan(next_gvel[i]));
            assert(!isnan(next_energy[i]));
        }
        for(int i = 0; i < re_damount; ++i)
        {
            next_dcoord[i] = found_next_coordinate(prev_dcoord[i], prev_dvel[i], problemParams);
            next_dvel[i] = found_next_dust_vel(prev_dvel[i], prev_drho[i], prev_gvel_in_dust[i], problemParams);
            assert(!isnan(next_dcoord[i]));
            assert(!isnan(next_dvel[i]));
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
        //gasDrag = found_gvel_drag(prev_grho, prev_gvel, prev_dvel_in_gas, problemParams.K, re_gamount);
        //dustDrag = found_dvel_drag(prev_drho, prev_dvel, prev_gvel_in_dust, problemParams.K, re_damount);

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

    }

    sprintf(filename, "%s/im_explShock_gas_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * gas_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(gas_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord[i], image_prev_grho[i], image_prev_gvel[i],
                image_prev_energy[i], image_prev_pressure[i]);
    }
    fclose(gas_out);

    sprintf(filename, "%s/im_explShock_dust_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * dust_out = fopen(filename, "w");
    for (int i = 0; i < all_gamount; i++) {
        fprintf(dust_out, "%lf %lf %lf\n", image_prev_dcoord[i], image_prev_drho[i], image_prev_dvel[i]);
    }
    fclose(dust_out);

    //fclose(moment);
}
