#include "monaghan_shock.h"

#include "explicit_shock.hpp"

static double found_next_gas_vel(double prev_gvel, double prev_grho, double prev_gpress, double prev_gcoord, double gmass,
                                 double dmass, std::vector<double> & image_prev_gpress, std::vector<double> & image_prev_grho,
                                 std::vector<double> & image_prev_gcoord, std::vector<double> & image_prev_gvel,
                                 std::vector<double> & image_prev_drho, std::vector<double> & image_prev_dvel,
                                 std::vector<double> & image_prev_dcoord, int image_gamount, int image_damount,
                                 ProblemParams params)


{
    double a_p = 0;
    double drag = 0;
    double result = 0;
    double r_ja = 0;

    for(uint j = 0; j < image_gamount; ++j)
    {
        if (fabs(prev_gcoord - image_prev_gcoord.at(j)) <= 2. * params.h)
        {
            a_p += (prev_gpress / pow(prev_grho, 2) + image_prev_gpress.at(j) / pow(image_prev_grho.at(j), 2) +
                    found_viscosity(prev_gpress, image_prev_gpress.at(j), prev_gvel, image_prev_gvel.at(j),
                                    prev_grho, image_prev_grho.at(j), prev_gcoord, image_prev_gcoord.at(j), params))
                   * spline_gradient(prev_gcoord, image_prev_gcoord.at(j), params);
        }
    }
    a_p *= gmass;

    for(uint j = 0; j < image_damount; ++j)
    {
        if (fabs(prev_gcoord - image_prev_dcoord.at(j)) <= 2. * params.h)
        {
            r_ja = image_prev_dcoord.at(j) - prev_gcoord;
            drag += (prev_gvel - image_prev_dvel.at(j)) * r_ja / (pow(r_ja, 2) + 0.001 * params.h * params.h) /
                    prev_grho /
                    image_prev_drho.at(j) * r_ja * spline_kernel(prev_gcoord, image_prev_dcoord.at(j), params);
        }
    }
    drag = drag * dmass * params.K;

    result = - params.tau * a_p - params.tau * drag + prev_gvel;

    return result;
}

static double found_next_dust_vel(double prev_dvel, double prev_drho, double prev_dcoord, double gmass,
                                  std::vector<double> & image_prev_gvel, std::vector<double> & image_prev_grho,
                                  std::vector<double> & image_prev_gcoord, int image_gamount, ProblemParams params)
{
    double result = 0;
    double r_ja = 0;
    for(uint j = 0; j < image_gamount; ++j)
    {
        if (fabs(prev_dcoord - image_prev_gcoord.at(j)) <= 2. * params.h)
        {
            r_ja = prev_dcoord - image_prev_gcoord.at(j);
            result += (image_prev_gvel.at(j) - prev_dvel) * r_ja / (pow(r_ja, 2) + 0.001 * params.h * params.h) /
                    prev_drho / image_prev_grho.at(j)  * r_ja * spline_kernel(prev_dcoord, image_prev_gcoord.at(j), params);
        }
    }

    return params.tau * gmass * params.K * result + prev_dvel;
}

void monaghan_shock(ParticleParams gasParams, ParticleParams dustParams, ProblemParams problemParams)
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

    //а теперь смотрим, когда количество частиц газа = количеству частиц пыли

    char filename[512];

    sprintf(filename, "%s/im_monaghanShock_gas_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * gas0_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(gas0_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord.at(i), image_prev_grho.at(i), image_prev_gvel.at(i),
                image_prev_energy.at(i), image_prev_pressure.at(i));
    }
    fclose(gas0_out);

    sprintf(filename, "%s/im_monaghanShock_dust_T0_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * dust0_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(dust0_out, "%lf %lf %lf\n", image_prev_dcoord.at(i), image_prev_drho.at(i), image_prev_dvel.at(i));
    }
    fclose(dust0_out);

    for(int frameId = 0; frameId < floor(problemParams.T / problemParams.tau); ++frameId)
    {
        printf("%d\n", frameId);

        for(uint i = 0; i < re_gamount; ++i)
        {
            next_gcoord.at(i) = found_next_coordinate(prev_gcoord.at(i), prev_gvel.at(i), problemParams);
            next_gvel.at(i) = found_next_gas_vel(prev_gvel.at(i), prev_grho.at(i), prev_pressure.at(i), prev_gcoord.at(i),
                                                 gmass, dmass, image_prev_pressure, image_prev_grho, image_prev_gcoord,
                                                 image_prev_gvel, image_prev_drho, image_prev_dvel, image_prev_dcoord,
                                                 all_gamount, all_damount, problemParams);
            next_energy.at(i) = found_next_energy(prev_energy.at(i), gmass, prev_gvel.at(i), prev_pressure.at(i),
                                                  image_prev_pressure, prev_grho.at(i), image_prev_grho, all_gamount,
                                                  image_prev_gvel, prev_gcoord.at(i), image_prev_gcoord, problemParams);
            assert(!isnan(next_gcoord.at(i)));
            assert(!isnan(next_gvel.at(i)));
            assert(!isnan(next_energy.at(i)));
        }
        for(uint i = 0; i < re_damount; ++i)
        {
            next_dcoord.at(i) = found_next_coordinate(prev_dcoord.at(i), prev_dvel.at(i), problemParams);
            next_dvel.at(i) = found_next_dust_vel(prev_dvel.at(i), prev_drho.at(i), prev_dcoord.at(i), gmass,
                                                  image_prev_gvel, image_prev_grho, image_prev_gcoord, all_gamount,
                                                  problemParams);
            assert(!isnan(next_dcoord.at(i)));
            assert(!isnan(next_dvel.at(i)));
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
        //gasDrag = found_gvel_drag(prev_grho, prev_gvel, prev_dvel_in_gas, problemParams.K, re_gamount);
        //dustDrag = found_dvel_drag(prev_drho, prev_dvel, prev_gvel_in_dust, problemParams.K, re_damount);

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

    }

    sprintf(filename, "%s/im_monaghanShock_gas_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * gas_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(gas_out, "%lf %lf %lf %lf %lf\n", image_prev_gcoord.at(i), image_prev_grho.at(i), image_prev_gvel.at(i),
                image_prev_energy.at(i), image_prev_pressure.at(i));
    }
    fclose(gas_out);

    sprintf(filename, "%s/im_monaghanShock_dust_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg_K%lg.dat", DATA_DIR, problemParams.T,
            problemParams.h, problemParams.tau, problemParams.alfa, problemParams.beta, dustParams.particles_amount,
            problemParams.nu_coef, problemParams.K);
    FILE * dust_out = fopen(filename, "w");
    for (uint i = 0; i < all_gamount; i++) {
        fprintf(dust_out, "%lf %lf %lf\n", image_prev_dcoord.at(i), image_prev_drho.at(i), image_prev_dvel.at(i));
    }
    fclose(dust_out);

}

