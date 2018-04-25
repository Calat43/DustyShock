#include "gas_shock.hpp"

static double found_next_vel(double prev_vel, std::vector<double> & prev_image_vel, double mass, double prev_pressure,
                             double prev_rho, uint image_amount, std::vector<double> & prev_image_pressure,
                             std::vector<double> & prev_image_rho, double prev_x, std::vector<double> & prev_image_x,
                             ProblemParams params)
{
    double result = 0;

    for(uint j = 0; j < image_amount; ++j)
    {
        if (fabs(prev_x - prev_image_x.at(j)) < 2.1 * params.h)
        {
            result += (prev_image_pressure.at(j) / pow(prev_image_rho.at(j), 2) + prev_pressure / pow(prev_rho, 2) +
                       found_viscosity(prev_pressure, prev_image_pressure.at(j), prev_vel, prev_image_vel.at(j),
                                       prev_rho, prev_image_rho.at(j), prev_x, prev_image_x.at(j), params))
                      * spline_gradient(prev_x, prev_image_x.at(j), params);
        }
    }

    return - params.tau * mass * result + prev_vel;
}

void gas_shock(ParticleParams gasParams, ProblemParams problemParams)
{
    uint re_amount = gasParams.particles_amount;
    uint all_amount = gasParams.particles_amount + gasParams.im_particles_amount;

    std::vector<double> prev_coord(re_amount);
    std::vector<double> image_prev_coord(all_amount);

    std::vector<double> next_coord(re_amount);
    std::vector<double> image_next_coord(all_amount);

    double l2r = gasParams.rho_left / gasParams.rho_right;

    // считаем количества реальных частиц
    uint real_right_p_num = (uint) round((double) gasParams.particles_amount / (l2r + 1));
    uint real_left_p_num = gasParams.particles_amount - real_right_p_num;

    // аналогично для виртуальных
    uint image_right_p_num = (uint) round((double) gasParams.im_particles_amount / (l2r + 1));
    uint image_left_p_num = gasParams.im_particles_amount - image_right_p_num;

    fill_initial_coord(prev_coord, image_prev_coord, real_left_p_num, real_right_p_num, image_left_p_num,
                       image_right_p_num, problemParams);

    double image_left_lenght = image_prev_coord.at(image_left_p_num + real_left_p_num) - image_prev_coord.at(0);
    double mass = found_mass(image_left_lenght, image_left_p_num + real_left_p_num, gasParams);


    /*
    // print coord so we can check it with gnuplot
    char filename[512];
    sprintf(filename, "%s/real_coord.dat", DATA_DIR);
    FILE * re_out = fopen(filename, "w");
    for (int i = 0; i < re_amount; i++) {
        fprintf(re_out, "%lf 2\n", coord.at(i));
    }
    fclose(re_out);

    sprintf(filename, "%s/all_coord.dat", DATA_DIR);
    FILE * all_out = fopen(filename, "w");
    for (int i = 0; i < all_amount; i++) {
        fprintf(all_out, "%lf 1\n", image_coord.at(i));
    }
    fclose(all_out);
    */

    std::vector<double> prev_rho(re_amount);
    std::vector<double> next_rho(re_amount);
    std::vector<double> prev_vel(re_amount);
    std::vector<double> next_vel(re_amount);
    std::vector<double> prev_energy(re_amount);
    std::vector<double> next_energy(re_amount);
    std::vector<double> prev_pressure(re_amount);
    std::vector<double> next_pressure(re_amount);

    std::vector<double> image_prev_rho(all_amount);
    std::vector<double> image_next_rho(all_amount);
    std::vector<double> image_prev_vel(all_amount);
    std::vector<double> image_next_vel(all_amount);
    std::vector<double> image_prev_energy(all_amount);
    std::vector<double> image_next_energy(all_amount);
    std::vector<double> image_prev_pressure(all_amount);
    std::vector<double> image_next_pressure(all_amount);

    fill_initial_gas_massives(prev_rho, prev_vel, prev_energy, prev_pressure, image_prev_rho, image_prev_vel,
                          image_prev_energy, image_prev_pressure, real_left_p_num, real_right_p_num,
                          image_left_p_num, image_right_p_num, gasParams);

    for (uint i = 0; i < re_amount; ++i)
    {
        prev_rho.at(i) = found_next_rho(all_amount, mass, prev_coord.at(i), image_prev_coord, problemParams);
    }
    for (uint j = 0; j < all_amount; ++j)
    {
        image_prev_rho.at(j) = found_next_rho(all_amount, mass, image_prev_coord.at(j), image_prev_coord, problemParams);
    }

    for(uint i = 0; i < re_amount; ++i)
    {
        prev_pressure.at(i) = found_pressure(prev_rho.at(i), prev_energy.at(i), problemParams);
    }
    for(uint j = 0; j < all_amount; ++j)
    {
        image_prev_pressure.at(j) = found_pressure(image_prev_rho.at(j), image_prev_energy.at(j), problemParams);
    }

    for(uint j = 0; j < all_amount; ++j)
    {
        image_next_coord.at(j) = image_prev_coord.at(j);
        image_next_vel.at(j) = image_prev_vel.at(j);
        image_next_rho.at(j) = image_prev_rho.at(j);
        image_next_energy.at(j) = image_prev_energy.at(j);
        image_next_pressure.at(j) = image_prev_pressure.at(j);
    }

    char filename[512];
    /*
    sprintf(filename, "%s/reals.dat", DATA_DIR);
    FILE * re_out = fopen(filename, "w");
    for (int i = 0; i < all_amount; i++)
    {
        fprintf(re_out, "%lf %lf %lf %lf %lf\n", image_prev_coord.at(i), image_prev_rho.at(i), image_prev_vel.at(i), image_prev_energy.at(i),
                image_prev_pressure.at(i));
    }
    fclose(re_out);
    */

    for (int frameId = 0; frameId < floor(problemParams.T / problemParams.tau); ++frameId)
    {
        printf("%d\n", frameId);

        for(uint i = 0; i < re_amount; ++i)
        {
            next_coord.at(i) = found_next_coordinate(prev_coord.at(i), prev_vel.at(i), problemParams);
            next_vel.at(i) = found_next_vel(prev_vel.at(i), image_prev_vel, mass, prev_pressure.at(i), prev_rho.at(i), all_amount,
                                         image_prev_pressure, image_prev_rho, prev_coord.at(i),
                                         image_prev_coord, problemParams);
            next_energy.at(i) = found_next_energy(prev_energy.at(i), mass, prev_vel.at(i),  prev_pressure.at(i), image_prev_pressure,
                                               prev_rho.at(i), image_prev_rho, all_amount, image_prev_vel, prev_coord.at(i),
                                               image_prev_coord, problemParams);
            assert(!isnan(next_coord.at(i)));
            assert(!isnan(next_vel.at(i)));
            assert(!isnan(next_energy.at(i)));
        }

        for(uint j = image_left_p_num; j < image_left_p_num + re_amount; ++j)
        {
            image_next_coord.at(j) = next_coord.at(j - image_left_p_num);
            image_next_vel.at(j) = next_vel.at(j - image_left_p_num);
            image_next_energy.at(j) = next_energy.at(j - image_left_p_num);

            assert(!isnan(image_next_coord.at(j)));
            assert(!isnan(image_next_vel.at(j)));
            assert(!isnan(image_next_energy.at(j)));
        }

        for(uint i = 0; i < re_amount; ++i)
        {
            next_rho.at(i) = found_next_rho(all_amount, mass, next_coord.at(i), image_next_coord, problemParams);
            assert(!isnan(next_rho.at(i)));
        }
        for(uint j = image_left_p_num; j < image_left_p_num + re_amount; ++j)
        {
            image_next_rho.at(j) = next_rho.at(j - image_left_p_num);
            assert(!isnan(image_next_rho.at(j)));
        }

        for(uint i = 0; i < re_amount; ++i)
        {
            next_pressure.at(i) = found_pressure(next_rho.at(i), next_energy.at(i), problemParams);
            assert(!isnan(next_pressure.at(i)));
        }
        for(uint j = image_left_p_num; j < image_left_p_num + re_amount; ++j)
        {
            image_next_pressure.at(j) = next_pressure.at(j - image_left_p_num);
            assert(!isnan(image_next_pressure.at(j)));
        }

        for(uint i = 0; i < re_amount; ++i)
        {
            prev_coord.at(i) = next_coord.at(i);
            prev_vel.at(i) = next_vel.at(i);
            prev_rho.at(i) = next_rho.at(i);
            prev_energy.at(i) = next_energy.at(i);
            prev_pressure.at(i) = next_pressure.at(i);
        }

        for(uint j = 0; j < all_amount; ++j)
        {
            image_prev_coord.at(j) = image_next_coord.at(j);
            image_prev_vel.at(j) = image_next_vel.at(j);
            image_prev_rho.at(j) = image_next_rho.at(j);
            image_prev_energy.at(j) = image_next_energy.at(j);
            image_prev_pressure.at(j) = image_next_pressure.at(j);
        }
    }
    sprintf(filename, "%s/shock_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg.dat", DATA_DIR, problemParams.T, problemParams.h,
            problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount, problemParams.nu_coef);
    FILE * res_out = fopen(filename, "w");
    for (uint i = 0; i < re_amount; i++) {
        fprintf(res_out, "%lf %lf %lf %lf %lf\n", prev_coord.at(i), prev_rho.at(i), prev_vel.at(i), prev_energy.at(i),
                prev_pressure.at(i));
    }
    fclose(res_out);

    sprintf(filename, "%s/im_shock_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg.dat", DATA_DIR, problemParams.T, problemParams.h,
            problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount, problemParams.nu_coef);
    FILE * imres_out = fopen(filename, "w");
    for (uint i = 0; i < all_amount; i++) {
        fprintf(imres_out, "%lf %lf %lf %lf %lf\n", image_prev_coord.at(i), image_prev_rho.at(i), image_prev_vel.at(i),
                image_prev_energy.at(i), image_prev_pressure.at(i));
    }
    fclose(imres_out);

}
