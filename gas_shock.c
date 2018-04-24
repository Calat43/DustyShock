#include "gas_shock.h"

static double found_next_vel(double prev_vel, double * prev_image_vel, double mass, double prev_pressure, double prev_rho, int image_amount,
                             double * prev_image_pressure, double * prev_image_rho,
                             double prev_x, double * prev_image_x, ProblemParams params)
{
    double result = 0;

    for(int j = 0; j < image_amount; ++j)
    {
        if (fabs(prev_x - prev_image_x[j]) < 2.1 * params.h)
        {
            result += (prev_image_pressure[j] / pow(prev_image_rho[j], 2) + prev_pressure / pow(prev_rho, 2) +
                       found_viscosity(prev_pressure, prev_image_pressure[j], prev_vel, prev_image_vel[j],
                                       prev_rho, prev_image_rho[j], prev_x, prev_image_x[j], params))
                      * spline_gradient(prev_x, prev_image_x[j], params);
        }
    }

    return - params.tau * mass * result + prev_vel;
}

void gas_shock(ParticleParams gasParams, ProblemParams problemParams)
{
    int re_amount = gasParams.particles_amount;
    int all_amount = gasParams.particles_amount + gasParams.im_particles_amount;
    double prev_coord[re_amount];
    double image_prev_coord[all_amount];

    double next_coord[re_amount];
    double image_next_coord[all_amount];

    double l2r = gasParams.rho_left / gasParams.rho_right;

    // считаем количества реальных частиц
    int real_right_p_num = (int) round((double) gasParams.particles_amount / (l2r + 1));
    int real_left_p_num = gasParams.particles_amount - real_right_p_num;

    // аналогично для виртуальных
    int image_right_p_num = (int) round((double) gasParams.im_particles_amount / (l2r + 1));
    int image_left_p_num = gasParams.im_particles_amount - image_right_p_num;

    fill_initial_coord(prev_coord, image_prev_coord, real_left_p_num, real_right_p_num, image_left_p_num,
                       image_right_p_num, problemParams);

    double image_left_lenght = image_prev_coord[image_left_p_num + real_left_p_num] - image_prev_coord[0];
    double mass = found_mass(image_left_lenght, image_left_p_num + real_left_p_num, gasParams);


    /*
    // print coord so we can check it with gnuplot
    char filename[512];
    sprintf(filename, "%s/real_coord.dat", DATA_DIR);
    FILE * re_out = fopen(filename, "w");
    for (int i = 0; i < re_amount; i++) {
        fprintf(re_out, "%lf 2\n", coord[i]);
    }
    fclose(re_out);

    sprintf(filename, "%s/all_coord.dat", DATA_DIR);
    FILE * all_out = fopen(filename, "w");
    for (int i = 0; i < all_amount; i++) {
        fprintf(all_out, "%lf 1\n", image_coord[i]);
    }
    fclose(all_out);
    */

    double prev_rho[re_amount];
    double next_rho[re_amount];
    double prev_vel[re_amount];
    double next_vel[re_amount];
    double prev_energy[re_amount];
    double next_energy[re_amount];
    double prev_pressure[re_amount];
    double next_pressure[re_amount];

    double image_prev_rho[all_amount];
    double image_next_rho[all_amount];
    double image_prev_vel[all_amount];
    double image_next_vel[all_amount];
    double image_prev_energy[all_amount];
    double image_next_energy[all_amount];
    double image_prev_pressure[all_amount];
    double image_next_pressure[all_amount];

    fill_initial_gas_massives(prev_rho, prev_vel, prev_energy, prev_pressure, image_prev_rho, image_prev_vel,
                          image_prev_energy, image_prev_pressure, real_left_p_num, real_right_p_num,
                          image_left_p_num, image_right_p_num, gasParams);

    for (int i = 0; i < re_amount; ++i)
    {
        prev_rho[i] = found_next_rho(all_amount, mass, prev_coord[i], image_prev_coord, problemParams);
    }
    for (int j = 0; j < all_amount; ++j)
    {
        image_prev_rho[j] = found_next_rho(all_amount, mass, image_prev_coord[j], image_prev_coord, problemParams);
    }

    for(int i = 0; i < re_amount; ++i)
    {
        prev_pressure[i] = found_pressure(prev_rho[i], prev_energy[i], problemParams);
    }
    for(int j = 0; j < all_amount; ++j)
    {
        image_prev_pressure[j] = found_pressure(image_prev_rho[j], image_prev_energy[j], problemParams);
    }

    for(int j = 0; j < all_amount; ++j)
    {
        image_next_coord[j] = image_prev_coord[j];
        image_next_vel[j] = image_prev_vel[j];
        image_next_rho[j] = image_prev_rho[j];
        image_next_energy[j] = image_prev_energy[j];
        image_next_pressure[j] = image_prev_pressure[j];
    }

    char filename[512];
    /*
    sprintf(filename, "%s/reals.dat", DATA_DIR);
    FILE * re_out = fopen(filename, "w");
    for (int i = 0; i < all_amount; i++)
    {
        fprintf(re_out, "%lf %lf %lf %lf %lf\n", image_prev_coord[i], image_prev_rho[i], image_prev_vel[i], image_prev_energy[i],
                image_prev_pressure[i]);
    }
    fclose(re_out);
    */

    for (int frameId = 0; frameId < floor(problemParams.T / problemParams.tau); ++frameId)
    {
        printf("%d\n", frameId);

        for(int i = 0; i < re_amount; ++i)
        {
            next_coord[i] = found_next_coordinate(prev_coord[i], prev_vel[i], problemParams);
            next_vel[i] = found_next_vel(prev_vel[i], image_prev_vel, mass, prev_pressure[i], prev_rho[i], all_amount,
                                         image_prev_pressure, image_prev_rho, prev_coord[i],
                                         image_prev_coord, problemParams);
            next_energy[i] = found_next_energy(prev_energy[i], mass, prev_vel[i],  prev_pressure[i], image_prev_pressure,
                                               prev_rho[i], image_prev_rho, all_amount, image_prev_vel, prev_coord[i],
                                               image_prev_coord, problemParams);
            assert(!isnan(next_coord[i]));
            assert(!isnan(next_vel[i]));
            assert(!isnan(next_energy[i]));
        }

        for(int j = image_left_p_num; j < image_left_p_num + re_amount; ++j)
        {
            image_next_coord[j] = next_coord[j - image_left_p_num];
            image_next_vel[j] = next_vel[j - image_left_p_num];
            image_next_energy[j] = next_energy[j - image_left_p_num];

            assert(!isnan(image_next_coord[j]));
            assert(!isnan(image_next_vel[j]));
            assert(!isnan(image_next_energy[j]));
        }

        for(int i = 0; i < re_amount; ++i)
        {
            next_rho[i] = found_next_rho(all_amount, mass, next_coord[i], image_next_coord, problemParams);
            assert(!isnan(next_rho[i]));
        }
        for(int j = image_left_p_num; j < image_left_p_num + re_amount; ++j)
        {
            image_next_rho[j] = next_rho[j - image_left_p_num];
            assert(!isnan(image_next_rho[j]));
        }

        for(int i = 0; i < re_amount; ++i)
        {
            next_pressure[i] = found_pressure(next_rho[i], next_energy[i], problemParams);
            assert(!isnan(next_pressure[i]));
        }
        for(int j = image_left_p_num; j < image_left_p_num + re_amount; ++j)
        {
            image_next_pressure[j] = next_pressure[j - image_left_p_num];
            assert(!isnan(image_next_pressure[j]));
        }

        for(int i = 0; i < re_amount; ++i)
        {
            prev_coord[i] = next_coord[i];
            prev_vel[i] = next_vel[i];
            prev_rho[i] = next_rho[i];
            prev_energy[i] = next_energy[i];
            prev_pressure[i] = next_pressure[i];
        }

        for(int j = 0; j < all_amount; ++j)
        {
            image_prev_coord[j] = image_next_coord[j];
            image_prev_vel[j] = image_next_vel[j];
            image_prev_rho[j] = image_next_rho[j];
            image_prev_energy[j] = image_next_energy[j];
            image_prev_pressure[j] = image_next_pressure[j];
        }
    }
    sprintf(filename, "%s/shock_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg.dat", DATA_DIR, problemParams.T, problemParams.h,
            problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount, problemParams.nu_coef);
    FILE * res_out = fopen(filename, "w");
    for (int i = 0; i < re_amount; i++) {
        fprintf(res_out, "%lf %lf %lf %lf %lf\n", prev_coord[i], prev_rho[i], prev_vel[i], prev_energy[i],
                prev_pressure[i]);
    }
    fclose(res_out);

    sprintf(filename, "%s/im_shock_T%lg_h%lg_tau%lg_alfa%lg_beta%lg_N%d_nu%lg.dat", DATA_DIR, problemParams.T, problemParams.h,
            problemParams.tau, problemParams.alfa, problemParams.beta, gasParams.particles_amount, problemParams.nu_coef);
    FILE * imres_out = fopen(filename, "w");
    for (int i = 0; i < all_amount; i++) {
        fprintf(imres_out, "%lf %lf %lf %lf %lf\n", image_prev_coord[i], image_prev_rho[i], image_prev_vel[i],
                image_prev_energy[i], image_prev_pressure[i]);
    }
    fclose(imres_out);

}
