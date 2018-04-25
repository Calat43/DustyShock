#include <time.h>
#include <assert.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include "common_use.hpp"
#include "gas_shock.hpp"
#include "explicit_shock.hpp"
#include "cells_shock.hpp"
#include <vector>

int main()
{
    ParticleParams gasParams;
    gasParams.particles_amount = 2250;  // should be devidable by (L2R + 1)
    gasParams.im_particles_amount = 450;  // should be devidable by (L2R + 1)
    gasParams.isGas = true;

    gasParams.rho_left = 1;
    gasParams.press_left = 1;
    gasParams.vel_left = 0;
    gasParams.energy_left = 2.5;

    gasParams.rho_right = 0.125;
    gasParams.press_right = 0.1;
    gasParams.vel_right = 0;
    gasParams.energy_right = 2;

    ParticleParams dustParams;
    dustParams.particles_amount = 2250;
    dustParams.im_particles_amount = 450;
    dustParams.isGas = false;

    dustParams.rho_left = 1;
    dustParams.vel_left = 0;
    dustParams.rho_right = 0.125;
    dustParams.vel_right = 0;

    ProblemParams problemParams;
    problemParams.T = 0.2;
    problemParams.h = 0.01;
    problemParams.tau = 0.0001;
    //problemParams.t_stop = 0.002;
    problemParams.K = 50000;

    problemParams.gamma = 1.4;
    problemParams.membrane = 0;
    problemParams.left = -0.5;
    problemParams.right = 0.5;
    problemParams.haveViscosity = true;
    problemParams.alfa = 1;
    problemParams.beta = 2;
    problemParams.nu_coef = 0.1;

    // check correctness of params
    double l2r = gasParams.rho_left / gasParams.rho_right;
    assert(is_close_to_int(fmod(gasParams.particles_amount, l2r + 1), __DOUBLE_ROUNDING_EPS));
    assert(is_close_to_int(fmod(gasParams.im_particles_amount, l2r + 1), __DOUBLE_ROUNDING_EPS));

    /*{  // этот код по данному предпочтительному кол-ву реальных частиц подбирает лучшие соотношения на все четыре части ({мнимые, реальные}x{левые, правые})
        double l2r = gasParams.rho_left / gasParams.rho_right;
        double r2i = gasParams.real2image_ratio;
        int best_im_right_num = (int) round(desired_particles_amount / (r2i + 1) / (l2r + 1));
        int best_re_right_num = (int) round(best_im_right_num * r2i);
        int best_im_left_num = (int) round(best_im_right_num * l2r);
        int best_re_left_num = (int) round(best_im_left_num * r2i);

        gasParams.particles_amount = best_re_left_num + best_re_right_num;
        printf("Using adjusted particles amount:\n");
        printf("Real amount    = %d\n", gasParams.particles_amount);
        printf("Overall amount = %d\n", gasParams.particles_amount + best_im_left_num + best_im_right_num);
        printf("  left image   = %d\n", best_im_left_num);
        printf("  left real    = %d\n", best_re_left_num);
        printf("  right real   = %d\n", best_re_right_num);
        printf("  right image  = %d\n", best_im_right_num);
        printf("  l2r          = %lf\n", l2r);
        printf("  r2i          = %lf\n", r2i);
   }*/


    clock_t startTime = clock();

    //gas_shock(gasParams, problemParams);
    cells_shock(0.005, gasParams, dustParams, problemParams);
    //explicit_shock(gasParams, dustParams, problemParams);

    clock_t finishTime = clock();

    double executionTime = (double)(finishTime - startTime) / CLOCKS_PER_SEC;
    printf("Finished in %lf seconds.\n", executionTime);

    return 0;
}