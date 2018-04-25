#pragma once

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "common_use.hpp"

void cells_shock(double cell_length, ParticleParams gasParams, ParticleParams dustParams, ProblemParams problemParams);
