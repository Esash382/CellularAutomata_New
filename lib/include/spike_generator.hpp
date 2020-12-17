#pragma once

#include <vector>
#include <random>

#include "util.hpp"

using std::vector;

class Spike_Generator
{
    public:
        static std::vector<float> generate_periodic_spikes(uint n, uint step, float value);
        static std::vector<float> generate_random_spikes(uint n, float value);
        static std::vector<float> generate_poisson_spikes(uint n, float value);
        static std::vector<float> generate_continuous_spikes(uint n, uint step, float value);
};
