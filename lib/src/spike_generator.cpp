#include "../include/spike_generator.hpp"

std::vector<float> Spike_Generator::generate_periodic_spikes(uint n, uint step, float value)
{
    if (step == 0.0) {
        std::vector<float> vec(n, value);
        return vec;
    }

    std::vector<float> vec;
    vec.reserve(n);

    for (uint i = 0; i < n; i++) {
        if ((i + 1) % step == 0)
            vec.push_back(value);
        else
            vec.push_back(0);
    }

    return vec;
}

std::vector<float> Spike_Generator::generate_random_spikes(uint n, float value)
{
    std::vector<float> vec;
    vec.reserve(n);

    for (uint i = 0; i < n; i++) {
        if ((rand() % 2) == 0)
            vec.push_back(value);
        else
            vec.push_back(0);
    }

    return vec;
}

std::vector<float> Spike_Generator::generate_poisson_spikes(uint n, float value)
{
    std::vector<float> vec;
    vec.reserve(n);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::poisson_distribution<int> pd(n/2);

	for (uint i = 0; i < n; i++) {
		int tmp = pd(gen);
		vec[tmp] = value;
	}

    return vec;
}
