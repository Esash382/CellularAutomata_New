#pragma once

#include <vector>

class Plot
{
public:
	static void plot_neurons_stats();
	static void plot_weight_matrix(std::vector<std::vector<float>> w_matrix);
};
