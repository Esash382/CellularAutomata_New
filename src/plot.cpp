#include "../include/plot.hpp"
#include <iostream>
#include <fstream>

// python plot.py

void Plot::plot_neurons_stats() {
    std::string command = "python scripts/plot.py";
    FILE* in = popen(command.c_str(), "r");
    pclose(in);
}
