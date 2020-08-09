// CellularAutomata Theta model project:
// main.cpp : Defines the entry point for the console application.

#include <cstdlib>
#include <time.h>
#include <vector>
#include <map>
#include <dlfcn.h>

#include "../lib/external/population.hpp"
#include "../include/plot.hpp"

typedef Population* (*populationCreatorFunction)(const char* filepath);

void init_and_process_networks(const char* filepath) {
	void* handle = dlopen("libca.so", RTLD_NOW);
	char* error;

	if (!handle) {
		fputs (dlerror(), stderr);
		exit(1);
    }

	populationCreatorFunction create = (populationCreatorFunction) dlsym(handle, "ext_create_population");
	if ((error = dlerror()) != NULL)  {
		fputs(error, stderr);
		exit(1);
	}

	clock_t start = clock();

	Population* pp_network = (*create)(filepath);

	if (!pp_network) {
		dlclose(handle);
		return;
	}

    pp_network->add_network(filepath);
	pp_network->process_networks();

	printf("Time taken: %.2f mins\n", (double)(clock() - start)/(CLOCKS_PER_SEC * 60));

	free(pp_network);
	dlclose(handle);
}

int main()
{
	const char* filepath = "/home/esash/Documents/Projects/CellularAutomata/CellularAutomata_Fast/";
	init_and_process_networks(filepath);

	Plot::plot_neurons_stats();

	return 0;
}
