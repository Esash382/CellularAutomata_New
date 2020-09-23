#pragma once

#include "log.hpp"

#include <stdlib.h>
#include <random>
#include <algorithm>
#include <memory>
#include <ctime>
#include "util.hpp"
#include "spike_generator.hpp"

typedef enum _NETWORK_TYPE_ {
	EXCITATORY = 0,
	INHIBITORY
} NETWORK_TYPE;

class Network
{
public:
	Network(std::string name, uint ntwk_id, uint N, 
			NETWORK_TYPE type, float threshold,
			float k, _time_t tau_ap, _time_t tau_ref, uint z,
			_time_t tau_del, _time_t tau_dur, _time_t tau_osc,
			EXTERNAL_INPUT e_in, float value, float step);

	// Logger
	Log* logger;

public:
	std::string m_ntwk_name;
	uint m_ntwk_id;
	uint m_N;	// Total number of neurons

	uint m_Z; // number of connections amongst themselves

	float m_K_strength;	// Strength of the connection

	// Time parameters
	_time_t tau_ap;
	_time_t tau_ref;
	_time_t tau_del;
	_time_t tau_dur;
	_time_t tau_osc;
	_time_t last_osc_run_time;

	float threshold;
	NETWORK_TYPE m_type;

	EXTERNAL_INPUT external_input;
	float external_input_value;
	float ext_step;

	std::vector<float> external_input_vec;
	std::vector<_time_t> time_vec;
    std::vector<double> ac_stats;
    std::vector<double> in_stats;
    std::vector<double> ref_stats;
    
    // starting row index in the matrix
    // which means, 
    // ntwk->start_from_row_index + i = required specific neuron in this population
    uint start_from_row_index; 
};