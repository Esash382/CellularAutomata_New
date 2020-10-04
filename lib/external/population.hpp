#pragma once

#include "../include/config.hpp"

#include <map>
#include <memory>

using std::unique_ptr;

typedef enum _RUN_TYPE_ {
	OSCILLATOR = 0,
	STIMULATED
} RUN_TYPE;

typedef enum _SYNAPSE_STATE_ {
    S_OFF = 0,
    S_ON = 1
} STYPE;

typedef enum _SYNAPSE_ACTIVATION_METHOD_ {
    RAND = 0,
    EINPUT
} SAM;

typedef std::multimap<_time_t, std::map<uint, std::vector<uint>>> mmap;

class Population {
public:
    Population();
    void get_time_vector(_time_t total_time, _time_t time_step);

    void add_network(std::string filepath);
    void add_network(std::string filepath, std::string foldername);
    void process_networks();

private:
    void calculate_number_of_synapses();
    void add_network_internal_func(Config *config);
    void calculate_network_size();
    void init_weight_matrix();
    void init_sstate_matrix();
    void init_sfiring_matrix();
    void init_nmatrices();
    void init_neuron_population(Config* config);
    void create_interneuronal_network_projections(Config* config);
    void create_recurrent_network_projections(Config* config);

    void synaptic_block(shared_ptr<Network> ntwk, uint n, uint i);
    void threshold_block(shared_ptr<Network> ntwk, uint n, uint i, MTYPE type);
    shared_ptr<Network> get_network(uint neuron_id);
    string get_network_name(uint neuron_id);

    void activate_single_neuron(shared_ptr<Network> ntwk, uint n, uint i);
    void activate_single_random_synapse(shared_ptr<Network> ntwk, uint n, uint i);
    void activate_synapses(uint n);
    void deactivate_synapses(uint n);
    void deactivate_neurons(uint n);
    double compute_weights(shared_ptr<Network> ntwk, uint n, uint i);

    void update_stats(uint n);
    void write_stats();

    _time_t get_noisy_delay(_time_t del);

    // getters and setters for tests
    // and some wrapper calls for private functions
public:
    uint get_network_size() const;
    std::vector<std::vector<double>> get_weight_matrix() const;

    void set_weight_matrix(std::vector<std::vector<double>> matrix);
    void set_state_matrix(uint row_start_index, uint i, STYPE state);
    void set_synapse_firing_time_matrix(uint row_start_index, uint i, double firing_time);
    void set_neuron_firing_time(shared_ptr<Network> ntwk, uint i, _time_t firing);
    void set_neuron_state(shared_ptr<Network> ntwk, uint i, MTYPE type);
    void call_activate_single_neuron(shared_ptr<Network> ntwk, uint n, uint i);
    void call_activate_synapses(uint n);
    void call_deactivate_synapses(uint n);
    void call_deactivate_neurons(uint n);
    void call_synaptic_block(shared_ptr<Network> ntwk, uint n, uint i);
    double call_compute_weights(shared_ptr<Network> ntwk, uint n, uint i);
    void call_threshold_block(shared_ptr<Network> ntwk, uint n, uint i, MTYPE type);
    void call_update_stats(uint n);

public:
	// Logger
	Log* logger;

    std::string filepath;
    std::string result_file_path;
    std::vector<_time_t> time_vec;
    _time_t total_time;
    _time_t time_step;
    _time_t last_osc_run_time;

    uint m_N;

    std::vector<shared_ptr<Network>> population_network;
    std::vector<std::vector<double>> m_w_matrix;

    std::vector<std::vector<STYPE>> m_s_matrix;
    std::vector<MTYPE> m_n_matrix;

    std::vector<std::vector<double>> m_sf_matrix;
    std::vector<double> m_nf_matrix;

    std::vector<double> weighted_sum;

     // Activate synapses: map <firing time + delay , map <neuron_id, vector<synapse_id >>>
    mmap s_activate;

    // Deactivate synapses: map <firing time + delay , map <neuron_id, vector<synapse_id >>>
    mmap s_deactivate;
    
    // Deactivate neuron: map <time , vector<neuron_id> >
    std::map<_time_t, std::vector<uint>> n_deactivate;
    
    // Refractory : map <time , vector<neuron_id> >
    std::map<_time_t, std::vector<uint>> n_refractory;
};

extern "C" Population* ext_create_population(const char* filepath) {
    auto popn = new Population();
    return popn;
}
