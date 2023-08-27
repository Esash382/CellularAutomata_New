#pragma once

#include "../include/config.hpp"

#include <set>
#include <memory>
#include <tuple>

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
    void set_firing_time_for_random_time_network();
    void init_bins();
    int get_bin_index(uint n);
    void init_p_rand_neurons();
    void init_intersecting_p_rand_neurons(uint p_rand_no_of_neurons,
                                          uint N,
                                          uint no_of_patterns,
                                          uint ntwk_id);
    void init_nonintersecting_p_rand_neurons(uint p_rand_no_of_neurons,
                                             uint N,
                                             uint no_of_patterns,
                                             uint start_from_row_index,
                                             uint ntwk_id);

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
    int get_random_pattern_index(std::vector<std::vector<uint>> recall_count);
    void get_random_pattern_index(std::vector<std::vector<uint>> recall_count, _time_t n);
    void get_recall_correlation(std::vector<std::vector<uint>> recall_count, uint pattern_size);
    void get_spurious_recall_count(_time_t n);

    _time_t get_noisy_delay(float del, float step);

    // getters and setters for tests
    // and some wrapper calls for private functions
public:
    uint get_network_size() const;
    std::vector<std::vector<double>> get_weight_matrix() const;

    void set_weight_matrix(std::vector<std::vector<double>> matrix);
    void set_p_rand_matrix(std::map<uint, std::vector<std::vector<uint>>> p_rand_neuron_ids);
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
    bool is_neuron_in_p_rand(uint i, uint ntwk_id, _time_t n);
    void learn(uint n);
	Log* logger; // Logger

    std::string filepath;
    std::string result_file_path;
    std::vector<_time_t> time_vec;
    _time_t total_time;
    _time_t time_step;
    uint bin_size;
    _time_t last_osc_run_time;

    const uint th_step;
    float del_step;

    uint m_N;

    std::vector<shared_ptr<Network>> population_network;
    std::vector<std::vector<double>> m_w_matrix;

    std::vector<std::vector<STYPE>> m_s_matrix;
    std::vector<MTYPE> m_n_matrix;

    std::vector<std::vector<double>> m_sf_matrix;
    std::vector<double> m_nf_matrix;

    std::vector<std::map<_time_t, std::set<uint>>> bin_values;

    std::vector<uint> m_spurious_recall_count;

     // Activate synapses: map <firing time + delay , map <neuron_id, vector<synapse_id >>>
    mmap s_activate;

    // Deactivate synapses: map <firing time + delay , map <neuron_id, vector<synapse_id >>>
    mmap s_deactivate;
    
    // Deactivate neuron: map <time , vector<neuron_id> >
    std::map<_time_t, std::vector<uint>> n_deactivate;
    
    // Refractory : map <time , vector<neuron_id> >
    std::map<_time_t, std::vector<uint>> n_refractory;

    // Activate neurons according to uniform distribution firing times
    std::map<uint, std::map<_time_t, std::vector<uint>>> n_rand_map;

    std::map<uint, std::vector<std::vector<uint>>> p_rand_neuron_ids;
    std::map<uint, std::tuple<uint, _time_t>> rand_pattern_times;

    // Store the current network id and the neuron in the network getting fired
    uint cur_ntwk_neuron;
    uint cur_ntwk_id;
    uint cur_ntwk_start_from_row_index;
    bool learning{false};
    std::pair<uint, uint> current_sensory_neuron; // <ntwk start index, ntwk neuron id>
    std::vector<std::tuple<uint, uint, uint, uint, uint>> ec_to_pyr_neurons; // std::vector<std::tuple<uint neuron_id, uint ntwk_id, uint start_index, learning_rate, unlearning_rate>>
    std::vector<std::tuple<uint, uint, uint>> ca3_neurons; // std::vector<std::tuple<uint neuron_id, uint ntwk_id, uint start_index>>
    float learning_inhibition_neuron_count{0};
    float learning_inhibition_total_neurons{0};
};

extern "C" Population* ext_create_population(const char* filepath) {
    auto popn = new Population();
    return popn;
}
