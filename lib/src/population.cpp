#include "../external/population.hpp"
#include "../include/corr_coeff.hpp"

#include <bits/c++config.h>
#include <functional>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <chrono>
#include <ostream>

#define ASSERT_WITH_MESSAGE(condition, message) do { \
if (!(condition)) { std::cout << ((message)) << std::endl; } \
assert ((condition)); } while(false)

const _time_t min_time_required_for_simulation = 100;

using namespace std;

Population::Population()
    : th_step(10)
    , del_step(5)
{
	logger = Log::getInstance();

	this->time_vec.reserve(0);
	this->filepath = "";
	this->result_file_path = "";
    this->bin_size = 0;

    logger->log("Population constructor");
}

void Population::get_time_vector(_time_t total_time, _time_t time_step)
{
    this->time_vec = _time_::getInstance(total_time, time_step);
    logger->log("Time vector: size : " + std::to_string(total_time) + " : step : " + std::to_string(time_step) + " : total time vector size : " + std::to_string(this->time_vec.size()));
}

uint Population::get_network_size() const
{
    return m_N;
}

std::vector<std::vector<double>> Population::get_weight_matrix() const
{
    return this->m_w_matrix;
}

void Population::set_weight_matrix(std::vector<std::vector<double>> matrix)
{
    this->m_w_matrix = matrix;
}

void Population::set_p_rand_matrix(std::map<uint, std::vector<std::vector<uint>>> p_rand_neuron_ids)
{
    this->p_rand_neuron_ids = p_rand_neuron_ids;
}

void Population::set_state_matrix(uint row_start_index, uint i, STYPE state)
{
    this->m_s_matrix[row_start_index][i] = state;
}

void Population::set_synapse_firing_time_matrix(uint row_start_index, uint i, double firing_time)
{
    this->m_sf_matrix[row_start_index][i] = firing_time;
}

void Population::calculate_network_size()
{
    uint sum = 0;
    for (auto cell : this->population_network) {
        sum += cell->m_N;
    }
    this->m_N = sum;

    logger->log("Total network size : " + std::to_string(this->m_N));
}

void Population::init_weight_matrix()
{
    for (uint i = 0; i < this->m_N; i++)
    {
        std::vector<double> vec(this->m_N);
        this->m_w_matrix.push_back(vec);
    }

    logger->log("Weight matrix initialized : size : " + std::to_string(this->m_w_matrix.size()));
}

void Population::init_sstate_matrix()
{
    for (uint i = 0; i < this->m_N; i++)
    {
        std::vector<STYPE> vec(this->m_N);
        this->m_s_matrix.push_back(vec);
    }

    logger->log("synapse state matrix initialized : size : " + std::to_string(this->m_s_matrix.size()));
}

void Population::init_sfiring_matrix()
{
    for (uint i = 0; i < this->m_N; i++)
    {
        std::vector<double> vec(this->m_N, 0.0);
        this->m_sf_matrix.push_back(vec);
    }

    logger->log("synapse firing matrix initialized : size : " + std::to_string(this->m_sf_matrix.size()));
}

void Population::init_nmatrices()
{
    this->m_n_matrix.reserve(this->m_N);
    this->m_nf_matrix.reserve(this->m_N);
    for (uint i = 0; i < this->m_N; i++) {
        this->m_n_matrix[i] = OFF;
        this->m_nf_matrix[i] = 0.0f;
    }

    logger->log("neuron firing matrix initialized : size : " + std::to_string(this->m_nf_matrix.size()));
    logger->log("neuron state matrix initialized : size : " + std::to_string(this->m_n_matrix.size()));
}

void Population::init_neuron_population(Config* config)
{
    auto unique_files = config->get_unique_files();
    uint prev_ntwk_end_index = 0;

    this->population_network.reserve(unique_files.size());

    for (uint i = 0; i < unique_files.size(); i++) {
        logger->log("files = " + unique_files[i]);
        shared_ptr<Network> ntwk = config->create_network(unique_files[i]);
        if (i == 0) {
            ntwk->start_from_row_index = 0;
        } else {
            ntwk->start_from_row_index = prev_ntwk_end_index;
        }
        prev_ntwk_end_index = ntwk->start_from_row_index + ntwk->m_N;
        this->population_network.push_back(std::move(ntwk));
    }

    logger->log("Neuron population initialized");
}

void Population::init_bins()
{
    //uint number_of_bins = this->total_time / bin_size;
    // For each network
    for (uint i = 0; i < this->population_network.size(); i++) {
        std::map<_time_t, std::set<uint>> tmap;

        // For each bin_size, put time in the map and the set
        for (uint j = bin_size; j <= this->total_time; j += bin_size) {
            std::set<uint> tset;
            tmap[j] = tset;
        }

        this->bin_values.push_back(tmap);
    }
}

void Population::create_interneuronal_network_projections(Config* config)
{
    logger->log("create_interneuronal_network_projections");

    auto network_files = config->get_files();
    for (auto cell : network_files) {
        logger->log("Network files : " + cell[0] + " - " + cell[1]);

        std::map<string, float> con = config->read_file(cell[0]+"_"+cell[1]);
        float ksd = con["ksd"];    // source to destination - number of projections
        float kds = con["kds"];    // destination to source - number of projections
        float zsd = con["zsd"];    // source to destination - connection weight
        float zds = con["zds"];    // source to destination - connection weight

        uint src_size = 0;
        uint src_start_from_row_index = 0;
        bool src_ext_pseudo_neuron = false;

        uint dst_size = 0;
        uint dst_start_from_row_index = 0;
        bool dst_ext_pseudo_neuron = false;

        for (auto n : this->population_network) {
            if (n->m_ntwk_name == cell[0]) {
                src_size = n->m_N;
                src_start_from_row_index = n->start_from_row_index;
                if (n->external_input == RANDOM_TIME)
                    src_ext_pseudo_neuron = true;
                zsd = (zsd * src_size) / 100.0;
            } else if (n->m_ntwk_name == cell[1]) {
                dst_size = n->m_N;
                dst_start_from_row_index = n->start_from_row_index;
                if (n->external_input == RANDOM_TIME)
                    dst_ext_pseudo_neuron = true;
                zds = (zds * dst_size) / 100.0;
            }

            if ((src_size != 0) && (dst_size != 0))
                break;
        }

//        std::cout << cell[0] << " = " << src_start_from_row_index << " : " << cell[1] << " = " << dst_start_from_row_index << std::endl;

        if (zsd != 0 || src_ext_pseudo_neuron) {
            for (uint i = src_start_from_row_index; i < (src_start_from_row_index + src_size); i++) {
                uint index = 0;
                std::vector<unsigned int> indices(dst_size);
                std::iota(indices.begin(), indices.end(), dst_start_from_row_index);
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));

                uint t_zsd = zsd;
                if (zsd == 0.0 && src_ext_pseudo_neuron)
                    t_zsd = get_random_number(1, zsd);

                for (uint k = 0; k < t_zsd; k++) {
//                    std::cout << "src = " << cell[0] << " EXT_id = " << i << " dst = " << cell[1] << " NEURON_id = " << k << "; z = " << t_zsd << std::endl;
                    uint j = indices[k];
                    if (this->m_w_matrix[i][j] == 0) {
                        std::string str = "There aren't enough number of neurons in the " + cell[0] + " to project to from " + cell[1];
                        ASSERT_WITH_MESSAGE((this->m_w_matrix[i][j] == 0), str);
                    }
                    this->m_w_matrix[i][j] = ksd;
                    this->m_s_matrix[i][j] = S_OFF;
                    this->m_sf_matrix[i][j] = 0.0f;
                    //std::cout << cell[0] << " -> " << cell[1] << " = " << zsd << " : " << index << " = " << k << " : this->m_w_matrix[" << i << "][" << j << "] = " << this->m_w_matrix[i][j] << std::endl;
                    index++;
                }

                /*
                string str = "";
                for (uint count = 0; count < indices.size(); count++) {
                    if (count == 0)
                        str += "[ " + std::to_string(indices[count]) + ", ";
                    else if (count == indices.size() - 1)
                        str += std::to_string(indices[count]) + " ]";
                    else
                        str += std::to_string(indices[count]) + " , ";
                }

                std::cout << cell[0] << " -> " << cell[1] << " = " << zsd << " : " << index << " = " << str << std::endl;
                */
            }
        }
        
        if (zds != 0 || dst_ext_pseudo_neuron) {
            for (uint i = dst_start_from_row_index; i < (dst_start_from_row_index + dst_size); i++) {
                uint index = 0;
                std::vector<unsigned int> indices(src_size);
                std::iota(indices.begin(), indices.end(), src_start_from_row_index);
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));

                uint t_zds = zds;
                if (zds == 0.0 && dst_ext_pseudo_neuron)
                    t_zds = get_random_number(1, zds);

                for (uint k = 0; k < t_zds; k++) {
//                    std::cout << "src = " << cell[1] << " EXT_id = " << i << " dst = " << cell[0] << " NEURON_id = " << k << "; z = " << t_zds << std::endl;
                    uint j = indices[k];
                    if (this->m_w_matrix[i][j] == 0) {
                        std::string str = "There aren't enough number of neurons in the " + cell[1] + " to project to from " + cell[0];
                        ASSERT_WITH_MESSAGE((this->m_w_matrix[i][j] == 0), str);
                    }
                    this->m_w_matrix[i][j] = kds;
                    this->m_s_matrix[i][j] = S_OFF;
                    this->m_sf_matrix[i][j] = 0.0f;
                    //std::cout << cell[1] << " -> " << cell[0] << " = " << zds << " : " << index << " = " << k << " : this->m_w_matrix[" << i << "][" << j << "] = " << this->m_w_matrix[i][j] << std::endl;
                    index++;
                }

                /*
                string str = "";
                for (uint count = 0; count < indices.size(); count++) {
                    if (count == 0)
                        str += "[ " + std::to_string(indices[count]) + ", ";
                    else if (count == indices.size() - 1)
                        str += std::to_string(indices[count]) + " ]";
                    else
                        str += std::to_string(indices[count]) + " , ";
                }

                std::cout << cell[1] << " -> " << cell[0] << " = " << zds << " : " << index << " = " << str << std::endl;
                */
            }
        }
    }
}

void Population::create_recurrent_network_projections(Config *config)
{
    logger->log("create_recurrent_network_projections");
    auto unique_files = config->get_unique_files();
    for (auto cell : this->population_network) {
        if (cell->m_Z != 0) {
            logger->log("recurrent connections are present for : " + cell->m_ntwk_name);
            uint size = cell->m_N;
            uint start_index = cell->start_from_row_index;

            for (uint i = start_index; i < (start_index + size); i++) {
                std::vector<unsigned int> indices(size);
                std::iota(indices.begin(), indices.end(), start_index);
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));

                for (uint k = 0; k < cell->m_Z; k++) {
                   uint j = indices[k];
                   if (i == j)
                       continue;
                   if (this->m_w_matrix[i][j] == 0) {
                       this->m_w_matrix[i][j] = cell->m_K_strength;
                       this->m_s_matrix[i][j] = S_OFF;
                       this->m_sf_matrix[i][j] = 0.0f;
                   } 
                }     
            }
            
        }
    }
}

// THUMB rule: number of external input neurons < size of time_vec array
/*
void Population::set_firing_time_for_random_time_network()
{
    for (auto ntwk : this->population_network) {
        if (ntwk->external_input == RANDOM_TIME) {
            std::map<time_t, uint> n_rand_activate;

            std::vector<unsigned int> t_indices(this->time_vec.size());
            std::iota(t_indices.begin(), t_indices.end(), this->time_vec[0]);
            unsigned tseed = std::chrono::system_clock::now().time_since_epoch().count();
            std::shuffle(t_indices.begin(), t_indices.end(), std::default_random_engine(tseed));

            uint total_possible_firing_repeat_times = this->time_vec.size() / ntwk->m_N;

            if (ntwk->number_of_firing_times > total_possible_firing_repeat_times) {
                std::cout << "Error: Number of firing times is greater than the length of time vector!" << std::endl;
                exit(1);
            }

            for (uint j= 0;
                 (j < ntwk->number_of_firing_times &&
                  j <= total_possible_firing_repeat_times);
                  j++) {
                std::vector<unsigned int> n_indices(ntwk->m_N);
                std::iota(n_indices.begin(), n_indices.end(), 0);
                unsigned nseed = std::chrono::system_clock::now().time_since_epoch().count();
                std::shuffle(n_indices.begin(), n_indices.end(), std::default_random_engine(nseed));

                if (ntwk->m_N < this->time_vec.size()) {
                    for (uint i = 0; i < ntwk->m_N; i++) {
                        if (j > 0)
                            n_rand_activate[this->time_vec[t_indices[i+j*ntwk->m_N]]] = n_indices[i];
                        else
                            n_rand_activate[this->time_vec[t_indices[i]]] = n_indices[i];
                    }
                } else {
                    for (uint i = 0; i < this->time_vec.size(); i++) {
                        if (j > 0)
                            n_rand_activate[this->time_vec[t_indices[i+j*ntwk->m_N]]] = n_indices[i];
                        else
                            n_rand_activate[this->time_vec[t_indices[i]]] = n_indices[i];
                    }
                }
            }

            n_rand_map[ntwk->m_ntwk_id] = n_rand_activate;
        }
    }

            // Each neuron can be activated nf number of times at different time steps
            for (uint neuron_id = 0; neuron_id < ntwk->m_N; neuron_id++) {
                uint nf = (ntwk->number_of_firing_times == 0) ? 
                                get_random_number(1, this->time_vec.size()) :
                                ntwk->number_of_firing_times;
                for (uint j = 0; j < nf; j++) {
                    _time_t rand_time = this->time_vec[get_random_number(0, this->time_vec.size())];
                    auto vec = n_rand_activate[rand_time];
                    vec.push_back(neuron_id);
                    n_rand_activate[rand_time] = vec;
                }
            }
}
*/

void Population::set_firing_time_for_random_time_network()
{
    for (auto ntwk : this->population_network) {
        if (ntwk->external_input == RANDOM_TIME) {
            std::map<_time_t, std::vector<uint>> n_rand_activate;

            uint time_index_between_patterns = 0;
            uint pattern_index = 0;
            _time_t starting_time = 0.0;
            _time_t time_between_patterns = min_time_required_for_simulation / ntwk->no_of_patterns;
            auto pattern_vec = this->p_rand_neuron_ids[ntwk->m_ntwk_id];

            for (uint j = 0; j < this->time_vec.size(); j++) {
                // insert p_rand_neurons ( patterns ) at specific time steps
                if (ntwk->enable_learning) {
                    if (fmod(this->time_vec[j], 100.0) == 0.0) {
                        starting_time = this->time_vec[j];
                        if (this->time_vec[j] != 0.0)
                            time_index_between_patterns = 1;
                    }

                    if ((this->time_vec[j] == 0.0) ||
                        (this->time_vec[j] == starting_time + (time_between_patterns * time_index_between_patterns) - 1)) {
                        if (!pattern_vec.empty()) {
//                            if ((int)(j > this->time_vec.size() / 2)) {
//                                auto pattern_tmp = pattern_vec[pattern_index];
//                                unsigned nseed = std::chrono::system_clock::now().time_since_epoch().count();
//                                std::shuffle(pattern_tmp.begin(), pattern_tmp.end(), std::default_random_engine(nseed));
//                                std::vector<uint> partial_pattern_vec(pattern_tmp.size()/2);
//                                partial_pattern_vec.reserve(pattern_tmp.size()/2);
//                                for (uint k = 0; k < pattern_tmp.size()/2; k++)
//                                    partial_pattern_vec[k] = pattern_tmp[k];
//                                n_rand_activate[this->time_vec[j]] = partial_pattern_vec;
//                            } else {
                                n_rand_activate[this->time_vec[j]] = pattern_vec[pattern_index];
//                            }

                            this->rand_pattern_times[this->time_vec[j]] = std::make_tuple(ntwk->m_ntwk_id, pattern_index);

                            // std::cout << __FUNCTION__ << "\ttime = " << this->time_vec[j] << " network id = " << ntwk->m_ntwk_id << " pattern index = " << pattern_index;

                            (pattern_index+1 < ntwk->no_of_patterns) ? pattern_index++ : pattern_index = 0;
                            time_index_between_patterns++;

                            continue;
                        }
                    }
                }

                std::vector<unsigned int> n_indices(ntwk->m_N);
                std::iota(n_indices.begin(), n_indices.end(), 0);
                unsigned nseed = std::chrono::system_clock::now().time_since_epoch().count();
                std::shuffle(n_indices.begin(), n_indices.end(), std::default_random_engine(nseed));

              // At every time step, activate 'n' number of pseudo-neurons
                std::vector<uint> n_indices_vec;
                uint rand = get_random_number(1, ntwk->m_N);
                for (uint i = 0; i < rand; i++) {
                    n_indices_vec.push_back(n_indices[i]);
                }
                n_rand_activate[this->time_vec[j]] = n_indices_vec;
            }

            n_rand_map[ntwk->m_ntwk_id] = n_rand_activate;
        }
    }
}

void Population::add_network(std::string filepath)
{
    logger->log("add_network");

    this->filepath = filepath;
    this->result_file_path = this->filepath + "/results/";

    Config config(this->filepath);

    add_network_internal_func(&config);
}

void Population::add_network(std::string filepath, std::string foldername)
{
    logger->log("add_network test framework");

    this->filepath = filepath;
    this->result_file_path = this->filepath + "/results/";

    Config config(this->filepath, foldername);

    add_network_internal_func(&config);
}

void Population::calculate_number_of_synapses()
{
    for (uint i = 0; i < m_N; i++) {
        uint sum = 0;
        for (uint j = 0; j < m_N; j++) {
            if (this->m_w_matrix[j][i] != 0)
                sum++;
        }
        //std::cout << get_network_name(i) << " = " << sum << std::endl;
    }

    for (uint i = 0; i < m_N; i++) {
        for (uint j = 0; j < m_N; j++) {
            //std::cout << this->m_w_matrix[i][j] << "\t";
        }
        //std::cout << std::endl;
    }
}

void Population::init_p_rand_neurons()
{
    for (auto ntwk : this->population_network) {
        if (ntwk->p_rand_no_of_neurons > 0) {
            uint percentage_p_rand_no_of_neurons = (ntwk->p_rand_no_of_neurons * ntwk->m_N) / 100;
            ntwk->intersecting_patterns ?
                init_intersecting_p_rand_neurons(percentage_p_rand_no_of_neurons,
                                                 ntwk->m_N,
                                                 ntwk->no_of_patterns,
                                                 ntwk->m_ntwk_id) :
                init_nonintersecting_p_rand_neurons(percentage_p_rand_no_of_neurons,
                                                 ntwk->m_N,
                                                 ntwk->no_of_patterns,
                                                 ntwk->start_from_row_index,
                                                 ntwk->m_ntwk_id);
        }
    }
}

void Population::init_intersecting_p_rand_neurons(uint p_rand_no_of_neurons,
                                                  uint N,
                                                  uint no_of_patterns,
                                                  uint ntwk_id)
{
    if (p_rand_no_of_neurons > 0) {
        std::vector<std::vector<uint>> tmp;
        tmp.reserve(0);

        std::vector<unsigned int> indices(N);
        std::iota(indices.begin(), indices.end(), 0);
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));

        std::vector<uint> pattern_vec;

        for (uint count = 0; count < no_of_patterns; count++) {
            uint prand_neuron = get_random_number(0, N);
            uint prand_neuron_start_index = prand_neuron;
            if (prand_neuron + p_rand_no_of_neurons > N)
                prand_neuron_start_index = prand_neuron - p_rand_no_of_neurons;

            for (uint index = prand_neuron_start_index;
                 index < (prand_neuron_start_index + p_rand_no_of_neurons);
                 index++) {
                pattern_vec.push_back(indices[index]);
            }

            tmp.push_back(pattern_vec);
            pattern_vec.clear();
        }
       
        this->p_rand_neuron_ids[ntwk_id] = tmp;
    }
}

void Population::init_nonintersecting_p_rand_neurons(uint p_rand_no_of_neurons,
                                                  uint N,
                                                  uint no_of_patterns,
                                                  uint start_from_row_index,
                                                  uint ntwk_id)
{
    if (p_rand_no_of_neurons > 0) {
//            if (ntwk->no_of_patterns > 0 &&
//                    (ntwk->m_N % ntwk->no_of_patterns != 0))
//                continue;

        std::vector<std::vector<uint>> tmp;
        tmp.reserve(0);

        std::vector<unsigned int> indices(N);
        std::iota(indices.begin(), indices.end(), 0);
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));

        std::vector<uint> pattern_vec;
        for (uint i = 0; i < N; i++) {
            if ((indices[i]+start_from_row_index) <
                    (start_from_row_index+N)) {
                pattern_vec.push_back(indices[i]);
            }
            if (i != 0 && (i+1) % p_rand_no_of_neurons == 0) {
                tmp.push_back(pattern_vec);
                pattern_vec.clear();
                if (tmp.size() >= no_of_patterns)
                    break;
            }
        }

        this->p_rand_neuron_ids[ntwk_id] = tmp;
    }
}

//void Population::set_p_rand_weight_matrix_without_learning(void)
//{
//    for (auto src_ntwk : population_network) {
//        for (auto dst_ntwk : population_network) {
//            if ((src_ntwk->m_type == PSEUDO_NEURON
//                    && dst_ntwk->m_type != PSEUDO_NEURON)
//                    && (src_ntwk->enable_learning == 1.0 
//                        && dst_ntwk->enable_learning == 1.0)) {
//                for (uint i = src_ntwk->start_from_row_index;
//                          i < src_ntwk->start_from_row_index + src_ntwk->m_N;
//                          i++) {
//                    for (uint j = dst_ntwk->start_from_row_index;
//                              j < dst_ntwk->start_from_row_index + dst_ntwk->m_N;
//                              j++) {
//                        if (m_w_matrix[i][j] != 0
//                            && is_neuron_in_p_rand(i, src_ntwk->m_ntwk_id)
//                            && is_neuron_in_p_rand(j, dst_ntwk->m_ntwk_id)) {
//                            m_w_matrix[i][j] += dst_ntwk->learning_rate;
//                        } else if (m_w_matrix[i][j] > 0) {
//                            m_w_matrix[i][j] -= dst_ntwk->unlearning_rate;
////                            m_w_matrix[i][j] = 0;
//                        }
//                    }
//                }
//            }
//        }
//    }
//}

void Population::add_network_internal_func(Config *config)
{
    this->total_time = config->ttime;
    this->time_step = (config->time_step == 0.0) ? 0.1 : config->time_step;
    this->bin_size = config->bin_size;

	get_time_vector(this->total_time, this->time_step);

    // Initialize neuron population one by one
    init_neuron_population(config);
    init_bins();

    calculate_network_size();

    init_weight_matrix();
    init_sstate_matrix();
    init_sfiring_matrix();
    init_nmatrices();
    init_p_rand_neurons();

    // create the random times for spiking population of neurons
    set_firing_time_for_random_time_network();

    // Create the network
    create_interneuronal_network_projections(config);

    // create another loop for recurrent collaterals
    create_recurrent_network_projections(config);

//    calculate_number_of_synapses();

//    set_p_rand_weight_matrix_without_learning();
}

void Population::set_neuron_firing_time(shared_ptr<Network> ntwk, uint i, _time_t firing)
{
    this->m_nf_matrix[ntwk->start_from_row_index + i] = firing;
}

void Population::set_neuron_state(shared_ptr<Network> ntwk, uint i, MTYPE type)
{
    this->m_n_matrix[ntwk->start_from_row_index + i] = type;
}

_time_t Population::get_noisy_delay(float del, float step)
{
    if (del == 0.0)
        return del;
    else if (step == 0.0)
        return del;
    return get_random_number(del - step, del + step);
}

void Population::activate_single_neuron(shared_ptr<Network> ntwk, uint n, uint i)
{
    logger->log("activate_single_neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]));

    uint cell_id = ntwk->start_from_row_index + i;

    if (this->m_n_matrix[cell_id] == ON) return;

    for (uint j = 0; j < this->m_N; j++) {
        if (this->m_w_matrix[j][cell_id] != 0) {
            this->m_sf_matrix[j][cell_id] = time_vec[n] + get_noisy_delay(ntwk->tau_del, this->del_step);
        }
    }
}

void Population::call_activate_single_neuron(shared_ptr<Network> ntwk, uint n, uint i)
{
    return activate_single_neuron(ntwk, n, i);
}

void Population::activate_single_random_synapse(shared_ptr<Network> ntwk, uint n, uint i)
{
    logger->log("activate_single_random_synapse: for neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]));

    uint cell_id = ntwk->start_from_row_index + i;

    if (this->m_n_matrix[cell_id] == ON) return;

    std::vector<uint> synapses;
    for (uint j = 0; j < this->m_N; j++) {
        if (this->m_w_matrix[j][cell_id] != 0) {
            synapses.push_back(j);
        }
    }

    uint rsyn = rand() % synapses.size();
    logger->log("activate_single_random_synapse: for neuron : " + std::to_string(i) + " : has #synapses : " + std::to_string(synapses.size()) + " : activate synapse : " + std::to_string(rsyn) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]));
    this->m_sf_matrix[rsyn][cell_id] = time_vec[n] + get_noisy_delay(ntwk->tau_del, this->del_step);
}

shared_ptr<Network> Population::get_network(uint neuron_id)
{
    for (auto ntwk : this->population_network) {
        uint size = ntwk->start_from_row_index + ntwk->m_N;
        if ((neuron_id >= ntwk->start_from_row_index) && (neuron_id < size)) {
            return ntwk;
        }
    }
    return nullptr;
}

string Population::get_network_name(uint neuron_id)
{
    for (auto ntwk : this->population_network) {
        uint size = ntwk->start_from_row_index + ntwk->m_N;
        if ((neuron_id >= ntwk->start_from_row_index) && (neuron_id < size)) {
            return ntwk->m_ntwk_name;
        }
    }
    return nullptr;
}



void Population::activate_synapses(uint n)
{
    logger->log("activate_synapses");

    std::vector<uint> neurons_to_compute;
    for (mmap::iterator itr = this->s_activate.begin(); itr != this->s_activate.end(); itr++) {
        if (itr->first == this->time_vec[n]) {
            std::map<uint, std::vector<uint>> neuron_map = itr->second;
            for (const auto &val : neuron_map) {
            // for(auto [neuron_id, syn_vec] : neuron_map) { // this is introduced in C++17. Commenting to avoid compiler warnings.
                auto neuron_id = val.first;
                auto syn_vec = val.second;
                for (uint i = 0; i < syn_vec.size(); i++) {
                    this->m_s_matrix[neuron_id][syn_vec[i]] = S_ON;
                    neurons_to_compute.push_back(syn_vec[i]);

                    logger->log("activate_synapses : switching on synapse : " + std::to_string(syn_vec[i]) + " for neuron : " + std::to_string(neuron_id)+ " : time : " + std::to_string(time_vec[n]));
                    std::string str = ("activate_synapses : switching on synapse : " + std::to_string(syn_vec[i]) + " for neuron : " + std::to_string(neuron_id)+ " : time : " + std::to_string(time_vec[n]));
//                    std::cout << str << std::endl;
                }
            }
        }
    }

    this->s_activate.erase(time_vec[n]);

    for (auto i : neurons_to_compute) {
        auto ntwk = get_network(i);
        if (ntwk) {
            threshold_block(ntwk, n, (i - ntwk->start_from_row_index), ON);
        }
    }
}

void Population::deactivate_synapses(uint n)
{
    logger->log("deactivate_synapses");

    if ( n == 0 ) return;

    /*
    std::cout << __FUNCTION__ << " : size of s map  = " << this->s_deactivate.size() << std::endl;

    std::vector<uint> neurons_to_compute;
    for (mmap::iterator itr = this->s_deactivate.begin(); itr != this->s_deactivate.end(); itr++) {
        if (itr->first == this->time_vec[n]) {
            std::map<uint, std::vector<uint>> neuron_map = itr->second;
            std::cout << __FUNCTION__ << " : size of neuron_map = " << neuron_map.size() << std::endl;
            for(auto [neuron_id, syn_vec] : neuron_map) {
                std::cout << __FUNCTION__ << " : neuron_id  = " << neuron_id << std::endl;
                std::cout << __FUNCTION__ << " : syn_vec.size = " << syn_vec.size() << std::endl;
                for (uint i = 0; i < syn_vec.size(); i++) {
                    this->m_s_matrix[neuron_id][syn_vec[i]] = S_OFF;
                    this->m_sf_matrix[neuron_id][syn_vec[i]] = 0.0f;
                    neurons_to_compute.push_back(syn_vec[i]);

                    logger->log("deactivate_synapses : switching off synapse : " + std::to_string(syn_vec[i]) + " for neuron : " + std::to_string(neuron_id)+ " : time : " + std::to_string(time_vec[n]));
                    std::string str = ("deactivate_synapses : switching off synapse : " + std::to_string(syn_vec[i]) + " for neuron : " + std::to_string(neuron_id)+ " : time : " + std::to_string(time_vec[n]));
                    std::cout << str << std::endl;
                }
            }
        }
    }

    this->s_deactivate.erase(time_vec[n]);

    std::cout << __FUNCTION__ << " : neurons_to_compute = " << neurons_to_compute.size() << std::endl;

    for (auto i : neurons_to_compute) {
        auto ntwk = get_network(i);
        if (ntwk) {
            std::cout << __FUNCTION__ << " : neuron_id - ntwk->start_from_row_index = " << i - ntwk->start_from_row_index << std::endl;
            threshold_block(ntwk, n, (i - ntwk->start_from_row_index), OFF);
        } else {
            std::cout << __FUNCTION__ << "ntwk = nullptr" << std::endl;
        }
    }
    */

    //std::cout << "SYNAPSE MATRIX = " << n << std::endl;
    /*
    for (uint i = 0; i < m_N; i++) {
        for (uint j = 0; j < m_N; j++)
            std::cout << this->m_s_matrix[i][j] << "\t";
        std::cout << std::endl;
    }
    */

    for (auto ntwk : this->population_network) {
        uint size = ntwk->start_from_row_index + ntwk->m_N;
        uint index = 0;
        for (uint i = ntwk->start_from_row_index; i < size; i++) {
            bool flag = 0;
            for (uint j = 0; j < m_N; j++) {
                if (this->m_s_matrix[i][j] == S_ON) {
                    if (std::abs((time_vec[n] - (this->m_sf_matrix[i][j] + ntwk->tau_dur))) < epsilon) {
                        this->m_s_matrix[i][j] = S_OFF;
                        this->m_sf_matrix[i][j] = 0.0f;
                        logger->log("deactivate_synapses : switching off synapse : " + std::to_string(j) + " for neuron : " + std::to_string(index)+ " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]));
                        std::string str = ("deactivate_synapses : switching off synapse : " + std::to_string(j) + " for neuron : " + std::to_string(index)+ " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]));
                        //std::cout << str << std::endl;

                        flag = 1;
                    }
                }
            }
            if (flag)
                threshold_block(ntwk, n, index, OFF);
            index++;
        }
    }
}

void Population::call_activate_synapses(uint n)
{
    activate_synapses(n);
}

void Population::call_deactivate_synapses(uint n)
{
    deactivate_synapses(n);
}

void Population::synaptic_block(shared_ptr<Network> ntwk, uint n, uint i)
{
    logger->log("synaptic_block");
    //std::cout << __FUNCTION__ << std::endl;
    
    deactivate_synapses(n);
    if (this->m_n_matrix[ntwk->start_from_row_index + i] == REF) {
        logger->log("synaptic_block: REF state neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
        std::string str = ("synaptic_block: REF state neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
        //std::cout << str << std::endl;
        threshold_block(ntwk, n, i, REF);
    } else {
        logger->log("synaptic_block: switch on neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
        std::string str = ("synaptic_block: switch on neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
        //std::cout << str << std::endl;
        threshold_block(ntwk, n, i, ON);
    }
    activate_synapses(n);
}

bool Population::is_neuron_in_p_rand(uint neuron_id, uint ntwk_id, _time_t n)
{
    auto patterns_vec = this->p_rand_neuron_ids[ntwk_id];

    for (uint p = 0; p < patterns_vec.size(); p++) {
        if (std::find(patterns_vec[p].begin(),
                      patterns_vec[p].end(),
                      neuron_id) != patterns_vec[p].end()) {
            return true;
        }
    }

    return false;
}

double Population::compute_weights(shared_ptr<Network> cell, uint n, uint i)
{
    logger->log("compute_weights with a given network for a particular cell");
    //std::cout << __FUNCTION__ << std::endl;

    double sum = 0.0f;
    uint cell_id = cell->start_from_row_index + i;

    if (this->m_n_matrix[cell_id] == ON) return 0.0f;

    for (auto ntwk : this->population_network) {
        uint size = ntwk->start_from_row_index + ntwk->m_N;
        for (uint j = ntwk->start_from_row_index; j < size; j++) {
            if ((this->m_w_matrix[j][cell_id] != 0) &&
                (this->m_s_matrix[j][cell_id] == S_ON)) {
                sum += this->m_w_matrix[j][cell_id];
            }
        }
    }
//    if (cell->m_ntwk_name == "bas")
//        std::cout << "network = " << cell->m_ntwk_name << " neuron ID = " << cell_id << " weight = " << sum << std::endl;

    // This is how weighted sum for external pseudo neurons is calculated
    if (cell->external_input == RANDOM_TIME)
        sum += cell->external_input_value;
    else
        sum += cell->external_input_vec[n];

    logger->log("weighted sum for neuron: " + std::to_string(i) + " in network : " + cell->m_ntwk_name + " at " + std::to_string(time_vec[n]) + " is : " + std::to_string(sum));
    std::string str = ("weighted sum for neuron: " + std::to_string(i) + " in network : " + cell->m_ntwk_name + " at " + std::to_string(time_vec[n]) + " is : " + std::to_string(sum));
    //std::cout << str << std::endl;

    return sum;
}

double Population::call_compute_weights(shared_ptr<Network> cell, uint n, uint i)
{
    return compute_weights(cell, n, i);
}

void Population::threshold_block(shared_ptr<Network> ntwk, uint n, uint i, MTYPE type)
{
    logger->log("threshold_block");
    //std::cout << __FUNCTION__ << std::endl;

    if (type == ON || type == REF) {
        if (this->m_n_matrix[ntwk->start_from_row_index + i] == ON ||
            this->m_n_matrix[ntwk->start_from_row_index + i] == REF) {
            return;
        }
    }

    double sum = compute_weights(ntwk, n, i);
    bool stop_learning = (int)(n > this->time_vec.size() / 2);

    if (type == ON) {
        //auto tmp = get_noisy_delay(ntwk->threshold, this->th_step);
        std::string str = ntwk->m_ntwk_name + " : time = " + std::to_string(this->time_vec[n]) + " : neuron_id = " + std::to_string(i) + " sum = " + std::to_string(sum) + "\n";
        logger->log(str);

        // If the neuron is not present in the selected p_rand, unlearn
        if (!stop_learning && ntwk->enable_learning && ntwk->p_rand_no_of_neurons > 0) {
            // If either of the neuron is not in the pattern and they share a projection, then unlearn
            if (!is_neuron_in_p_rand(i, ntwk->m_ntwk_id, this->time_vec[n]) ||
                !is_neuron_in_p_rand(this->cur_ntwk_neuron, this->cur_ntwk_id, this->time_vec[n])) {
                if (this->m_w_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron]
                                    [ntwk->start_from_row_index + i] != 0) {
                    this->m_w_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron]
                                    [ntwk->start_from_row_index + i]
                                                -= ntwk->unlearning_rate;
                }
            }
        }

        if (sum > ntwk->threshold) {
            //std::cout << __FUNCTION__ << " " << str << std::endl;
            this->m_n_matrix[ntwk->start_from_row_index + i] = ON;
            this->m_nf_matrix[ntwk->start_from_row_index + i] = time_vec[n];

            _time_t off_time = this->m_nf_matrix[ntwk->start_from_row_index + i] + ntwk->tau_ap + ntwk->tau_ref;
            _time_t ref_time = this->m_nf_matrix[ntwk->start_from_row_index + i] + ntwk->tau_ap;

            n_deactivate[off_time].push_back(ntwk->start_from_row_index + i);
            n_refractory[ref_time].push_back(ntwk->start_from_row_index + i);

            logger->log("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : off time = " + std::to_string(off_time)); 
            std::string str = ("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : off time = " + std::to_string(off_time)); 
//            std::cout << str << std::endl;

            logger->log("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : ref time = " + std::to_string(ref_time)); 
            str = ("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : ref time = " + std::to_string(ref_time));
//            std::cout << str << std::endl;

            // activate neuron i 's synapses including delay
            uint cur = ntwk->start_from_row_index + i;
            std::vector<uint> syn_vec;
            std::map<uint, std::vector<uint>> neuron_map;
            _time_t del = 0.0f;
            if (ntwk->del_step != 0.0f)
                del = get_noisy_delay(ntwk->tau_del, ntwk->del_step);
            else
                del = get_noisy_delay(ntwk->tau_del, this->del_step);
            //_time_t del = ntwk->tau_del;
            for (uint k = 0; k < this->m_N; k++) {
                if (this->m_w_matrix[cur][k] != 0) {
                    syn_vec.push_back(k);
                    if (ntwk->tau_del == 0) {
                        del = ntwk->tau_del;
                        logger->log("switch on synapses of neuron " + std::to_string(cur) + " synapse s = " + std::to_string(k) + "with no delay and duration = " + std::to_string(ntwk->tau_dur) + " at time " + std::to_string(time_vec[n]));
                        this->m_s_matrix[cur][k] = S_ON;
                    }
                    logger->log("switch on synapses of ntwk = " + ntwk->m_ntwk_name + " neuron " + std::to_string(cur) + " synapse s = " + std::to_string(k) + "with delay = " + std::to_string(del) + " and duration = " + std::to_string(ntwk->tau_dur) + " at time " + std::to_string(time_vec[n]));
                    std::string str = ("switch on synapses of ntwk = " + ntwk->m_ntwk_name + " neuron " + std::to_string(cur) + " synapse s = " + std::to_string(k) + " with delay = " + std::to_string(del) + " and duration = " + std::to_string(ntwk->tau_dur) + " at time " + std::to_string(time_vec[n]));
 //                   std::cout << str << std::endl;
                    this->m_sf_matrix[cur][k] = time_vec[n] + del;
                }
            }
            neuron_map[cur] = syn_vec;
            this->s_activate.insert(pair<_time_t, std::map<uint, std::vector<uint>>> ((time_vec[n] + del), neuron_map));
            this->s_deactivate.insert(pair<_time_t, std::map<uint, std::vector<uint>>> ((time_vec[n] + del + ntwk->tau_dur), neuron_map));

            logger->log("threshold_block: switching on neuron : " + std::to_string(ntwk->start_from_row_index + i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
            str = ("threshold_block: switching on neuron : " + std::to_string(ntwk->start_from_row_index + i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
            // std::cout << str << std::endl;

            // Logic to add learning : STDP rule
            if (!stop_learning && ntwk->enable_learning && ntwk->p_rand_no_of_neurons > 0) {
                // If both the neurons in the projection are in the pattern, and they share a projection, then learn
                if (is_neuron_in_p_rand(i, ntwk->m_ntwk_id, this->time_vec[n]) &&
                    is_neuron_in_p_rand(this->cur_ntwk_neuron, this->cur_ntwk_id, this->time_vec[n])) {
                    if (this->m_w_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron]
                                   [ntwk->start_from_row_index + i] != 0) { // if the 2 neurons are connected
                        // Check if the pre-synaptic neuron ( in this case, the external pseudo neuron is active
                        // and check if the pre-synaptic neuron is active within the last 20ms
                        bool cur_neuron_param_check = false;
                        if ((this->m_n_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron] == ON) &&
                                (this->m_nf_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron] <= time_vec[n]) &&
                                (this->m_nf_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron] > (time_vec[n]-20))) {
                            cur_neuron_param_check = true;
                        } else {
                            cur_neuron_param_check = false;
                        }

                        if (cur_neuron_param_check) {
                            this->m_w_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron]
                                       [ntwk->start_from_row_index + i]
                                                    += ntwk->learning_rate;
                        } else { // if cur neuron is not active in the last 20ms, unlearn the projection weight
                            this->m_w_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron]
                                       [ntwk->start_from_row_index + i]
                                                    -= ntwk->learning_rate;
                        }
//                        std::cout << "STDP: time : " << time_vec[n] << "\tneuron : " << (this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron)
//                                  << " - " << (ntwk->start_from_row_index + i) << " = "
//                                  << this->m_w_matrix[this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron][ntwk->start_from_row_index + i]
//                                  << std::endl;
                    }
                }
            }

//            if (stop_learning && ntwk->m_type != PSEUDO_NEURON) {
//                    std::cout << "time : " << time_vec[n] << "\tneuron : " << (this->cur_ntwk_start_from_row_index + this->cur_ntwk_neuron)
//                              << " - " << ntwk->start_from_row_index + i << std::endl;
//            }
        }
    } else {
        if (this->m_n_matrix[ntwk->start_from_row_index + i] != ON)
            return;

        //auto tmp = get_noisy_delay(ntwk->threshold, this->th_step);
        if (sum <= ntwk->threshold) {
            this->m_n_matrix[ntwk->start_from_row_index + i] = OFF;
            this->m_nf_matrix[ntwk->start_from_row_index + i] = 0.0f;

            logger->log("threshold_block: switching off neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
            std::string str = ("threshold_block: switching off neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
            //std::cout << str << std::endl;
        }
    }
}

void Population::deactivate_neurons(uint n)
{
    if ( n == 0 ) return;

    std::vector<uint> vec2 = n_refractory[time_vec[n]];
    for (uint i = 0; i < vec2.size(); i++) {
        this->m_n_matrix[vec2[i]] = REF;
        logger->log("deactivate_neurons : switching to REF neuron: " + std::to_string(vec2[i]) + " : time : " + std::to_string(time_vec[n]));
        std::string str = ("deactivate_neurons : switching to REF neuron: " + std::to_string(vec2[i]) + " : time : " + std::to_string(time_vec[n]));
        //std::cout << str << std::endl;
    }
    n_refractory.erase(time_vec[n]);

    std::vector<uint> vec1 = n_deactivate[time_vec[n]];
    for (uint i = 0; i < vec1.size(); i++) {
        this->m_n_matrix[vec1[i]] = OFF;
        this->m_nf_matrix[vec1[i]] = 0.0f;
        logger->log("deactivate_neurons : switching off neuron: " + std::to_string(vec1[i]) + " : time : " + std::to_string(time_vec[n]));
        std::string str = ("deactivate_neurons : switching off neuron: " + std::to_string(vec1[i]) + " : time : " + std::to_string(time_vec[n]));
//        std::cout << str << std::endl;
    }
    n_deactivate.erase(time_vec[n]);
}

void Population::call_deactivate_neurons(uint n)
{
    deactivate_neurons(n);
}

void Population::call_synaptic_block(shared_ptr<Network> ntwk, uint n, uint i)
{
    synaptic_block(ntwk, n, i);
}

void Population::call_threshold_block(shared_ptr<Network> ntwk, uint n, uint i, MTYPE type)
{
    threshold_block(ntwk, n, i, type);
}

/*
void Population::process_networks()
{
    logger->log("process_networks");

    srand(NULL);
    for (uint n = 0; n < time_vec.size(); n++) {
        deactivate_neurons(n);
        uint ntwk_id = rand() % this->population_network.size();
        uint k = (rand() % this->population_network[ntwk_id]->m_N);
        synaptic_block(this->population_network[ntwk_id], n, k);
        update_stats(this->time_vec[n]);
    }

    write_stats();

    logger->log("DONE");
}
*/

/*
double Population::get_random_pattern_index(uint ntwk_id,
                                            std::vector<uint> active_neurons_vec,
                                            uint &index)
{
    auto patterns_vec = this->p_rand_neuron_ids[ntwk_id];
    uint patterns_vec_size = 0;
    std::vector<uint> false_positive_count(patterns_vec.size(), 0);

    for (uint p = 0; p < patterns_vec.size(); p++) {
        bool found = true;
        for (uint active_neuron_index : active_neurons_vec) {
            if (std::find(patterns_vec[p].begin(),
                          patterns_vec[p].end(),
                          active_neuron_index) == patterns_vec[p].end()) {
                found = false;
                false_positive_count[p]++;
            }
        }
        if (found) { // If 100% recall
            index = p;
            return (patterns_vec[p].size() / patterns_vec[p].size())*100;
        }
        patterns_vec_size = patterns_vec[p].size();
    }

    // If no 100% recall, return that pattern which has the maximum number
    // of neurons that are recalled.
    
    auto it = std::max_element(false_positive_count.begin(),
                               false_positive_count.end());
    index = (it - false_positive_count.begin());
//    std::cout << "index = " << index << "    *it = " << *it << "    % = " << ((*it/patterns_vec_size)*100) << std::endl;
    return (*it/patterns_vec_size)*100;
}
*/


int Population::get_random_pattern_index(std::vector<std::vector<uint>> recall_count)
{
    for (uint i = 0; i < recall_count.size(); i++) {
        uint j = 0;
        for (; j < recall_count[i].size(); j++) {
            if (recall_count[i][j] == 0) {
                break;
            }
        }

        if (j == recall_count[i].size())
            return i;
    }

    return -1;
}

// Correlation
// (B . B') / sqrt(sumof(B) * sumof(B'))
// B = recalled pattern , B' = original pattern
// B and B' are binary : 1 = neuron is recalled, 0 = neuron is not recalled
// sumof(B') = always sizeof B'
// sumof(B) = how many neurons are recalled
// B . B' = dot product = how many neurons are recalled

void Population::get_recall_correlation(std::vector<std::vector<uint>> recall_count, uint pattern_size)
{
    std::ofstream out(this->result_file_path + "/recall_correlation.csv", std::ios_base::binary);
    std::vector<uint> expected_pattern(pattern_size, 1);
    for (uint i = 0; i < recall_count.size(); i++) {
        uint count = 0;
        uint bb = 0;
        for (uint j = 0; j < recall_count[i].size(); j++) {
            if (recall_count[i][j] > 0) {
                count++;
                bb += (recall_count[i][j] > 0) ? 1 * expected_pattern[j] : 0 * expected_pattern[j];
            }
        }
        auto corr = bb / sqrt(count * pattern_size);
        std::cout << " Pattern = " << (i+1) << " Recall correlation = " << corr << std::endl;

        if (out.good()) {
            out << (i+1) << "\t" << corr << "\n";
        }
    }

    out.close();
}

void Population::get_random_pattern_index(std::vector<std::vector<uint>> recall_count, _time_t n)
{
    std::ofstream out(this->result_file_path + "/recall_percentage.csv", std::ios_base::binary);

    for (uint i = 0; i < recall_count.size(); i++) {
        uint count = 0;
        for (uint j = 0; j < recall_count[i].size(); j++) {
            if (recall_count[i][j] != 0) {
                count++;
            }
        }
        auto percent = (count * 100.0 / recall_count[i].size());
        std::cout << "Time = " << (n+1) << " Pattern = " << (i+1) << " Recall quality = " << percent << "%"<< std::endl;

        if (out.good() && n == this->time_vec.back()) {
            out << (n+1) << "\t" << (i+1) << "\t" << percent  << "\n";
        }
    }

    out.close();
}

void Population::get_spurious_recall_count(_time_t n)
{
    uint count = 0;
    for (uint i = 0; i < this->m_spurious_recall_count.size(); i++) {
        if (this->m_spurious_recall_count[i]) count++;
    }

    auto percent = (double)count/this->m_spurious_recall_count.size();

    std::cout << "Time = " << (n+1) << " Spurious recall % = " << (percent * 100.0) << "%" << std::endl;
    std::cout << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}

void Population::process_networks()
{
    logger->log("process_networks");

    srand(0);

    for (uint n = 0; n < time_vec.size(); n++) {
        deactivate_neurons(n);

        bool flag = false;
        for (auto ntwk : this->population_network) {
            auto it = n_rand_map.find(ntwk->m_ntwk_id);
            if (it != n_rand_map.end()) {
//                if (ntwk->m_ntwk_name == "ext" || ntwk->m_ntwk_name == "ec")
//                    std::cout << "time = " << this->time_vec[n] << "    network = " << ntwk->m_ntwk_name << std::endl;
                auto n_rand_activate = it->second;
                auto it2 = n_rand_activate.find(time_vec[n]);
                if (it2 != n_rand_activate.end()) {
                    auto neuron_vec = it2->second;
                    for (uint j = 0; j < neuron_vec.size(); j++) {
                        cur_ntwk_neuron = neuron_vec[j];
                        cur_ntwk_id = ntwk->m_ntwk_id;
                        cur_ntwk_start_from_row_index = ntwk->start_from_row_index;
                        std::string str = ("process_networks: current external network id = " + std::to_string(ntwk->m_ntwk_id) + " : neuron id = " + std::to_string(neuron_vec[j]) + " : time = " + std::to_string(time_vec[n]));
                        logger->log(str);
//                        std::cout << str << std::endl;
                        synaptic_block(ntwk, n, neuron_vec[j]);
                        flag = true;
                    }
                }
            }
        }

        if (!flag) {
            deactivate_synapses(n);
            activate_synapses(n);
            flag = false;
        }

        update_stats(this->time_vec[n]);
    }

    // Find the final weights of ext->E neurons
//    uint pseudo_neuron_starting_row_index = 0;
//    uint no_of_pseudo_neurons = 0;
//    uint pseudo_neuron_ntwk_id = 0;
//    uint excitatory_population_starting_row_index = 0;
//    uint no_of_excitatory_neurons = 0;
//    uint excitatory_neuron_ntwk_id = 0;
//
//    for (auto ntwk : this->population_network) {
//        if (ntwk->enable_learning && ntwk->p_rand_no_of_neurons > 0) {
//            if (ntwk->m_type == PSEUDO_NEURON) {
//                pseudo_neuron_starting_row_index = ntwk->start_from_row_index;
//                no_of_pseudo_neurons = ntwk->m_N;
//                pseudo_neuron_ntwk_id = ntwk->m_ntwk_id;
//            } else if (ntwk->m_type == EXCITATORY) {
//                excitatory_population_starting_row_index = ntwk->start_from_row_index;
//                no_of_excitatory_neurons = ntwk->m_N;
//                excitatory_neuron_ntwk_id = ntwk->m_ntwk_id;
//            }
//        }
//    }

//    for (auto i = pseudo_neuron_starting_row_index;
//         i < (no_of_pseudo_neurons + pseudo_neuron_starting_row_index);
//         i++) {
//        if (is_neuron_in_p_rand((i-pseudo_neuron_starting_row_index), pseudo_neuron_ntwk_id, 0)) {
//            for (auto j = excitatory_population_starting_row_index;
//                 j < (no_of_excitatory_neurons + excitatory_population_starting_row_index);
//                 j++) {
//                if (is_neuron_in_p_rand((j-excitatory_population_starting_row_index), excitatory_neuron_ntwk_id, 0)) {
//                    if (this->m_w_matrix[i][j] > 0) {
//                        std::cout << "[" << i << "][" << j << "] = " << this->m_w_matrix[i][j] << std::endl;
//                    }
//                }
//            }
//        }
//    }

    // Recall metrics
    for (auto ntwk : this->population_network) {
        if (ntwk->m_type != PSEUDO_NEURON && ntwk->enable_learning == 1.0) {
            std::vector<uint> spurious_recall(ntwk->m_N, 0);
            std::vector<std::vector<uint>> recall_count;
            std::map<_time_t, int> pattern_recalled_at_time_t;
            auto patterns_vec = this->p_rand_neuron_ids[ntwk->m_ntwk_id];

            // initialize recall_count to zeroes
            for (uint p = 0; p < patterns_vec.size(); p++) {
                std::vector<uint> tmp_p(patterns_vec[p].size(), 0);
                recall_count.push_back(tmp_p);
            }

            std::vector<_time_t> random_pattern_times_after_900ms;
//            for (auto const & t : rand_pattern_times)
//                if (t.first >= 900) random_pattern_times_after_900ms.push_back(t.first);
            for (uint t = 899; t < this->time_vec.size();  t+=20)
                random_pattern_times_after_900ms.push_back(t);

            uint index = 1;
            for (uint i = 900; i < this->time_vec.size(); i++) {
                if (i == random_pattern_times_after_900ms[index]) {
                    get_random_pattern_index(recall_count, this->time_vec[i]);
                    get_spurious_recall_count(this->time_vec[i]);
                    index++;
                }

                auto active_neurons_vec = ntwk->active_neuron_stats[this->time_vec[i]];
                if (active_neurons_vec.empty()) continue;

                for (uint active_neuron_index : active_neurons_vec) {
                    bool found = false;
                    for (uint p = 0; p < patterns_vec.size(); p++) {
                        auto flag = std::find(patterns_vec[p].begin(), patterns_vec[p].end(), active_neuron_index);
                        if (flag != patterns_vec[p].end()) { // if the neuron in pattern is found to be active
                            uint pattern_index = flag - patterns_vec[p].begin();
                            if (pattern_index != patterns_vec[p].size()) {
                                recall_count[p][pattern_index] += 1;
                                found = true;
                            }
                        }
                    }

                    if (!found) { // if it is a spurious neuron recall
                        std::cout << "spurious recall of " << active_neuron_index << std::endl;
                        spurious_recall[active_neuron_index] += 1;
                    }
                }

                this->m_spurious_recall_count = spurious_recall;
            }

            get_recall_correlation(recall_count, patterns_vec[0].size());
        }
    }

    write_stats();

    logger->log("DONE");
}

int Population::get_bin_index(uint n)
{
    for (uint j = bin_size; j <= this->total_time; j += bin_size) {
        if (this->time_vec[n] < j) {
            return j;
        }
    }
    return -1;
}

void Population::update_stats(uint n)
{
    logger->log("update_stats");

    for (uint i = 0; i < this->population_network.size(); i++) { // i loop in init_bins
        double ac_stats = 0;
        double in_stats = 0;
        double ref_stats = 0;
        uint size = this->population_network[i]->start_from_row_index +
                                            this->population_network[i]->m_N;

        auto& tmap = this->bin_values[i]; // map of the ntwk -> this->population_network[i]
        std::vector<uint> a_n_stats;
        uint index = 0;

        for (uint j = this->population_network[i]->start_from_row_index; j < size; j++) { // j loop in init_bins

            if (this->m_n_matrix[j] == ON) {
                ++ac_stats;
                tmap[get_bin_index(n)].insert(j);
                a_n_stats.push_back(index);
            } else if (this->m_n_matrix[j] == OFF)
                ++in_stats;
            else if (this->m_n_matrix[j] == REF)
                 ++ref_stats;
            index++;
        }
        this->population_network[i]->active_neuron_stats.push_back(a_n_stats);

        logger->log("update_stats: network : " + this->population_network[i]->m_ntwk_name + " : time : " + std::to_string(this->time_vec[n]) + " : #active : " + std::to_string(ac_stats));

        this->population_network[i]->ac_stats.push_back(ac_stats / this->population_network[i]->m_N);
        this->population_network[i]->in_stats.push_back(in_stats / this->population_network[i]->m_N);
        this->population_network[i]->ref_stats.push_back(ref_stats / this->population_network[i]->m_N);
    }
}

void Population::call_update_stats(uint n)
{
    update_stats(n);
}

void Population::write_stats()
{
    logger->log("write_stats");

    std::ofstream out(this->result_file_path + "/ca_stats.csv", std::ios_base::binary);
    if (out.good()) {
        std::string heading = "time\t";
        for (auto cell : this->population_network) {
            heading += cell->m_ntwk_name + "_active\t" + cell->m_ntwk_name + "_inactive\t" +cell->m_ntwk_name + "_ref\t";
        }

        out << heading << "\n";
        for (uint i = 0; i < this->time_vec.size(); i++) {
            out << this->time_vec[i] << "\t";
            for (auto cell : this->population_network) {
                out << cell->ac_stats[i] << "\t";
                out << cell->in_stats[i] << "\t";
                out << cell->ref_stats[i] << "\t";
            }
            out << "\n";
        }
    }
    out.close();

    std::ofstream out2(this->result_file_path + "/ca_bin_stats.csv", std::ios_base::binary);
    if (out2.good()) {
        for (uint i = 0; i < this->bin_values.size(); i++) { // for each ntwk
            auto tmap = this->bin_values[i];
            out2 << "bins\t";
            for (auto it = tmap.begin(); it != tmap.end(); it++) {
                out2 << it->first << "\t";
            }
            out2 << "\n";
            out2 << this->population_network[i]->m_ntwk_name << "\t";

            for (auto it = tmap.begin(); it != tmap.end(); it++) {
                out2 << it->second.size() << "\t";
            }
            out2 << "\n";
        }
    }
    out2.close();

    for (auto cell : this->population_network) {
        std::ofstream out(this->result_file_path + "/"+cell->m_ntwk_name+".csv", std::ios_base::binary);
        if (out.good()) {
            for (uint i = 0; i < cell->active_neuron_stats.size(); i++) { // size of time_vec
                out << this->time_vec[i];
                uint index = 0;
                for (uint j = 0; j < cell->m_N; j++) { // for each neuron at each time step
                    if (!cell->active_neuron_stats[i].empty() && cell->active_neuron_stats[i][index] == j) {
                        out << "\t" << j;
                        index++;
                    } else {
                        out << "\t -1";
                    }
                }
                out << "\n";
            }
        }
        out.close();
    }

    std::ofstream out3(this->result_file_path + "/ca_p_rand_stats.csv", std::ios_base::binary);
    if (out3.good()) {
        for (auto ntwk : this->p_rand_neuron_ids) {
            if (ntwk.first > 101) continue;
            for (auto i : ntwk.second) {
                for (auto j : i) {
                    out3 << j << "\t";
                }
                out3 << "\n";
            }
        }
    }
    out3.close();

    std::ofstream out4(this->result_file_path + "/ca_weight_matrix.csv", std::ios_base::binary);
    if (out4.good()) {
        for (uint i = 0; i < this->m_w_matrix.size(); i++) {
            for (uint j = 0; j < this->m_w_matrix.size(); j++) {
                out4 << this->m_w_matrix[i][j] << "\t";
            }
            out4 << "\n";
        }
    }
    out4.close();

    std::ofstream out5(this->result_file_path + "/ca_could_be_spurious_recalled_neurons.csv", std::ios_base::binary);
    if (out5.good()) {
        for (auto ntwk : this->population_network) {
            if (ntwk->m_ntwk_id != 101) continue;

            std::vector<uint> sp_n;
            std::vector<uint> n_ids;
            std::vector<uint> prand_neurons;

            for (uint k = 0; k < ntwk->m_N; k++) n_ids.push_back(k);

            auto prand_ntwks = this->p_rand_neuron_ids[101];
            for (auto const &prand_ntwk : prand_ntwks) {
                prand_neurons.insert(prand_neurons.end(), prand_ntwk.begin(), prand_ntwk.end());
            }
            std::sort(prand_neurons.begin(), prand_neurons.end());

            std::set_difference(n_ids.begin(), n_ids.end(), prand_neurons.begin(), prand_neurons.end(), std::inserter(sp_n, sp_n.begin()));

            for (auto elem : sp_n)
                out5 << elem << "\t";
        }
    }
    out5.close();

    
    std::ofstream out6(this->result_file_path + "/ca_spurious_recalled_neurons.csv", std::ios_base::binary);
    if (out6.good()) {
        for (uint i = 0; i < this->m_spurious_recall_count.size(); i++) {
            out6 << i << "\t" << this->m_spurious_recall_count[i] << "\n";
        }
    }
    out6.close();
}
