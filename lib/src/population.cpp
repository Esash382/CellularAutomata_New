#include "../external/population.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

using namespace std;

Population::Population()
    : th_step(10)
    , del_step(2)
{
	logger = Log::getInstance();

	this->time_vec.reserve(0);
	this->filepath = "";
	this->result_file_path = "";

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
    this->weighted_sum.reserve(this->m_N);

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

        uint dst_size = 0;
        uint dst_start_from_row_index = 0;

        for (auto n : this->population_network) {
            if (n->m_ntwk_name == cell[0]) {
                src_size = n->m_N;
                src_start_from_row_index = n->start_from_row_index;
            } else if (n->m_ntwk_name == cell[1]) {
                dst_size = n->m_N;
                dst_start_from_row_index = n->start_from_row_index;
            }

            if ((src_size != 0) && (dst_size != 0))
                break;
        }

//        std::cout << cell[0] << " = " << src_start_from_row_index << " : " << cell[1] << " = " << dst_start_from_row_index << std::endl;

        if (zsd != 0) {
            for (uint i = src_start_from_row_index; i < (src_start_from_row_index + src_size); i++) {
                uint index = 0;
                std::vector<unsigned int> indices(dst_size);
                std::iota(indices.begin(), indices.end(), dst_start_from_row_index);
                std::random_shuffle(indices.begin(), indices.end());

                for (uint k = 0; k < zsd; k++) {
                    uint j = indices[k];
                    assert(this->m_w_matrix[i][j] == 0);
                    this->m_w_matrix[i][j] = ksd;
                    this->m_s_matrix[i][j] = S_OFF;
                    this->m_sf_matrix[i][j] = 0.0f;
//                    std::cout << cell[0] << " -> " << cell[1] << " = " << zsd << " : " << index << " = " << k << " : this->m_w_matrix[" << i << "][" << j << "]" << this->m_w_matrix[i][j] << std::endl;
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
        
        if (zds != 0) {
            for (uint i = dst_start_from_row_index; i < (dst_start_from_row_index + dst_size); i++) {
                uint index = 0;
                std::vector<unsigned int> indices(src_size);
                std::iota(indices.begin(), indices.end(), src_start_from_row_index);
                std::random_shuffle(indices.begin(), indices.end());

                for (uint k = 0; k < zds; k++) {
                    uint j = indices[k];
                    assert(this->m_w_matrix[i][j] == 0);
                    this->m_w_matrix[i][j] = kds;
                    this->m_s_matrix[i][j] = S_OFF;
                    this->m_sf_matrix[i][j] = 0.0f;
//                    std::cout << cell[1] << " -> " << cell[0] << " = " << zds << " : " << index << " = " << k << " : this->m_w_matrix[" << i << "][" << j << "]" << this->m_w_matrix[i][j] << std::endl;
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
                std::random_shuffle(indices.begin(), indices.end());

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

void Population::add_network_internal_func(Config *config)
{
    this->total_time = config->time;
    this->time_step = (config->time_step == 0.0) ? 0.1 : config->time_step;

	get_time_vector(this->total_time, this->time_step);

    // Initialize neuron population one by one
    init_neuron_population(config);

    calculate_network_size();

    init_weight_matrix();
    init_sstate_matrix();
    init_sfiring_matrix();
    init_nmatrices();

    // Create the network
    create_interneuronal_network_projections(config);

    // create another loop for recurrent collaterals
    create_recurrent_network_projections(config);

    calculate_number_of_synapses();
}

void Population::set_neuron_firing_time(shared_ptr<Network> ntwk, uint i, _time_t firing)
{
    this->m_nf_matrix[ntwk->start_from_row_index + i] = firing;
}

void Population::set_neuron_state(shared_ptr<Network> ntwk, uint i, MTYPE type)
{
    this->m_n_matrix[ntwk->start_from_row_index + i] = type;
}

_time_t Population::get_noisy_delay(_time_t del, uint step)
{
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
                    //std::cout << str << std::endl;
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
//                        std::cout << str << std::endl;

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
    
    deactivate_synapses(n);
    activate_synapses(n);
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
}

double Population::compute_weights(shared_ptr<Network> cell, uint n, uint i)
{
    logger->log("compute_weights with a given network for a particular cell");

    double sum = 0.0f;
    uint cell_id = cell->start_from_row_index + i;

    if (this->m_n_matrix[cell_id] == ON) return 0.0f;

    for (auto ntwk : this->population_network) {
        uint size = ntwk->start_from_row_index + ntwk->m_N;
        for (uint j = ntwk->start_from_row_index; j < size; j++) {
            if (cell->m_ntwk_name == "hippocamposeptal" && ntwk->m_ntwk_name == "pyramidal") {
                //std::cout << n << " : " << ntwk->m_ntwk_name << " -> " << cell->m_ntwk_name << " = " << this->m_w_matrix[j][cell_id]
                  //          << " : " << ntwk->m_ntwk_name << " [" << j << "] = " << this->m_s_matrix[j][cell_id] << std::endl;
            }
            if ((this->m_w_matrix[j][cell_id] != 0) &&
                (this->m_s_matrix[j][cell_id] == S_ON)) {
                sum += this->m_w_matrix[j][cell_id];
            }
            // In case of pseudo neuronal synapses for external inputs
            if ((this->m_w_matrix[j][cell_id] != 0) &&
                (this->m_s_matrix[j][cell_id] != S_ON) &&
                (this->m_n_matrix[j] == ON)) {
                sum += this->m_w_matrix[j][cell_id];
            }
        }
    }

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

    if (type == ON || type == REF) {
        if (this->m_n_matrix[ntwk->start_from_row_index + i] == ON ||
            this->m_n_matrix[ntwk->start_from_row_index + i] == REF) {
            /*
            if (type == ON)
                std::cout << n << "\t" << ntwk->m_ntwk_name << "\t" << i << "\t ON" << "\t" << std::endl;
            else if (type == REF)
                std::cout << n << "\t" << ntwk->m_ntwk_name << "\t" << i << "\t REF" << "\t" << std::endl;
            */
            return;
        }
    }

    double sum = compute_weights(ntwk, n, i);

    /*
    if (type == ON)
        std::cout << n << "\t" << ntwk->m_ntwk_name << "\t" << i << "\t ON" << "\t" << sum << std::endl;
    else if (type == OFF)
        std::cout << n << "\t" << ntwk->m_ntwk_name << "\t" << i << "\t OFF" << "\t" << sum << std::endl;
    */

    if (type == ON) {
        if (sum > get_noisy_delay(ntwk->threshold, this->th_step)) {
            this->m_n_matrix[ntwk->start_from_row_index + i] = ON;
            this->m_nf_matrix[ntwk->start_from_row_index + i] = time_vec[n];

            _time_t off_time = this->m_nf_matrix[ntwk->start_from_row_index + i] + ntwk->tau_ap + ntwk->tau_ref;
            _time_t ref_time = this->m_nf_matrix[ntwk->start_from_row_index + i] + ntwk->tau_ap;

            n_deactivate[off_time].push_back(ntwk->start_from_row_index + i);
            n_refractory[ref_time].push_back(ntwk->start_from_row_index + i);

            logger->log("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : off time = " + std::to_string(off_time)); 
            std::string str = ("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : off time = " + std::to_string(off_time)); 
            //std::cout << str << std::endl;

            logger->log("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : ref time = " + std::to_string(ref_time)); 
            str = ("threshold_block: neuron: " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n]) + " : ref time = " + std::to_string(ref_time));
            //std::cout << str << std::endl;

            // activate neuron i 's synapses including delay
            uint cur = ntwk->start_from_row_index + i;
            std::vector<uint> syn_vec;
            std::map<uint, std::vector<uint>> neuron_map;
            _time_t del = get_noisy_delay(ntwk->tau_del, this->del_step);
            for (uint k = 0; k < this->m_N; k++) {
                if (this->m_w_matrix[cur][k] != 0) {
                    syn_vec.push_back(k);
                    if (ntwk->tau_del == 0) {
                        this->m_s_matrix[cur][k] = S_ON;
                    }
                    this->m_sf_matrix[cur][k] = time_vec[n] + del;
                }
            }
            neuron_map[cur] = syn_vec;
            this->s_activate.insert(pair<_time_t, std::map<uint, std::vector<uint>>> ((time_vec[n] + del), neuron_map));
            this->s_deactivate.insert(pair<_time_t, std::map<uint, std::vector<uint>>> ((time_vec[n] + del + ntwk->tau_dur), neuron_map));

            logger->log("threshold_block: switching on neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
            str = ("threshold_block: switching on neuron : " + std::to_string(i) + " : network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(time_vec[n])); 
            //std::cout << str << std::endl;
        }
    } else {
        if (sum <= get_noisy_delay(ntwk->threshold, this->th_step)) {
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

    std::vector<uint> vec1 = n_deactivate[time_vec[n]];
    for (uint i = 0; i < vec1.size(); i++) {
        this->m_n_matrix[vec1[i]] = OFF;
        this->m_nf_matrix[vec1[i]] = 0.0f;
        logger->log("deactivate_neurons : switching off neuron: " + std::to_string(i) + " : time : " + std::to_string(time_vec[n]));
        std::string str = ("deactivate_neurons : switching off neuron: " + std::to_string(i) + " : time : " + std::to_string(time_vec[n]));
        //std::cout << str << std::endl;
    }
    n_deactivate.erase(time_vec[n]);

    std::vector<uint> vec2 = n_refractory[time_vec[n]];
    for (uint i = 0; i < vec2.size(); i++) {
        this->m_n_matrix[vec2[i]] = REF;
        logger->log("deactivate_neurons : switching to REF neuron: " + std::to_string(i) + " : time : " + std::to_string(time_vec[n]));
        std::string str = ("deactivate_neurons : switching to REF neuron: " + std::to_string(i) + " : time : " + std::to_string(time_vec[n]));
        //std::cout << str << std::endl;
    }
    n_refractory.erase(time_vec[n]);
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

    /*
    uint index = 0;
    for (uint n = 0; n < time_vec.size(); n++) {
        deactivate_neurons(n);
        uint k = (rand() % this->population_network[index]->m_N);
        synaptic_block(this->population_network[index], n, k);
        update_stats(this->time_vec[n]);

        index++;
        if (index == this->population_network.size())
            index = 0;
    }
    */

    write_stats();

    logger->log("DONE");
}

void Population::update_stats(uint n)
{
    logger->log("update_stats");

    for (auto ntwk : this->population_network) {
        double ac_stats = 0;
        double in_stats = 0;
        double ref_stats = 0;
        uint size = ntwk->start_from_row_index + ntwk->m_N;
        for (uint i = ntwk->start_from_row_index; i < size; i++) {
            if (this->m_n_matrix[i] == ON)
                ++ac_stats;
            else if (this->m_n_matrix[i] == OFF)
                ++in_stats;
            else if (this->m_n_matrix[i] == REF)
                 ++ref_stats;
        }

        logger->log("update_stats: network : " + ntwk->m_ntwk_name + " : time : " + std::to_string(this->time_vec[n]) + " : #active : " + std::to_string(ac_stats));

//        std::cout << n << " : " << ntwk->m_ntwk_name << " : " << ac_stats << std::endl;

        ntwk->ac_stats.push_back(ac_stats / ntwk->m_N);
        ntwk->in_stats.push_back(in_stats / ntwk->m_N);
        ntwk->ref_stats.push_back(ref_stats / ntwk->m_N);
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
}
