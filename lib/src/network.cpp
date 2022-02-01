#include "../include/network.hpp"

Network::Network(std::string name, uint ntwk_id, uint N, 
			NETWORK_TYPE type, float threshold,
			float k, _time_t tau_ap, _time_t tau_ref, uint z,
			_time_t tau_del, _time_t del_step, _time_t tau_dur,
            _time_t tau_osc, uint nf, EXTERNAL_INPUT e_in,
            float value, float step, uint enable_learning,
            float learning_rate, float unlearning_rate,
            uint p_rand_no_of_neurons)
{
	logger = Log::getInstance();

    this->m_ntwk_name = name;
	this->m_ntwk_id = ntwk_id;
	this->m_N = N;

	this->time_vec = _time_::getInstance(0.0, 0.0);

	this->m_Z = z;
	this->m_K_strength = k;
	this->m_type = type;

	this->external_input = e_in;
	this->external_input_value = value;

	this->threshold = threshold;
	this->tau_ap = tau_ap;
	this->tau_ref = tau_ref;
	this->tau_del = tau_del;
	this->del_step = del_step;
	this->tau_dur = tau_dur;
	this->tau_osc = tau_osc;
	this->last_osc_run_time = 0;

    this->number_of_firing_times = nf;

	uint time_vec_size = this->time_vec.size();
	this->ext_step = step;

    this->ac_stats.reserve(time_vec_size);
    this->in_stats.reserve(time_vec_size);
    this->ref_stats.reserve(time_vec_size);

	if (e_in == PERIODIC)
		external_input_vec = Spike_Generator::generate_periodic_spikes(time_vec_size, ext_step, this->external_input_value);
	else if (e_in == RANDOM)
		external_input_vec = Spike_Generator::generate_random_spikes(time_vec_size, this->external_input_value);
	else if (e_in == POISSON)
		external_input_vec = Spike_Generator::generate_poisson_spikes(time_vec_size, this->external_input_value);
    else if (e_in == CONTINUOUS_RANDOM)
        external_input_vec = Spike_Generator::generate_continuous_spikes(time_vec_size, ext_step, this->external_input_value);
    else {
        std::vector<float> vec(time_vec_size, 0);
        external_input_vec = vec;
    }

    this->enable_learning = enable_learning;
    this->learning_rate = learning_rate;
    this->unlearning_rate = unlearning_rate;
    this->p_rand_no_of_neurons = p_rand_no_of_neurons;
}
