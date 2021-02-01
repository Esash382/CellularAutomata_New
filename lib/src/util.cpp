#include "../include/util.hpp"

uint get_random_number(uint low, uint high)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint> d(low, high-1);

    return d(gen);
}

_time_t get_random_real_number(float low, float high)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> d(low, high-1);

    return d(gen);
}

_time_* _time_::m_instance = 0;
_time_t _time_::m_total_time = 0.0;
_time_t _time_::m_time_step = 0.0;

std::vector<_time_t> _time_::getInstance(_time_t ttime, _time_t tstep)
{
    if (m_instance == 0) {
        m_instance = new _time_();
        m_total_time = ttime;
        m_time_step = tstep;
    }

    return m_instance->set_linspace_time(m_total_time, m_time_step);
}

_time_::_time_() { }

std::vector<_time_t> _time_::set_linspace_time(_time_t total_time, _time_t time_step)
{
	std::vector<_time_t> time_vec;
    for (_time_t i = 0.0; i < total_time; i += time_step) {
		time_vec.push_back(i);
    }
	
	return time_vec;
}
