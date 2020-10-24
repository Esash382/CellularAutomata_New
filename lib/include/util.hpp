#pragma once

#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#ifdef _WIN32
#include <conio.h>
#endif

using std::string;
using namespace std;

typedef unsigned int uint;
typedef double _time_t;

typedef enum _MTYPE {
    OFF = 0,
    ON,
    REF
} MTYPE;

typedef enum _EXTERNAL_INPUT {
    PERIODIC = 0,
    RANDOM,
    POISSON,
    CONTINUOUS_RANDOM
} EXTERNAL_INPUT;

uint get_random_number(uint low, uint high);

const double epsilon = 0.000001;

class _time_ {
    public:
        static std::vector<_time_t> getInstance(_time_t total_time = 0.0, _time_t time_step = 0.0);

    private:
        _time_();
        std::vector<_time_t> set_linspace_time(_time_t total_time, _time_t time_step);

    public:
        static _time_t m_total_time;

    private:
        static _time_ *m_instance;
        static _time_t m_time_step;        
};
