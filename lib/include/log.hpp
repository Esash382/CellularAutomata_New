#pragma once

#include "util.hpp"

// #define DEBUG

class Log
{
    private:
        static Log *m_instance;
        static std::ofstream out;

    public:
        static Log * getInstance();
        static void freeInstance();
        static bool isLoggerActive();

        void close();

        void log(string msg);
};
