#include "../include/log.hpp"
#include <unistd.h>

Log* Log::m_instance = nullptr;
std::ofstream Log::out;

Log* Log::getInstance()
{
    std::string username = getlogin();
    if (!m_instance) {
        m_instance = new Log();
        out.open("/home/" + username + "/Documents/Notes/CellularAutomata_Fast/results/log.txt", ios::out | ios::binary);
    }

    return m_instance;
}

void Log::freeInstance() {
    delete m_instance;
    m_instance = nullptr;
}

bool Log::isLoggerActive() {
    if (!m_instance)
        return false;

    return true;
}

void Log::close() {
    out.close();
    freeInstance();
}

void Log::log(string msg) {
    #ifdef DEBUG
        out << msg << endl;
    #endif
}
