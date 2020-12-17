#include "../include/log.hpp"

Log* Log::m_instance = nullptr;
std::ofstream Log::out;

Log* Log::getInstance() {
    if (!m_instance) {
        m_instance = new Log();
        out.open("/home/ashraya/Documents/Projects/CellularAutomata/CellularAutomata_Fast/results/log.txt", ios::out | ios::binary);
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
