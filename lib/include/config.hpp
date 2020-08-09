#pragma once

#include <map>
#include <vector>
#include <fstream>

#include "network.hpp"
#include "test.hpp"

class Config {
public:
    Config(std::string filepath);
    Config(std::string filepath, std::string foldername);

    std::map<string, float> read_file(std::string filename);
    std::vector<std::vector<std::string>> get_files();
    std::vector<std::string> get_unique_files();

    shared_ptr<Network> create_network(std::string name);

public:
    double time;
    double time_step;

private:
    void get_time();
    void read_network_file();

private:
    std::string filepath;
    std::vector<std::vector<std::string>> files;
    std::vector<std::string> unique_files;
};
