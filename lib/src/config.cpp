#include "../include/json.hpp"
#include "../include/config.hpp"

#include <string.h>

using json = nlohmann::json;

Config::Config(std::string filepath)
{
    this->filepath = filepath + "/config";
    get_time();
    read_network_file();
}

Config::Config(std::string filepath, std::string foldername)
{
    this->filepath = filepath + "/config/" + foldername;
    get_time();
    read_network_file();
}

void Config::get_time()
{
	Log* logger = Log::getInstance();

    string file = this->filepath + "/time.json";
    std::map<std::string, float> con;
    std::ifstream ifs(file);

    if (ifs.good()) {
        json j;
        ifs >> j;
        for (auto it: j.items()) {
            con[it.key()] = it.value();
        }
    } else {
        logger->log("Error: json file time.json is not readable.");
    }

    ifs.close();

    this->ttime = con["time"];
    this->time_step = con["step"];
    this->bin_size = con["binsize"];
}

std::map<string, float> Config::read_file(std::string filename)
{
	Log* logger = Log::getInstance();

    string file = this->filepath + "/" + filename + ".json";
    std::map<std::string, float> con;
    std::ifstream ifs(file);

    if (ifs.good()) {
        json j;
        ifs >> j;
        for (auto it: j.items()) {
            con[it.key()] = it.value();
        }
    } else {
        logger->log("Error: json file "+filename+".json is not readable.");
    }

    ifs.close();

    return con;
}

void Config::read_network_file()
{
	std::string path = filepath + "/networks.dat";
    std::ifstream infile(path);

    std::string line;

    while (std::getline(infile, line))
    {
        std::vector<string> file;
        char* line_chr = const_cast<char*>(line.c_str());
        char* s_file = strtok (line_chr, " ");
        while (s_file != NULL) {
            file.push_back(s_file);
            this->unique_files.push_back(s_file);
            s_file = strtok (NULL, " ");
        }
        this->files.push_back(file);
    }

    std::sort( this->unique_files.begin(), this->unique_files.end() );
    this->unique_files.erase( unique( this->unique_files.begin(), this->unique_files.end() ), this->unique_files.end() );

    return;
}

std::vector<std::vector<std::string>> Config::get_files()
{
    return this->files;
}

std::vector<std::string> Config::get_unique_files()
{
    return this->unique_files;
}

shared_ptr<Network> Config::create_network(std::string name)
{
	Log* logger = Log::getInstance();
	std::map<string, float> con = read_file(name);
    if (con.empty()) {
        return NULL;
    }
    logger->log("read_file = " + name);

    float ext_val = 0.0f;
    uint num_ext_input = con["num_ext_input"];

    if (num_ext_input == 0) {
        return std::make_shared<Network>(name, uint(con["id"]), uint(con["N"]), 
                        NETWORK_TYPE(con["type"]),
                        con["threshold"], con["k"], con["ap"],
                        con["ref"], uint(con["z"]), con["del"], 
                        con["dur"], con["osc"], con["nf"],
                        PERIODIC, 0.0f, con["ext_step"]);
    }

    for (uint i = 1; i <= num_ext_input; i++) {
        ext_val += con["ext_val_"+std::to_string(i)];
    }

    return std::make_shared<Network>(name, uint(con["id"]), uint(con["N"]), 
                        NETWORK_TYPE(con["type"]),
                        con["threshold"], con["k"], con["ap"],
                        con["ref"], uint(con["z"]), con["del"], 
                        con["dur"], con["osc"], con["nf"],
                        EXTERNAL_INPUT(con["ext_type"]), ext_val, con["ext_step"]);
}
