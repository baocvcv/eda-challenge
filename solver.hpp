#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__

#include "common.h"
#include "node.hpp"
#include "net.hpp"
#include "fpga.hpp"

typedef std::pair<std::map<std::string, Node*>, std::vector<Net*>> Graph;

class Solver {
public:
    Solver(std::string input_dir, const Config &config):
            _input_dir(input_dir), _config(config) {
        _parse();
    }

    bool run() { return _run(); }
    
    void save_results() { _save_results(); }

private:
    std::string _input_dir;
    Config _config;

    std::vector<Graph> _graphs;
    std::vector<FPGA> _fpgas;

    void _parse() {
        std::map<std::string, Node*> nodes = _parse_nodes();
        std::vector<Net*> nets = _parse_nets(nodes);
        _fpgas = _parse_fpgas();
        _parse_constraints(nodes);
        _graphs.emplace_back(nodes, nets);
    }

    bool _run() {
        return false;
    }

    void _save_results() {

    }

    /**
     * @brief Parse nodes. Assume `node_num` appear sequentially.
     * 
     * @return std::vector<Node> 
     */
    std::map<std::string, Node*> _parse_nodes() {
        std::map<std::string, Node*> nodes;
        std::ifstream in(_input_dir + "/" + node_def_file);
        if (!in.is_open()) {
            std::cerr << "Error reading " << node_def_file << std::endl;
        }

        std::string line;
        while (std::getline(in, line)) {
            std::istringstream iss(line);

            // int node_num;
            std::string node_name;
            int resources[10];
            std::vector<std::string> clk_names;
            bool is_ff = false;

            iss >> node_name;
            // if (node_name.substr(0, 2) == "gp") // ignore gp nodes
            //     continue;
            // node_num = std::stoi(node_name.substr(1));
            for (int i = 0; i < 10; i++) {
                iss >> resources[i];
            }
            // parse attributes if any
            std::string rest_of_line;
            std::getline(iss, rest_of_line);
            rest_of_line = _trim(rest_of_line);
            if (rest_of_line.length() > 0) {
                rest_of_line = _trim(rest_of_line, "{}");
                std::istringstream iss2(rest_of_line);
                std::string tmp;
                while (iss2 >> tmp) {
                    tmp = _trim(tmp);
                    if (tmp == "ff") {
                        is_ff = true;
                    } else {
                        clk_names.push_back(tmp);
                    }
                }
            }

            // add node
            nodes[node_name] = new Node(node_name, resources, clk_names, is_ff);
        }
        return nodes;
    }

    std::vector<Net*> _parse_nets(std::map<std::string, Node*> nodes) {
        std::vector<Net*> nets;
        std::ifstream in(_input_dir + "/" + net_def_file);
        if (!in.is_open()) {
            std::cerr << "Error reading " << net_def_file << std::endl;
        }

        Net *net = nullptr;
        std::string line;
        while (std::getline(in, line)) {
            int pos = line.find_first_of(' ');
            std::string node_name = line.substr(0, pos);
            // int node_num = std::stoi(node_name.substr(1));
            line = line.substr(pos + 1);

            pos = line.find_first_of(' ');
            if (pos != std::string::npos) { // type "s"
                line = line.substr(pos + 1);
                int weight = std::stoi(line);
                nodes[node_name]->is_driver = true;

                if (net != nullptr) {
                    nets.push_back(net);
                }
                net = new Net(nodes[node_name], weight);
            } else { // type "l"
                net->node_set.insert(nodes[node_name]);
            }
        }
        return nets;
    }

    std::vector<FPGA> _parse_fpgas() {
        std::vector<FPGA> fpgas;
        std::ifstream in(_input_dir + "/" + fpga_res_file);
        if (!in.is_open()) {
            std::cerr << "Error reading " << fpga_res_file << std::endl;
        }

        std::string line;
        while (std::getline(in, line)) {
            int resources[10];
            std::vector<int> int_list;

            int pos = line.find_first_of(' ');
            line = line.substr(pos + 1);

            std::istringstream iss(line);
            for (int i = 0; i < 10; i++) {
                iss >> resources[i];
            }

            std::string rest_of_line;
            std::getline(iss, rest_of_line);
            rest_of_line = _trim(rest_of_line);
            if (rest_of_line.length() > 0) {
                rest_of_line = _trim(rest_of_line, "{}");
                std::istringstream iss2(rest_of_line);
                std::string tmp;
                while (iss2 >> tmp) {
                    int_list.push_back(std::stoi(tmp));
                }
            }
            fpgas.emplace_back(resources, int_list);
        }
        return fpgas;
    }

    void _parse_constraints(std::map<std::string, Node*> nodes) {
        std::ifstream in(_input_dir + "/" + constraint_file);
        if (!in.is_open()) {
            std::cerr << "Error reading " << constraint_file << std::endl;
        }

        std::string line;
        while (std::getline(in, line)) {
            int pos = line.find_first_of(':');
            std::string fpga_type = line.substr(0, pos);
            std::string node_list = line.substr(pos + 1);

            int type;
            sscanf(fpga_type.c_str(), "FPGA TYPE %d", &type);
            type -= 1; // fpga numbering start with 0

            std::istringstream iss(node_list);
            std::string node_name;
            while (iss >> node_name) {
                // int node_num = std::stoi(node_name.substr(1));
                _config.fixed_assignment[node_name] = type;
                nodes[node_name]->set_fixed(type);
            }
        }
    }

    static std::string _trim(const std::string &s, const std::string &delim = WHITE_SPACES) {
        std::string result = s;

        // left trim
        size_t start = result.find_first_not_of(delim);
        if (start != std::string::npos)
            result = result.substr(start);
        else
            return "";

        // right trim
        size_t end = result.find_last_not_of(delim);
        if (end != std::string::npos)
            result = result.substr(0, end + 1);

        return result;
    }

    // static std::vector<std::string> _split(const std::string &s, const char c) {
    //     std::vector<std::string> result;
    //     std::string tmp = s;
    //     while (tmp.size() > 0) {
    //         int pos = tmp.find_first_of(c);
    //         if (pos != )
    //         result.push_back(tmp.substr(pos));
    //         tmp = tmp.substr
    //     }
    // }
};

#endif