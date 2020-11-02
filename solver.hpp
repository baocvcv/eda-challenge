#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__

#include "common.h"
#include "node.hpp"
#include "net.hpp"
#include "fpga.hpp"

#define COARSEN_THRESH 20

typedef std::pair<std::map<std::string, Node*>, std::vector<Net*>> Graph;

class Solver {
public:
    Solver(std::string input_dir, const Config &config):
            _input_dir(input_dir), _config(config) {
        _parse();
    }

    bool run() {
        _graphs.resize(1);
        _boundaries.clear();
        return _run(_graphs[0], _fpgas.size());
    }
    
    void save_results() { _save_results(); }

private:
    std::string _input_dir;
    Config _config;

    // fpga list
    std::vector<FPGA> _fpgas;

    // tree of partitioned graphs
    std::vector<Graph> _graphs;

    std::vector<std::vector<Net*>> _boundaries;

    void _parse() {
        std::map<std::string, Node*> nodes = _parse_nodes();
        std::vector<Net*> nets = _parse_nets(nodes);
        _fpgas = _parse_fpgas();
        _parse_constraints(nodes);
        _graphs.emplace_back(nodes, nets);
    }

    bool _run(Graph &graph, int num_parts) {
        if (num_parts == 1) // base
            return true;

        auto res = _bisect(graph, num_parts);
        auto bisection = res.first;
        auto boundary = res.second;
        if (bisection.size() == 2) {
            _graphs.push_back(bisection[0]);
            _graphs.push_back(bisection[1]);
            _boundaries.push_back(boundary);
            int len = _graphs.size();
            return _run(_graphs[len-2], num_parts / 2) && _run(_graphs[len-1], num_parts / 2);
        }
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
        std::ifstream in(_input_dir + "/" + net_def_file);
        if (!in.is_open()) {
            std::cerr << "Error reading " << net_def_file << std::endl;
        }

        // map from net list string to weigh
        std::map<std::string, int> net_raw;
        std::string net = "";
        std::string line;
        while (std::getline(in, line)) {
            int pos = line.find_first_of(' ');
            std::string node_name = line.substr(0, pos);
            line = line.substr(pos + 1);

            pos = line.find_first_of(' ');
            if (pos != std::string::npos) { // type "s"
                line = line.substr(pos + 1);
                int weight = std::stoi(line);

                if (net != "") {
                    net_raw[net] = weight; // duplicate net will be overwritten
                }
                net = node_name;
            } else { // type "l"
                net = net + " " + node_name;
            }
        }

        std::vector<Net*> nets;
        for (auto e: net_raw) {
            int weight = e.second;
            std::istringstream iss(e.first);

            std::string driver;
            iss >> driver;
            Node *driver_node = nodes[driver];
            Net *net = new Net(driver_node, weight);
            net->node_set.insert(driver_node);
            driver_node->net_set.insert(net);
            driver_node->drives.insert(net);

            std::string tmp;
            while (iss >> tmp) {
                Node *node = nodes[tmp];
                net->node_set.insert(node);
                node->net_set.insert(net);
            }
            nets.push_back(net);
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

    /**
     * @brief Bisect the graph into 2 parts
     * 
     * @param graph 
     * @param target_num_parts the number of parts the graph will eventually be divided into
     * @return std::pair<std::vector<Graph>, std::vector<Net*>> 
     */
    std::pair<std::vector<Graph>, std::vector<Net*>> _bisect(Graph &graph, int target_num_parts) {
        std::vector<Graph> graph_seq = _coarsen(graph, target_num_parts);
        Graph &last = graph_seq.back();
        _do_bisect(last);
        _uncoarsen(graph, graph_seq);

        // interpret result

        return std::make_pair(std::vector<Graph>(), std::vector<Net*>());
    }

    // coarsen
    std::vector<Graph> _coarsen(Graph &graph, int target_num_parts) {
        std::vector<Graph> res;
        Graph &cur_graph = graph;
        int size = cur_graph.first.size();
        while (size > COARSEN_THRESH) {
            // new graph
            std::map<std::string, Node*> nodes;
            std::vector<Net*> nets;

            // vertice address -> hypernode name
            std::map<Node*, std::string> matched_vertices;

            std::vector<Net*> original_nets(cur_graph.second);
            std::sort(original_nets.begin(), original_nets.end(), [](const Net* &a, const Net* &b) {
                return (a->cost > b->cost) || (a->cost == b->cost && a->node_set.size() < b->node_set.size());
            });

            // first pass
            for (auto net: original_nets) {
                std::vector<Node*> nodes_to_merge;

                // prepare
                std::set<int> board_constraints;
                for (auto _n: net->node_set) {
                    if (matched_vertices.count(_n)) { // already matched
                        continue;
                    }
                    if (_n->is_fixed) {
                        board_constraints.insert(_n->assigned_fpga.begin(), _n->assigned_fpga.end());
                    }
                    nodes_to_merge.push_back(_n);
                }
                if (board_constraints.size() > target_num_parts) {
                    // the number of nodes with different fixed allocations
                    // exceeds `target_num_parts`, so nodes cannot be merged
                    continue;
                }

                // merge the node
                Node *merged_node = new Node(nodes_to_merge);
                for (auto _n: nodes_to_merge) { // add to matched
                    matched_vertices[_n] = merged_node->node_name;
                }
            }

            // second pass


            
            // add to list
            res.push_back(std::make_pair(nodes, nets));
            cur_graph = res.back();
            size = cur_graph.first.size();
        }

        return res;
    }

    // split
    bool _do_bisect(Graph &graph) {
        return false;
    }

    void _uncoarsen(Graph &base, std::vector<Graph> &graph_seq) {

        // release memory
        for (auto g: graph_seq) {
            for (auto n: g.first) {
                delete n.second;
            }
            for (auto n: g.second) {
                delete n;
            }
        }
    }
};

#endif