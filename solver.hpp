#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include "common.h"
#include "node.hpp"
#include "net.hpp"
#include "fpga.hpp"
#include "graph.hpp"
#include "lib/patoh.h"

#define COARSEN_THRESH 0.05

class Solver {
public:
    Solver(std::string input_dir, const Config &config):
            _input_dir(input_dir), _config(config) {
        _parse();
    }

    Error check_input() { return _check_input(); }

    enum Method {
        Baseline,
        Multilevel
    };

    bool run(Method mode = Baseline) {
        srand(time(nullptr));

        std::set<int> available_fpgas;
        for (int i = 0; i < _fpgas.size(); i++) {
            available_fpgas.insert(i);
        }
        GraphResource res{_original_graph, available_fpgas};
        switch (mode) {
            case Baseline:
                _cur_result = _run_baseline(res);
                break;
            case Multilevel:
                _cur_result = _run_multilevel(res);
                break;
        }
        return _cur_result.size() == _original_graph.nodes.size();
    }

    void validate_and_save() {
        Result_T result_T;
        for (const auto &e: _cur_result) {
            result_T[e.second].push_back(e.first);
        }
        _cur_result_T = result_T;

        _evaluate_result();
        Error err = _validate_result();
        if (err.first) {
            std::cerr << "Result incorrect: " << err.second << std::endl;
        } else {
            _save_results(_cur_result_T, _cur_evaluation);
        }
    }

private:
    // partition result: node_name -> fpga
    typedef std::map<std::string, int> Result;
    typedef std::map<int, std::vector<std::string>> Result_T;

    // <Graph, set of available fpgas>
    typedef std::pair<Graph, std::set<int>> GraphResource;
    // tree of partitioned graphs
    // std::vector<GraphResource> _graphs;

    // input
    std::string _input_dir;
    Config _config;

    // fpga list
    std::vector<FPGA> _fpgas;

    // Original graph
    Graph _original_graph;

    // multisection result
    Result _cur_result;
    Result_T _cur_result_T;
    std::vector<FPGA> _cur_evaluation;

    void _parse() {
        auto nodes = _parse_nodes();
        auto nets = _parse_nets(nodes);
        auto needed_fpgas = _parse_constraints(nodes);
        _fpgas = _parse_fpgas(needed_fpgas);

        _original_graph.nodes = std::move(nodes);
        _original_graph.nets = std::move(nets);
    }

    Error _check_input() {
        // check resource constraint
        auto err = _check_fpga_resource_constraints();
        if (err.first)
            return err;

        // TODO: check other conditions?

        return {false, ""};
    }

    Result _run_baseline(GraphResource &gr) {
        //TODO: available for any number of fpgas
        Result res;

        std::stack<GraphResource> st;
        st.push(gr);
        while (!st.empty()) {
            GraphResource tmp = st.top(); st.pop();
            auto tmp_res = _split(tmp.first, tmp.second);

            // parse tmp_res
            for (int i = 0; i < tmp_res.first.size(); i++) {
                auto &nodes = tmp_res.first[i];
                auto &assignment = tmp_res.second[i];
                std::cout << assignment.size() << std::endl;
                if (nodes.size() > 0) {
                    if (assignment.size() > 1) { // need further partition
                        Graph g = tmp.first.get_subgraph(nodes);
                        st.push(GraphResource(g, assignment));
                    } else if (assignment.size() == 1) {
                        int ass = *assignment.begin();
                        for (const auto &n: nodes) // update result
                            res[n] = ass;
                    }
                }
            }
        }

        return res;
    }

    Result _run_patoh(GraphResource &gr) {
        // pack gr into format accepted by patoh
        auto &nodes = gr.first.nodes;
        auto &nets = gr.first.nets;
        int pin_total = 0;
        for (const auto& net: nets) {
            pin_total += net.node_set.size();
        }

        PaToH_Parameters args;
        PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_QUALITY);
        // TODO: can change this to gr.second.size()
        args._k = _config.fixed_fpgas.size();

        int _c = nodes.size();
        int _n = nets.size();
        int _nconst = 9;
        int useFixCells = 1;
        int *cwghts = new int[_c*_nconst];
        int *nwghts = new int[_n];
        int *xpins = new int[_n+1];
        int *pins = new int[pin_total];
        float *targetweights = new float[args._k];
        int *partvec = new int[_c];
        int *partweights = new int[args._k*_nconst];
        // TODO: int *cut =

        // run patoh_part

        // unpack result
    }

    Result _run_multilevel(GraphResource &gr) {
        auto &resources = gr.second;
        if (gr.first.nodes.size() == 0) { // base1
            return Result();
        }
        if (resources.size() == 1) { // base2
            // parse result
            int fpga = *resources.begin();
            Result res;
            for (const auto &e: gr.first.nodes) {
                res[e.first] = fpga;
                if (e.second.is_fixed() && e.second.assigned_fpga != fpga) {
                    // not satisfiable
                    return Result();
                }
            }
            return res;
        }

        auto res = _bisect(gr.first, resources);
        if (res.size() == 2) {
            Result r0 = _run_multilevel(res[0]);
            Result r1 = _run_multilevel(res[1]);
            if (r0.size() == 0)
                return r1;
            else if (r1.size() == 0)
                return r0;
            r0.insert(r1.begin(), r1.end());
            return r0;
        }
        return Result();
    }


    /*
     * Parsing helpers
     */

    std::map<std::string, Node> _parse_nodes() {
        std::map<std::string, Node> nodes;
        std::string _fn = _input_dir + "/" + node_def_file;
        std::ifstream in(_fn);
        if (!in.is_open()) {
            std::cerr << "Error reading " << _fn << std::endl;
        }

        std::string line;
        while (std::getline(in, line)) {
            std::istringstream iss(line);

            // int node_num;
            std::string node_name;
            int resources[10];
            std::set<std::string> clk_names;
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
                        clk_names.insert(tmp);
                    }
                }
            }

            // add node
            nodes[node_name] = Node(node_name, resources, clk_names, is_ff);
        }
        in.close();
        return nodes;
    }

    std::vector<Net> _parse_nets(std::map<std::string, Node> &nodes) {
        std::string _fn = _input_dir + "/" + net_def_file;
        std::ifstream in(_fn);
        if (!in.is_open()) {
            std::cerr << "Error reading " << _fn << std::endl;
        }

        // map from net list string to weigh
        std::map<std::string, int> net_raw;
        std::string net = "";
        int weight = 0; // weight of the `net`
        std::string line;
        while (std::getline(in, line)) {
            int pos = line.find_first_of(' ');
            std::string node_name = line.substr(0, pos);
            line = line.substr(pos + 1);

            pos = line.find_first_of(' ');
            if (pos != std::string::npos) { // type "s"
                line = line.substr(pos + 1);
                weight = std::stoi(line);

                if (net != "") {
                    // duplicate net will have their weights summed up
                    net_raw[net] += weight;
                }
                net = node_name;
            } else { // type "l"
                net = net + " " + node_name;
            }
        }
        if (net != "") { // the last net
            // duplicate net will have their weights summed up
            net_raw[net] += weight;
        }
        in.close();

        std::vector<Net> nets;
        for (auto &e: net_raw) {
            int idx = nets.size();
            int weight = e.second;

            std::istringstream iss(e.first);
            std::string driver;
            iss >> driver;
            nets.emplace_back(idx, driver, weight);
            Net &net = nets.back();

            Node &driver_node = nodes[driver];
            driver_node.net_set.insert(idx);
            driver_node.drives.insert(idx);

            std::string tmp;
            while (iss >> tmp) {
                net.node_set.insert(tmp);
                nodes[tmp].net_set.insert(idx);
            }
        }
        return nets;
    }

    std::vector<FPGA> _parse_fpgas(int max_num) {
        std::vector<FPGA> fpgas;
        std::string _fn = _input_dir + "/" + fpga_res_file;
        std::ifstream in(_fn);
        if (!in.is_open()) {
            std::cerr << "Error reading " << _fn << std::endl;
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
        in.close();
        return fpgas;
    }

    int _parse_constraints(std::map<std::string, Node> &nodes) {
        std::string _fn = _input_dir + "/" + constraint_file;
        std::ifstream in(_fn);
        if (!in.is_open()) {
            std::cerr << "Error reading " << _fn << std::endl;
        }

        std::string line;
        int max_num = 0;
        while (std::getline(in, line)) {
            int pos = line.find_first_of(':');
            std::string fpga_type = line.substr(0, pos);
            std::string node_list = line.substr(pos + 1);

            int type;
            sscanf(fpga_type.c_str(), "FPGA TYPE %d", &type);
            max_num = std::max(max_num, type);
            type -= 1; // fpga numbering start with 0

            std::istringstream iss(node_list);
            std::string node_name;
            while (iss >> node_name) {
                node_name = _trim(node_name);
                _config.fixed_assignment[node_name] = type;
                _config.fixed_fpgas.insert(type);
                nodes[node_name].set_fixed(type);
            }
        }
        in.close();
        for (unsigned i = 0; i < _fpgas.size(); i++) {
            if (_config.fixed_fpgas.count(i) == 0)
                _config.free_fpgas.insert(i);
        }

        return max_num;
    }


    /*
     * Input checking
     */
    
    /**
     * @brief Check if sum of required resources exceed the sum of available resources
     * 
     * @return Error 
     */
    Error _check_fpga_resource_constraints() {
        int required_resources[10], available_resources[10];
        for (int i = 0; i < 10; i++) {
            required_resources[i] = 0;
            available_resources[i] = 0;
        }

        for (const auto &e: _original_graph.nodes) {
            for (int i = 0; i < 10; i++) {
                required_resources[i] += e.second.resources[i];
            }
        }

        for (const auto &f: _fpgas) {
            for (int i = 0; i < 10; i++) {
                available_resources[i] += f.resources[i];
            }
        }

        for (int i = 0; i < 10; i++) {
            if (i == 1)
                continue;
            if (available_resources[i] < required_resources[i]) {
                std::ostringstream oss;
                oss << "Error: Not enough " << resource_names[i] << " resources." << std::endl;
                return {true, oss.str()};
            }
        }

        return {false, ""};
    }

    /*
     * Graph multi-section helpers
     */

    // <nodes_list, fpga_list>
    typedef std::pair<std::vector<std::set<std::string>>, std::vector<std::set<int>>> TmpResult;

    // Bisect the graph into 2 parts
    std::vector<GraphResource> _bisect(Graph &graph, const std::set<int> &available_fpgas) {
        std::vector<Graph> graph_seq = _coarsen(graph, available_fpgas);
        Graph &last = graph_seq.back();
        TmpResult result = _split(last, available_fpgas);
        return _uncoarsen(graph, graph_seq, result);
    }

    // coarsen
    std::vector<Graph> _coarsen(Graph &graph, std::set<int> available_fpgas) {
        std::vector<Graph> res;
        Graph cur_graph = graph; //TODO: better options?
        int size = cur_graph.nodes.size(); // num vertices
        float rate_of_change = 1.0;
        while (rate_of_change > COARSEN_THRESH) {
            // new graph
            Graph new_graph = cur_graph.make_copy();
            auto &nodes = new_graph.nodes;
            auto &original_nets = cur_graph.nets;
            std::sort(original_nets.begin(), original_nets.end(), [](const Net &a, const Net &b) {
                return (a.cost > b.cost) || (a.cost == b.cost && a.node_set.size() < b.node_set.size());
            });

            // vertex name -> hypernode name
            std::map<std::string, std::string> matched_vertices;

            // first pass
            // std::set<Net*> processed_nets;
            for (auto &net: original_nets) {
                std::vector<std::string> nodes_to_merge;

                if (net.node_set.size() == 0) {
                    std::cout << net.driver << std::endl;
                }

                // prepare
                // std::set<int> board_constraints;
                int board_constraint = -1;
                bool net_mergable = true;
                for (const auto &node_name: net.node_set) {
                    if (matched_vertices.count(node_name)) { // any node already matched
                        net_mergable = false;
                        break;
                    }

                    Node &n = nodes[node_name];
                    if (n.is_fixed()) {
                        if (board_constraint == -1) {
                            board_constraint = n.assigned_fpga;
                        } else if (board_constraint != n.assigned_fpga) { // fixed constraint in conflict
                            net_mergable = false;
                            break;
                        }
                    }
                    nodes_to_merge.push_back(node_name);
                }
                if (!net_mergable) {
                    // the number of nodes with different fixed allocations
                    // exceeds `target_num_parts`, so nodes cannot be merged
                    continue;
                }

                // merge the node
                std::string new_name = new_graph.merge(nodes_to_merge);
                for (auto _n: nodes_to_merge) { // add to matched
                    matched_vertices[_n] = new_name;
                }
            }

            // second pass
            // process the nodes that have not been merged
            for (auto &net: original_nets) {
                std::vector<std::string> nodes_to_merge;
                int board_constraint = -1;
                for (const auto &node_name: net.node_set) {
                    if (matched_vertices.count(node_name))
                        continue;

                    Node &n = nodes[node_name];
                    if (n.is_fixed()) {
                        if (board_constraint == -1) {
                            board_constraint = n.assigned_fpga;
                        } else if (board_constraint != n.assigned_fpga) {
                            // cannot satisfy fixed constraint
                            // merge current nodes
                            std::string new_name = new_graph.merge(nodes_to_merge);
                            for (auto _n: nodes_to_merge) { // add to matched
                                matched_vertices[_n] = new_name;
                            }
                            nodes_to_merge.clear();
                            board_constraint = n.assigned_fpga;
                        }
                    }
                    nodes_to_merge.push_back(node_name);
                }

                if (nodes_to_merge.size() > 0) {
                    std::string new_name = new_graph.merge(nodes_to_merge);
                    for (auto _n: nodes_to_merge) { // add to matched
                        matched_vertices[_n] = new_name;
                    }
                }
            }
            new_graph.update_net();

            if (cur_graph.nodes.size() != matched_vertices.size()) {
                for (auto &n: cur_graph.nodes) {
                    if (matched_vertices.count(n.first) == 0) {
                        auto name = new_graph.merge(std::vector<std::string>{n.first});
                        matched_vertices[n.first] = name;
                    }
                }
            }

            // std::set<std::string> all_nodes;
            // for (auto &node: new_graph.nodes)
            //     for (auto &n: node.second.node_set)
            //         all_nodes.insert(n);
            // for (auto &n: all_nodes) {
            //     if (cur_graph.nodes.count(n) == 0) {
            //         std::cout << n << " ";
            //     }
            // }
            // std::cout << std::endl;
            // std::cout << matched_vertices.size() << " " << all_nodes.size() << " " << cur_graph.nodes.size() << std::endl;
            
            // add to list
            res.push_back(new_graph);
            cur_graph = res.back();

            // update size
            int tmp = cur_graph.nodes.size();
            rate_of_change = 1.0 - (float)tmp / size;
            size = tmp;
        }

        return res;
    }

    // split
    TmpResult _split(Graph &graph, const std::set<int> &available_fpgas) {
        auto &nodes = graph.nodes;
        auto &nets = graph.nets;

        if (nodes.size() == 0) return TmpResult();

        std::set<std::string> selected_nodes;
        int resources_used[10]; 
        int resources_needed[10];
        for (int i = 0; i < 10; i++) { resources_used[i] = resources_needed[i] = 0; }
        for (const auto &e: graph.nodes) {
            for (int i = 0; i < 10; i++) {
                resources_needed[i] += e.second.resources[i];
            }
        }
        // a subset of `available_fpgas`, which are preferred in node allocation
        std::set<int> usable_fpgas = _set_intersection(available_fpgas, _config.fixed_fpgas);
        // fpgas already used
        std::set<int> assigned_fpgas;
        int total_available_resources[10];
        while (true) { // enlarge usable_fpgas if needed
            for (int i = 0; i < 10; i++) total_available_resources[i] = 0;
            for (auto f: usable_fpgas) {
                for (int i = 0; i < 10; i++) total_available_resources[i] += _fpgas[f].resources[i];
            }
            bool flag = true;
            for (int i = 0; i < 10; i++) {
                if (i != 1 && resources_needed[i] > total_available_resources[i]) {
                    flag = false;
                    break;
                }
            }
            if (!flag) { // not enough resources
                for (auto f: available_fpgas) {
                    if (usable_fpgas.count(f) == 0) {
                        usable_fpgas.insert(f);
                        break;
                    }
                }
                // cannot allocate more boards
                std::cerr << "Error: not enough resources" << std::endl;
                return TmpResult();
            } else {
                break;
            }
        }
        int num_fpga = usable_fpgas.size();
        int balanced_allocation[10];
        for (int i = 0; i < 10; i++)
            balanced_allocation[i] = resources_needed[i] / num_fpga;
        int available_resources[2][10]; // available resources for each split
        for (int i = 0; i < 20; i++) available_resources[i/10][i%10] = 0;

        /* priority queue */
        std::map<std::string, int> node2cost;
        std::vector<std::string> pq;

        auto comp = [&](const std::string &n1, const std::string &n2) { return node2cost[n1] > node2cost[n2]; };
        // push pq
        auto push_pq = [&](std::string n) {
            pq.push_back(n); std::push_heap(pq.begin(), pq.end(), comp);
        };
        // pop the top element from pq
        auto pop_pq = [&]() -> std::string {
            std::pop_heap(pq.begin(), pq.end(), comp); auto cur_node_name = pq.back(); pq.pop_back();
            node2cost.erase(cur_node_name);
            return cur_node_name;
        };
        // calculates the FM gain
        auto calc_cost = [&](const std::string &node_name) -> int {
            // if a node cannot be selected, then the cost is INT_MAX
            if (assigned_fpgas.size() >= usable_fpgas.size()/2 && nodes[node_name].is_fixed()) {
                return INT_MAX;
            }

            //TODO: is this correct?
            //TODO: more accurate cost function
            int cost = 0;
            for (auto net_num: nodes[node_name].net_set) {
                bool all_other_nodes_selected = true;
                bool all_nodes_unselected = true;
                for (const auto &name: nets[net_num].node_set) {
                    if (selected_nodes.count(name) > 0) {
                        // some name selected
                        all_nodes_unselected = false;
                    } else {
                        all_other_nodes_selected = false;
                    }
                }
                if (all_other_nodes_selected) { // about to become an inside net
                    cost -= nets[net_num].cost;
                } else if (all_nodes_unselected) { // about to become a cut net
                    cost += nets[net_num].cost;
                }
            }
            return cost;
        };
        // add edges to pq
        auto add_related_edges = [&](const std::vector<std::string> &nodes_to_add) {
            // check if there are enough resources
            int resources_d[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (const auto &n: nodes_to_add) {
                for (int i = 0 ; i < 10; i++) resources_d[i] += nodes[n].resources[i];
            }
            for (int i = 0; i < 10; i++) {
                if (i == 1) continue;
                if (resources_used[i] + resources_d[i] > available_resources[0][i])
                    return;
            }
            for (int i = 0; i < 10; i++) resources_used[i] += resources_d[i];
            // add nodes
            selected_nodes.insert(nodes_to_add.begin(), nodes_to_add.end());
            std::set<std::string> nodes_to_update;
            for (const auto &n: nodes_to_add) { // every newly added node
                // for (const auto &net_num: nodes[n].net_set) { // every related net
                //     for (const auto &nn: nets[net_num].node_set) { // every node in that net
                for (const auto &nn: graph.neighbors(n)) {
                    bool is_selected = selected_nodes.count(nn);
                    bool is_in_queue = node2cost.count(nn);
                    if (!is_selected && !is_in_queue) { // add new node
                        node2cost[nn] = calc_cost(nn);
                        push_pq(nn);
                    } else if (is_in_queue) { // update node already in queue
                        nodes_to_update.insert(nn);
                    }
                }
            }

            // update pq
            bool queue_updated = false; // whether elements in queue are updated
            //TODO: this is too slow
            // #pragma omp parallel for
            // for (const auto &nn: nodes_to_update) {
            //     int new_cost = calc_cost(nn);
            //     int& prev_cost = node2cost[nn];
            //     if (new_cost != prev_cost) {
            //         prev_cost = new_cost;
            //         queue_updated = true;
            //     }
            // }
            if (queue_updated)
                std::make_heap(pq.begin(), pq.end(), comp);
        };

        /* init pq */
        //TODO: start with a random vertex??
        // start with a random initial vertex with fixed assignment
        std::map<int, std::vector<std::string>> fixed_assignment;
        for (const auto &e: nodes) {
            const Node &n = e.second;
            if (n.is_fixed()) {
                fixed_assignment[n.assigned_fpga].push_back(n.node_name);
            }
        }

        std::vector<std::string> nodes_to_add;
        if (fixed_assignment.size() > 0) { // graph contains fixed assignments
            std::vector<int> fixes;
            for (const auto &e: fixed_assignment) {
                fixes.push_back(e.first);
            }
            int _start = fixes[rand() % fixes.size()];
            assigned_fpgas.insert(_start);
            nodes_to_add = fixed_assignment[_start];
            // add_related_edges(fixed_assignment[_start]);
        } else {
            // select a board
            int fpga = *available_fpgas.begin();
            assigned_fpgas.insert(fpga);
            // select random vertex
            int _start = rand() % nodes.size();
            auto it = nodes.begin();
            while (_start-- > 0) it++;
            // add_related_edges(std::vector<std::string>{it->first});
            nodes_to_add.push_back(it->first);
        }
        for (auto fpga: usable_fpgas) { // calcualte available resources
            int idx = 1 - assigned_fpgas.count(fpga);
            for (int i = 0; i < 10; i++)
                available_resources[idx][i] += _fpgas[fpga].resources[i];
        }
        add_related_edges(nodes_to_add);
        
        /* split */
        // select the best and update the queue
        // until the resource constraint and balanced conditions are met
        auto is_acceptable = [&]() -> bool {
            int num_fpga_used = assigned_fpgas.size();
            /* check resource constraints */
            int num_balanced = 0;
            for (int i = 0; i < 10; i++) {
                if (i == 1) continue;
                int res_other = resources_needed[i] - resources_used[i];
                // resource limit
                if (resources_used[i] > available_resources[0][i] || res_other > available_resources[1][i]) {
                    return false;
                }
                // resource balance
                if (balanced_allocation[i] != 0) {
                    //TODO: refine balanced conditoin
                    double ratio = fabs(1 - resources_used[i]/((double)balanced_allocation[i] * num_fpga_used)) * 100;
                    if (ratio <= _config.mean_percent)
                        num_balanced++;
                }
                // TODO: add constraints for int-list in --mean-mincut?
            }

            //TODO: add other constraints
            return num_balanced > 0;
            // return true;
        };
        std::set<std::string> node_visited;
        while (!is_acceptable() && pq.size() > 0) {
            std::string cur_node_name = pop_pq();
            node_visited.insert(cur_node_name);
            Node &cur_node = nodes[cur_node_name];
            if (cur_node.is_fixed() && assigned_fpgas.size() < usable_fpgas.size()/2) {
                // include all the nodes of the same assignment
                assigned_fpgas.insert(cur_node.assigned_fpga);
                add_related_edges(fixed_assignment[cur_node.assigned_fpga]);
                for (auto fpga: usable_fpgas) { // calcualte available resources
                    int idx = 1 - assigned_fpgas.count(fpga);
                    for (int i = 0; i < 10; i++)
                        available_resources[idx][i] += _fpgas[fpga].resources[i];
                }
            } else if (!cur_node.is_fixed()) {
                // add a single node
                add_related_edges(std::vector<std::string>{cur_node_name});
            }
            if (pq.size() <= 0) { // queue empty, try to restart
                std::vector<std::string> vt;
                for (const auto &e: nodes) {
                    if (selected_nodes.count(e.first) == 0 && node_visited.count(e.first) == 0) {
                        vt.push_back(e.first);
                    }
                }
                if (vt.size() > 0) {
                    std::string n = vt[rand() % vt.size()];
                    node2cost[n] = calc_cost(n);
                    push_pq(n);
                }
            }
        }

        /* parse result */
        // assign fpga
        // assigned_fpga = assign_fpga();
        std::set<int> complement;
        for (auto i: available_fpgas) {
            if (assigned_fpgas.count(i) == 0)
                complement.insert(i);
        }
        std::set<std::string> nodes_partition[2];
        for (auto &e: nodes) {
            const auto &name = e.first;
            int idx = 1 - selected_nodes.count(name);
            nodes_partition[idx].insert(name);
        }

        return {
            std::vector<std::set<std::string>>{nodes_partition[0], nodes_partition[1]},
            std::vector<std::set<int>>{assigned_fpgas, complement}
        };
    }

    std::vector<GraphResource> _uncoarsen(Graph &base, std::vector<Graph> &graph_seq, TmpResult result) {
        for (int i = graph_seq.size() - 1; i >= 0; i--) {
            if (i > 0)
                result = _level_up(graph_seq[i], graph_seq[i-1], result);
            else
                result = _level_up(graph_seq[i], base, result);
        }

        std::vector<GraphResource> grs;
        auto &nodes_list = result.first;
        auto &fpga_list = result.second;
        for (int i = 0; i < nodes_list.size(); i++) {
            auto sub_graph = base.get_subgraph(nodes_list[i]);
            grs.emplace_back(std::move(sub_graph), fpga_list[i]);
        }
        return grs;
    }

    TmpResult _level_up(Graph &child, Graph &parent, const TmpResult &result) {
        // project up
        auto &nodes_list = result.first;
        auto &fpga_list = result.second;

        assert(nodes_list.size() == 2);
        
        auto nodes1 = child.unwrap(nodes_list[0]);
        auto nodes2 = child.unwrap(nodes_list[1]);
        
        // adjustments
        //TODO

        return {
            std::vector<std::set<std::string>>{nodes1, nodes2},
            fpga_list
        };
    }

    std::set<int> _set_intersection(const std::set<int> &s1, const std::set<int> &s2) {
        std::set<int> res;
        for (auto i: s1) {
            if (s2.count(i) > 0)
                res.insert(i);
        }
        return res;
    }


    /*
     * Validate, evaluate and save result
     */
    Error _validate_result() {

        /* basic evaluation */
        // fix assignment condition is met
        for (const auto &e: _config.fixed_assignment) {
            if (_cur_result[e.first] != e.second) {
                std::ostringstream oss;
                oss << "Fixed assignment " << e.first << " to " << e.second << " not satisfied.";
                return {true, oss.str()};
            }
        }

        // every node is accounted for and only appear once
        for (const auto &e: _original_graph.nodes) {
            if (_cur_result.count(e.first) == 0) {
                std::ostringstream oss;
                oss << "Node " << e.first << " not allocated to a board.";
                return {true, oss.str()};
            }
        }

        // fpga resource constraint
        for (const auto &e: _cur_result_T) {
            int board = e.first;
            auto &nodes = e.second;
            int resources_used[10];
            for (int i = 0; i < 10; i++) resources_used[i] = 0;
            for (const auto &name: nodes) {
                Node &n = _original_graph.nodes[name];
                for (int i = 0; i < 10; i++)
                    resources_used[i] += n.resources[i];
            }
            FPGA &fpga = _fpgas[board];
            for (int i = 0; i < 10; i++) {
                if (i == 1) continue;
                if (resources_used[i] > fpga.resources[i]) {
                    std::ostringstream oss;
                    oss << "Resource limit exceeded on FPGA TYPE " << board << ". ";
                    oss << resource_names[i] << " used " << resources_used[i];
                    oss << " exceeds maximum allowed resource of " << fpga.resources[i];
                    return {true, oss.str()};
                }
            }
        }

        return {false, ""};
    }

    void _evaluate_result() {
        //TODO: only calculate stats for used fpga
        std::vector<FPGA> result(_fpgas.size(), FPGA(_fpgas.size()));
        // sum up resources
        for (const auto &e: _cur_result) {
            auto &n_name = e.first;
            auto &n = _original_graph.nodes[n_name];
            auto fpga = e.second;
            for (int i = 0; i < 10; i++)
                result[fpga].resources[i] += n.resources[i];
        }
        // calculate int list
        for (const auto &net: _original_graph.nets) {
            int driver_board = _cur_result[net.driver];

            bool *contains_node_in = new bool[result.size()];
            for (unsigned i = 0; i < result.size(); i++) // init
                contains_node_in[i] = false;
            for (const auto &n_name: net.node_set)
                contains_node_in[_cur_result[n_name]] = true;
            contains_node_in[driver_board] = false;

            for (unsigned i = 0; i < result.size(); i++) {
                if (contains_node_in[i]) {
                    result[driver_board].int_list[i] += net.cost;
                    result[i].int_list[driver_board] += net.cost;
                }
            }
        }
        _cur_evaluation = result;
    }

    void _save_results(const Result_T &allocation_result, const std::vector<FPGA> &resource_result) {
        // save allocation results
        std::string _fn = _input_dir + "/" + result_file;
        std::ofstream out(_fn);
        if (!out.is_open()) {
            std::cerr << "Cannot open " << _fn << " to write" << std::endl;
            return;
        } else {
            std::cout << "Saving allocation result to " << _fn << std::endl;
        }

        for (const auto &e: allocation_result) {
            int fpga = e.first;
            const auto &ns = e.second;
            out << "FPGA TYPE " << fpga+1 << " : ";
            for (int i = 0; i < ns.size(); i++) {
                out << ns[i];
                out << (i % 20 == 19 ? "\n" : " ");
            }
            if (ns.size() % 20 != 19)
                out << std::endl;
        }
        out.close();

        // save resource results
        _fn = _input_dir + "/" + eval_file;
        out = std::ofstream(_fn);
        if (!out.is_open()) {
            std::cerr << "Cannot open " << _fn << " to write" << std::endl;
            return;
        } else {
            std::cout << "Saving evaluation result to " << _fn << std::endl;
        }

        for (unsigned i = 0; i < resource_result.size(); i++) {
            const FPGA &fpga = resource_result[i];
            out << "FPGA TYPE " << i+1 << " : ";
            for (int i = 0; i < 10; i++)
                out << fpga.resources[i] << " ";
            // if (result.size() > 2) {
                out << "{ ";
                for (auto t: fpga.int_list)
                    out << t << " ";
                out << "}";
            // }
            out << std::endl;
        }
        out.close();
    }


    /*
     * Utility helpers
     */
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
};

#endif