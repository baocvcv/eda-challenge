#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__

#include "common.h"
#include "node.hpp"
#include "net.hpp"

struct Graph {
    // map from node name to `Node`
    std::map<std::string, Node> nodes;

    // array of Net
    std::vector<Net> nets;

    Graph() {}

    Graph(const std::map<std::string, Node> &_nodes, const std::vector<Net> &_nets):
        nodes(_nodes), nets(_nets) {}

    // Convert to higher level graph
    Graph make_copy() const {
        Graph res(*this);
        for (auto &n: res.nodes) {
            n.second.node_set.clear();
            n.second.node_set.insert(n.first);
            n.second.net_set.clear();
            n.second.is_leaf = false;
        }
        return res;
    }

    // generate a sub graph containing only nodes in `node_set`
    Graph get_subgraph(const std::set<std::string> &node_set) const {
        // select nodes
        std::map<std::string, Node> nodes_n;
        for (const auto &e: nodes) {
            if (node_set.count(e.first) > 0) {
                nodes_n[e.first] = e.second;
            }
        }

        // shrink nets
        std::vector<Net> nets_n(nets.begin(), nets.end());
        auto it = nets_n.begin();
        while (it != nets_n.end()) {
            std::set<std::string> tmp;
            for (const auto &n: it->node_set)
                if (node_set.count(n) > 0)
                    tmp.insert(n);
            if (tmp.size() > 1) {
                it->node_set = tmp;
                if (node_set.count(it->driver) == 0)
                    it->driver = "";
                it++;
            } else {
                it = nets_n.erase(it);
            }
        }

        if (nets_n.size() != nets.size()) { // numbering has changed
            for (auto &e: nodes_n) { // reset
                e.second.net_set.clear();
                e.second.drives.clear();
            }
            for (unsigned i = 0; i < nets_n.size(); i++) {
                for (const auto &name: nets_n[i].node_set) {
                    Node &n = nodes_n[name];
                    n.net_set.insert(i);
                }
                const auto &driver = nets_n[i].driver;
                if (driver != "")
                    nodes_n[driver].drives.insert(i);
            }
        }
        return Graph(nodes_n, nets_n);
    }

    // unwrap to get nodes in parent graph
    std::set<std::string> unwrap(const std::set<std::string> &node_set) {
        std::set<std::string> result;
        for (const auto &node_name: node_set) {
            const auto &node = nodes[node_name];
            result.insert(node.node_set.begin(), node.node_set.end());
        }
        return result;
    }

    // return the border nodes corresponding to a split represented by `node_set`
    // std::pair<std::set<std::string>, std::set<std::string>> get_border_nodes(const std::set<std::string> &node_set) const {
    //     //TODO:
    // }

    std::string merge(const std::vector<std::string> &nodes_to_merge) {
        if (nodes_to_merge.size() == 0) {
            std::cerr << "No nodes to be merged!" << std::endl;
            return "";
        }

        Node &merged_node = nodes[nodes_to_merge[0]];
        for (unsigned int i = 1; i < nodes_to_merge.size(); i++) {
            const std::string &name = nodes_to_merge[i];
            merged_node.merge_with(nodes[name]);
            nodes.erase(name);
        }
        return merged_node.node_name;
    }

    void merge(std::string src, std::string target) {
        nodes[target].merge_with(nodes[src]);
        nodes.erase(src);
    }

    void update_net(std::map<std::string, std::string> mapping) {
        // update nets' node list
        auto it = nets.begin();
        while (it != nets.end()) {
            std::set<std::string> new_node_set;
            for (const auto &node: it->node_set) {
                new_node_set.insert(mapping[node]);
            }
            if (new_node_set.size() == 1) {
                it = nets.erase(it);
            } else {
                it->node_set = new_node_set; 
                it++;
            }
        }
        // update nodes' net list
        for (unsigned i = 0; i < nets.size(); i++) {
            Net &net = nets[i];
            for (auto node_name: net.node_set) {
                Node &node = nodes[node_name];
                node.net_set.insert(i);
                if (node.node_set.count(net.driver) > 0) {
                    node.drives.insert(i);
                    net.driver = node_name;
                }
            }
        }
    }

    bool is_connected(std::string n1, std::string n2) {
        std::set<std::string> visited;
        std::stack<std::string> st;
        st.push(n1);
        while (!st.empty()) {
            std::string n = st.top(); st.pop();
            visited.insert(n);
            for (const auto nn: neighbors(n)) {
                if (nn == n2)
                    return true;
                if (visited.count(nn) == 0)
                    st.push(nn);
            }
        }
        return false;
    }

    std::set<std::string> neighbors(std::string n) {
        std::set<std::string> result;
        for (auto net_num: nodes[n].net_set) {
            for (const auto &nn: nets[net_num].node_set) {
                result.insert(nn);
            }
        }
        return result;
    }
};

#endif