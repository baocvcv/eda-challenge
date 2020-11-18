#ifndef __NODE_HPP__
#define __NODE_HPP__

#include "common.h"
#include "net.hpp"

#include <set>
#include <string>

/**
 * @brief (hyper) graph node
 * 
 */
struct Node {
    // name of the node
    std::string node_name;

    // whether node is a leaf node
    bool is_leaf;

    // nodes contained in this hyper node
    std::set<std::string> node_set;

    // the set of net connected to this node
    std::set<int> net_set;

    // resources used by this node
    int resources[10];

    /* attributes */
    std::set<int> drives;
    std::set<std::string> clk_names;
    bool is_clk;
    bool is_ff;

    /* fixed assignment */
    int assigned_fpga;

    Node():
        is_leaf(true),
        is_clk(false),
        is_ff(false),
        assigned_fpga(-1)
    {
        for (int i = 0; i < 10; i++)
            resources[i] = 0;
    }

    Node(std::string _node_name, int *_resources, const std::set<std::string> &_clk_names, bool _is_ff):
        is_leaf(true),
        node_name(_node_name),
        clk_names(_clk_names),
        is_clk(_clk_names.size() > 0),
        is_ff(_is_ff),
        assigned_fpga(-1)
    {
        if (_resources) {
            for (int i = 0; i < 10; i++)
                resources[i] = _resources[i];
        }
    }

    void set_fixed(int assignment) { assigned_fpga = assignment; }

    bool is_fixed() const { return assigned_fpga >= 0; }

    void merge_with(const Node &rhs) {
        // net related merge is done in Graph.update_net()
        node_set.insert(rhs.node_set.begin(), rhs.node_set.end());
        net_set.clear();
        for (int i = 0; i < 10; i++) {
            resources[i] += rhs.resources[i];
        }
        drives.clear();
        clk_names.insert(rhs.clk_names.begin(), rhs.clk_names.end());
        is_clk |= rhs.is_clk;
        is_ff |= rhs.is_clk;
        if (assigned_fpga < 0 || rhs.assigned_fpga < 0) {
            // at least one is not fixed
            assigned_fpga = std::max(assigned_fpga, rhs.assigned_fpga);
        } else if (assigned_fpga >= 0 && rhs.assigned_fpga >= 0 && assigned_fpga != rhs.assigned_fpga) {
            // both are fixed and not of the same type
            std::cerr << "Failed merging two nodes: " << node_name << " " << rhs.node_name;
            exit(1);
        }
    }
};


#endif