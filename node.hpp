#ifndef __NODE_HPP__
#define __NODE_HPP__

#include "net.hpp"

#include <set>
#include <string>

/**
 * @brief (hyper) graph node
 * 
 */
struct Node
{
    bool is_leaf;
    // leaf node
    std::string node_name;

    // nodes contained in this hyper node
    std::set<Node*> node_set;

    // net set
    std::set<Net*> net_set;

    // total resources
    int resources[10];

    /* attributes */
    std::set<Net*> drives;
    std::vector<std::string> clk_names;
    bool is_clk;
    bool is_ff;

    /* assignments */
    bool is_fixed;
    std::set<int> assigned_fpga;

    Node(std::string _node_name, int *_resources, const std::vector<std::string> &_clk_names, bool _is_ff):
        is_leaf(true),
        node_name(_node_name),
        clk_names(_clk_names),
        is_clk(_clk_names.size() > 0),
        is_ff(_is_ff),
        is_fixed(false),
    {
        for (int i = 0; i < 10; i++)
            resources[i] = _resources[i];
    }

    Node(const std::vector<Node*> &nodes):
        is_leaf(false)
        // node_name(_node_name),
        // clk_names(_clk_names),
        // is_clk(_clk_names.size() > 0),
        // is_ff(_is_ff),
        // is_fixed(false),
        // assigned_fpga(-1)
    {
        // for (int i = 0; i < 10; i++)
        //     resources[i] = _resources[i];

        // merge the nodes
        // modify the nets involved
    }

    void set_fixed(int assignment) {
        assigned_fpga.insert(assignment);
        is_fixed = true;
    }
};


#endif