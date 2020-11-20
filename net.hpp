#ifndef __NET_HPP__
#define __NET_HPP__

#include <set>
#include <cinttypes>

struct Node;

struct Net {
    int index;

    std::set<std::string> node_set;

    std::string driver;

    uint32_t cost;

    Net (int _index, const std::string &_driver, int _cost): index(_index), driver(_driver), cost(_cost) {
        node_set.insert(_driver);
    }
};

#endif
