#ifndef __NET_HPP__
#define __NET_HPP__

#include <set>

struct Node;

struct Net {
    std::set<Node*> node_set;

    Node* driver;

    int cost;

    Net (Node *_driver, int _cost): driver(_driver), cost(_cost) {
        node_set.insert(_driver);
    }
};

#endif
