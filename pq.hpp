#ifndef __PR_HPP__
#define __PR_HPP__

#include <vector>
#include <map>

template <class T>
struct my_pq {
    std::vector<T> vec;
    std::map<T, int> cost;

    int size() { return vec.size(); }

    bool contains(const T &x) { return cost.count(x) > 0; }

    void insert(const T &key, int c) {
        vec.push_back(key);
        cost[key] = c;
        percolate_up();
    }

    T top() {
        if (vec.size() == 0)
            return T();
        return vec[0];
    }

    T pop() {
        if (vec.size() == 0)
            return T();
        T tmp = vec[0];
        vec[0] = vec.back();
        vec.pop_back();
        cost.erase(tmp);
        percolate_down();
        return tmp;
    }

    void update(const T &key, int c) {
        if (cost[key] == c)
            return;

        cost[key] = c;
        int i;
        for (i = 0; i < vec.size(); i++)
            if (vec[i] == key)
                break;
        _swap(i, vec.size()-1);
        percolate_down(i);
        percolate_up();
    }

    void percolate_up(int idx = -1) {
        int c = (idx == -1) ? (vec.size()-1) : idx;
        int p = parent(c);
        while (p >= 0 && cost[vec[p]] < cost[vec[c]]) {
            _swap(p, c);
            c = p;
            p = parent(c);
        }
    }

    void percolate_down(int idx = 0) {
        int root = idx;
        int l = l_child(root);
        int r = r_child(root);
        int size = vec.size();
        while (l < size) {
            if (r < size) { // has right child
                int max_n = max(root, l, r);
                if (max_n == root) {
                    break;
                } else {
                    _swap(max_n, root);
                    root = max_n;
                }
            } else { // no right child
                if (cost[vec[root]] >= cost[vec[l]]) {
                    break;
                } else {
                    _swap(root, l);
                    root = l;
                }
            }
            l = l_child(root);
            r = r_child(root);
        }
    }

    void _swap(int a, int b) {
        T tmp = vec[a];
        vec[a] = vec[b];
        vec[b] = tmp;
    }

    int max(int a, int b, int c) {
        T ta = vec[a];
        T tb = vec[b];
        T tc = vec[c];
        return (cost[ta] >= cost[tb]) ? (cost[ta] >= cost[tc] ? a : c) : (cost[tb] >= cost[tc] ? b : c);
    }

    int l_child(int x) { return 2 * x + 1; }
    int r_child(int x) { return 2 * x + 2; }
    int parent(int x) { return (x - 1) / 2; }
};

#endif