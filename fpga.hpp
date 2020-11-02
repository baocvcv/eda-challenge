#ifndef __FPGA_HPP__
#define __FPGA_HPP__

#include <vector>

struct FPGA {
    int resources[10];
    std::vector<int> int_list;

    FPGA(int *_resources, const std::vector<int> &_int_list): int_list(_int_list) {
        for (int i = 0; i < 10; i++) {
            resources[i] = _resources[i];
        }
    }
};

#endif