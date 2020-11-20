#ifndef __FPGA_HPP__
#define __FPGA_HPP__

#include <vector>
#include <cinttypes>

struct FPGA {
    uint64_t resources[10];
    std::vector<uint64_t> int_list;

    FPGA(int int_list_size): int_list(int_list_size, 0) {
        for (int i = 0; i < 10; i++)
            resources[i] = 0;
    }

    FPGA(uint64_t *_resources, const std::vector<uint64_t> &_int_list): int_list(_int_list) {
        for (int i = 0; i < 10; i++) {
            resources[i] = _resources[i];
        }
    }
};

#endif
