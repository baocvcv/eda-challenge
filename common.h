#ifndef __COMMON_H__
#define __COMMON_H__

#include <bits/stdc++.h>

const std::string WHITE_SPACES = " \n\t\r\f\v";

enum Mode {
    FIX_MINCUT = 0x1,
    INT_MINCUT = 0x2,
    FFD_MINCUT = 0x4,
    CLK_MINCUT = 0x8,
    MEAN_MINCUT = 0x10
};

int operator&(Mode a, Mode b) {
    return static_cast<int>(a) & static_cast<int>(b);
}

int operator|(Mode a, Mode b) {
    return static_cast<int>(a) | static_cast<int>(b);
}

struct Config {
    Mode mode;
    int clk_num;
    std::vector<int> clk_list;
    int mean_percent;

    // fixed assignment info
    std::map<std::string, int> fixed_assignment;
    std::set<int> fixed_fpgas;
    std::set<int> free_fpgas;

    Config(): mode(FIX_MINCUT), clk_num(-1), mean_percent(20) {}

    void add_mode(Mode m) {
        mode = static_cast<Mode>(static_cast<int>(mode) | static_cast<int>(m));
    }
};

const std::string node_def_file = "design.are";
const std::string net_def_file = "design.net";
const std::string fpga_res_file = "design.info";
const std::string constraint_file = "design.fix";
const std::string result_file = "design.output";
const std::string eval_file = "design.rpt";

typedef std::pair<bool, std::string> Error;

std::string resource_names[] = {"PIO", "INT", "FF", "LUT", "BUFG", "TBUF", "DCM", "BRAM", "DSP", "PPC"};

#endif