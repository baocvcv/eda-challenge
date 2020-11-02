#include "common.h"
#include "cxxopts.hpp"
#include "solver.hpp"

using namespace std;

int main(int argc, const char *argv[]) {
    Config config;
    string input_dir;

    try {
        cxxopts::Options options("Graph split", "Description");
        options.add_options()
            ("i,input", "Input file directory", cxxopts::value<string>())
            ("fix-mincut", "Default mode", cxxopts::value<bool>()->default_value("true"))
            ("int-mincut", "Use int-mincut", cxxopts::value<bool>()->default_value("false"))
            ("ffd-mincut", "Use ffd-mincut", cxxopts::value<bool>()->default_value("false"))
            ("clk-mincut", "Use clk-mincut", cxxopts::value<bool>()->default_value("false"))
            ("num", "Max clk num", cxxopts::value<int>())
            ("clks", "List of clk names", cxxopts::value<vector<string>>())
            ("mean-mincut", "Use mean-mincut", cxxopts::value<bool>()->default_value("false"))
            ("percentage", "Mean deviation percentage", cxxopts::value<int>()->default_value("20"))
            ("h,help", "Print help")
            ;

        options.parse_positional({"input"});

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            cout << options.help({""}) << endl;
            exit(0);
        }

        if (result.count("input")) {
            input_dir = result["input"].as<string>();
        } else {
            cerr << "Missing input file directory" << endl;
            cout << options.help({""}) << endl;
            exit(0);
        }

        if (result.count("int-mincut"))
            config.add_mode(INT_MINCUT);
        if (result.count("ffd-mincut"))
            config.add_mode(FFD_MINCUT);
        if (result.count("clk-mincut")) {
            config.add_mode(CLK_MINCUT);
            if (result.count("num")) 
                config.clk_num = result["num"].as<int>();
            if (result.count("clks"))
                config.clk_list = result["clks"].as<vector<int>>();
        }
        if (result.count("mean-mincut")) {
            config.add_mode(MEAN_MINCUT);
            if (result.count("percentage"))
                config.mean_percent = result["percentage"].as<int>();
        }

    } catch (const cxxopts::OptionException& e) {
        cerr << "error parsing options: " << e.what() << endl;
        exit(1);
    }

    Solver solver(input_dir, config);

    return 0;
}