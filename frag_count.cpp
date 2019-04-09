#include <iostream>
#include "bam_stat.hpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;
po::options_description desc;

void print_help() {
        std::cout << desc << "\n";
    std::cout << "Usage: seed_gen --infile <sam/bam>"
        "\n\n";
}

int main(int argc, char** argv) {

    std::string infile_str;
    
    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input sam/bam file.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        print_help();
        return 0;
    } else {
    }

    if (vm.count("infile")) {
        std::cout << "Infile is set to: " << infile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: infile is not set.\n";
    }

    if (!all_set) {
        print_help();
        return 0;
    }

    bam_stat in_stat(infile_str);
    unsigned long frag_count = in_stat.get_frag_count();
    std::cout << "@frag_count " << frag_count << "\n";

}
