#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <random>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <boost/program_options.hpp>

#include "bam_reader.hpp"
#include "bam_stat.hpp"

namespace po = boost::program_options;

// This class would generate a table with list of top seeds along with the
// percentage of reads with replicates. The basic idea is that it would 
// produce the top seeds from a main seed in a pseudorandom (deterministic)
// fashion.


class args_c {
    public:
        po::options_description desc;
        std::string infile_str;
        std::string outfile_str;
        std::string run_type;
        long main_seed;
        int steps_num;
        double depth_p;
        unsigned int repeat_num;
        void print_help();
        bool parse_args(int argc, char* argv[]);

};

class seed_genc {

    public:
        seed_genc(args_c args_o);
        void initialize();
        void main_func();
        void free_vars();

    private:
        std::string infile_str;
        std::string outfile_str;
        std::string run_type;
        long main_seed;
        int steps_num;
        double depth_p;
        unsigned int repeat_num;
        unsigned long infile_frag_count;
        bam_reader in_reader;        
        bam_stat in_stat;
        long depth;
};

seed_genc::seed_genc(args_c args_o)
    : in_reader(args_o.infile_str),
    in_stat(args_o.infile_str) {
    infile_str = args_o.infile_str;
    outfile_str = args_o.outfile_str;
    run_type = args_o.run_type;
    main_seed = args_o.main_seed;
    steps_num = args_o.steps_num;
    depth_p = args_o.depth_p;
    repeat_num = args_o.repeat_num;
    
}


void args_c::print_help() {
        std::cout << desc << "\n";
    std::cout << "Usage: seed_gen --infile <sam/bam> --outfile <seed_table>"
        " --run_type <steps/depth> [ --main_seed <main_seed_number>"
        " --repeat_num <number of repeats> --steps_num <number of steps>"
        " --depth_p <depth percentage>]"
        "\n\n";
}

bool args_c::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input sam/bam file.")
        ("outfile,o", po::value<std::string>(&outfile_str), "output file.")
        ("main_seed,m", po::value(&main_seed)->default_value(12345),
            "Man seed for generating the top seeds.")
        ("run_type,r", po::value(&run_type), "Type of run: steps/depth.")
        ("repeat_num", po::value(&repeat_num)->default_value(1), "Number of repeats.")
        ("steps_num,s", po::value(&steps_num)->default_value(-1),
            "Number of steps.")
        ("depth_p,d", po::value(&depth_p)->default_value(-1),
            "Percentage of depth")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        return 0;
    } else {
    }

    if (vm.count("infile")) {
        std::cout << "Infile is set to: " << infile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: infile is not set.\n";
    }

    if (vm.count("outfile")) {
        std::cout << "Outfile is set to " << outfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: outfile is not set.\n";
    }

   if (vm.count("main_seed")) {
        std::cout << "main_seed is set to: " << main_seed << "\n";
    } else {
        all_set = false;
        std::cout << "Error: main_seed is not set.\n";
    }

    if (vm.count("run_type")) {
        std::cout << "run_type is set to: " << run_type << "\n";
        if (0 == run_type.compare("steps")) {
            if (steps_num <= 0) {
                all_set = false;
                std::cout << "Error: steps_num need to be a positive integer.\n";
            } else {    
                std::cout << "steps_num is set to: " << steps_num << "\n";
            }
        } else if (0 == run_type.compare("depth_p")) {
            if (depth_p <= 0) {
                all_set = false;
                std::cout << "Error: depth need to be a positive integer.\n";
            } else {
                std::cout << "depth is set to: " << depth_p << "\n";
            }
        } else {
            std::string lstr = "Error: Illegal run_type: " + run_type;
        }

    } else {
        all_set = false;
        std::cout << "Error: run_type is not set.\n";
    }


    if (repeat_num <= 0) {
        all_set = false;
        std::cout << "Error: repeat_num need to be a positive number.\n";
    } else {
        std::cout << "repeat_num is set to: " << repeat_num << "\n";
    }

    return all_set;
}

void seed_genc::initialize() {
    infile_frag_count = in_stat.get_frag_count();
    
}

void seed_genc::main_func() {
    std::mt19937_64 r_engine(main_seed);

    bam1_t* read1 = bam_init1();
    bam1_t* read2 = bam_init1();


    unsigned long total_frag = 0;
      while(in_reader.sam_read(read1) >= 0) {
        // Increment for read1

        if (in_reader.sam_read(read2) < 0) {
            std::string lstr = "Failed to read second read.\n";
            throw std::runtime_error(lstr);
        } else {
            // Increment for read2
        }
        total_frag++;
    }

    std::cout << "total_frag: " << total_frag << "\n";
     
    bam_destroy1(read1);
    bam_destroy1(read2);
}


void seed_genc::free_vars() {

}

int main(int argc, char** argv) {



    args_c args_o;
    bool all_set = true;

    try {
        all_set = args_o.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    } catch(...) {
        return 0;
    }

    if (!all_set) {
        args_o.print_help();
        return 0;
    }

    seed_genc sgc(args_o);
    try {
        sgc.initialize();
        sgc.main_func();
        sgc.free_vars();
    } catch(const std::runtime_error& e) {
        std::cerr << "error: "  << e.what() << "\n";
        return 1;
    }

    return 0;
}
