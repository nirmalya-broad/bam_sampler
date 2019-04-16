#include <iostream>
#include <fstream>
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
        double start_read_p;
        double end_read_p;
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
        double start_read_p;
        double end_read_p;
        std::vector<double> sample_p_vec;
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
    start_read_p = args_o.start_read_p;
    end_read_p = args_o.end_read_p;
    
}


void args_c::print_help() {
        std::cout << desc << "\n";
    std::cout << "Usage: seed_gen --infile <sam/bam> --outfile <seed_table>"
        " --run_type <steps/depth> [ --main_seed <main_seed_number>"
        " --repeat_num <number of repeats> --steps_num <number of steps>"
        " --depth_p <depth percentage> --start_read_p --end_read_p]"
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
        ("start_read_p", po::value(&start_read_p), "Percentage of reads for start breakpoint")
        ("end_read_p", po::value(&end_read_p)->default_value(100.0), "Percentage of reads for end breakpoint")
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
        } else if (0 == run_type.compare("depth")) {
            if (depth_p <= 0) {
                all_set = false;
                std::cout << "Error: depth need to be a positive number.\n";
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
    std::ofstream outf(outfile_str);

    // Let's get the number of fragments in the file.
    std::cout << "infile_frag_count: " << infile_frag_count << "\n";

    
    if (0 == run_type.compare("steps")) {
        
        // get the percentage of reads to sample 
        // If steps equal to 1; then it would indicate 100% of the reads
        // A typical step size would be more than 10 and may be less than 20.

        double end_p = end_read_p;
        double start_log = log2(start_read_p);
        double end_log = log2(end_p);
        double diff_start_end = end_log - start_log;
        double step_val = diff_start_end / (steps_num - 1);


        for (int j = 0; j < steps_num; j++) {
            double l_sample_p_log = start_log + j * step_val;
            double l_sample_p = exp2(l_sample_p_log);
            sample_p_vec.push_back(l_sample_p);
        }


        for (int j = 0 ; j < steps_num; j++) {
            for (int k = 0; k < repeat_num; k++) {
                int l_step = j + 1;
                int l_repeat = k + 1;
                if (j == (steps_num -1) && end_read_p == 100.0)  {
                    outf << "S" << l_step << "\t" << "R1" << "\t" <<  sample_p_vec[j] << "\t" << r_engine() << "\n";
                    break;
                } else {
                    outf << "S" << l_step << "\t" << "R" << l_repeat 
                        << "\t" <<  sample_p_vec[j] << "\t" << r_engine() << "\n";
                }
            }
        }

    } else if (0 == run_type.compare("depth")) {
        // Get the depth in terms of percentage of reads
        for (int k = 0; k < repeat_num; k++) {
            int l_repeat = k + 1;
            int l_step = 1;
            outf << "S" << l_step << "\t" << "R" << l_repeat
                << "\t" << depth_p << "\t" << r_engine() << "\n";
        }
    }
    
    // Now create the output file that would be used for generating the 
    // sampling distribution of the bam file.    
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
