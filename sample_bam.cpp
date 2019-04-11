#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <random>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

class sample_bamc {

    public:
        void print_help();
        bool parse_args(int argc, char* argv[]); 
        void initialize();
        void main_func();
        void free_vars();
        void bam_handle_init(std::string fname_str, const bam_hdr_t* lhdr,
            htsFile*& lfp_out, bam_hdr_t*& lhdr_out, std::string open_mode);
        std::string get_format(std::string lfname_str, std::string open_mode);
        htsFile* get_hts_handle(std::string lfname_str, std::string format);
        bool has_suffix(const std::string &str, const std::string &suffix);
        unsigned long get_frag_count(std::string& infile_str);

    private:
        po::options_description desc;
        std::string infile_str;
        std::string outfile_str;
        unsigned long long top_seed;
        unsigned long interval;
        double sample_p;
        unsigned long infile_frag_count;
        htsFile *fp_in = NULL;
        htsFile *fp_out = NULL;
        bam_hdr_t *hdr_in = NULL;
        bam_hdr_t *hdr_out = NULL;
        long read_limit;

};

void sample_bamc::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: sample_bam --infile <sam/bam> --outfile <sam/bam>"
        " --top_seed <top_seed_number>  --sample_p <sample percentage>"
        "\n\n";
}

bool sample_bamc::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input sam/bam file.")
        ("outfile,o", po::value<std::string>(&outfile_str), "output file.")
        ("top_seed,t", po::value(&top_seed),
            "Top seed for sampling.")
        ("sample_p,s", po::value(&sample_p),
            "Percentage of reads to return.")
        ("interval", po::value(&interval)->default_value(1000000),
            "Interval for shuffle.")
        ("read_limit,r", po::value(&read_limit)->default_value(-1),
            "Pair of reads to process")
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

 
    if (vm.count("top_seed")) {
        std::cout << "top_seed is set to: " << top_seed << "\n";
    } else {
        all_set = false;
        std::cout << "Error: top_seed is not set.\n";
    }

    if (vm.count("sample_p")) {
        std::cout << "sample_p is set to: " << sample_p << "\n";
    } else {
        std::cout << "Error: sample_p is not set.\n";
    }

    std::cout << "Interval is set to: " << interval << "\n";

    return all_set;

}

bool sample_bamc::has_suffix(const std::string &str, const std::string &suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}


// Note that this function does a dulication; so there mighe be memory
// issue that need to be taken care of.
// Also, it is passed a lhdr as input, if that is null the function open
// the input file and assign the header to it. Otherwise, the function 
// copies the header to the output header.

std::string sample_bamc::get_format(std::string lfname_str, std::string open_mode) {

    std::string format;
    if (has_suffix(lfname_str, "sam")) {
        format = open_mode;
    } else if (has_suffix(lfname_str, "bam")) {
        format = open_mode + "b";
    } else {
        std::string lstr = "File with illegal suffix: " + lfname_str + "\n";
        throw std::runtime_error(lstr);
    }
    return format;
}


htsFile* sample_bamc::get_hts_handle(std::string lfname_str, std::string format) {
    htsFile* lfp = NULL;
    const char* fname_cstr = lfname_str.c_str();
    const char* format_cstr = format.c_str();
    printf("%s\t%s\n", fname_cstr, format_cstr);
    std::cout << "format from get_hts_handle: " << format << "\n";
    if (!(lfp = sam_open(fname_cstr, format_cstr))) {
        std::cout << "Error in opening the file: " << lfname_str << "\n";
    } else {
        std::cout << "Successfully opened the file: " << lfname_str << "\n";
    }
    return lfp;

}

// If the input header is NOT null, it duplicates the input header to the 
// output header, otherwise it directly copies the header from the
// input file.

void sample_bamc::bam_handle_init(std::string fname_str, const bam_hdr_t* lhdr,
        htsFile*& lfp_out, bam_hdr_t*& lhdr_out, std::string open_mode) {

    std::string format = get_format(fname_str, open_mode);
    htsFile* lfp = get_hts_handle(fname_str, format);
    lfp_out = lfp;

    if (open_mode == "r") {
        // Create one and return
        lhdr_out = sam_hdr_read(lfp);
    } else {
        // Create a duplicate and return
        lhdr_out = bam_hdr_dup(lhdr);
    }


    if (open_mode == "w") {
        if (sam_hdr_write(lfp_out, lhdr_out) < 0 ) {
            std::cout << "Error in sam_hdr_write." << "\n";
        } else {
            std::cout << "Wrote the header successfully." << "\n";
        }

    }
}


// Returns the number of reads in a file
unsigned long sample_bamc::get_frag_count(std::string& infile_str) {

    htsFile *lfp = NULL;
    bam_hdr_t *lhdr = NULL;

    bam_handle_init(infile_str, NULL, lfp, lhdr, "r");

    bam1_t* read1 = bam_init1();
    bam1_t* read2 = bam_init1();
    unsigned long frag_count = 0;
    while(sam_read1(lfp, lhdr, read1) >= 0) {
        // Increment for read1

        if (sam_read1(lfp, lhdr, read2) < 0) {
            bam_hdr_destroy(lhdr);
            sam_close(lfp);
            bam_destroy1(read1);
            bam_destroy1(read2);

            std::string lstr = "Failed to obtain second read from get_read_count.\n";
            throw std::runtime_error(lstr);
        } else {
            // Increment for read2
        }

        frag_count++;
    }

    bam_hdr_destroy(lhdr);
    sam_close(lfp);
    bam_destroy1(read1);
    bam_destroy1(read2);

    return frag_count;

}

void sample_bamc::initialize() {

    bam_handle_init(infile_str, NULL, fp_in, hdr_in, "r");
    bam_handle_init(outfile_str, hdr_in, fp_out, hdr_out, "w");

    // Initialize the number of reads in the input file.
    infile_frag_count = get_frag_count(infile_str);
    std::cout << "Infile_frag_count: " << infile_frag_count << "\n";

}

void sample_bamc::main_func() {

    // Open the input sam/bam file; the assumption is that the file is 
    // sorted with respect to query name. So pair of reads would be one after
    // another. 

    std::vector<unsigned long> seq_vec;
    std::set<unsigned long> write_set;
    std::mt19937_64 r_engine(top_seed);

    bam1_t* read1 = bam_init1();
    bam1_t* read2 = bam_init1();

    // Total number of reads loaded/read till this point.
    unsigned long total_frag_ind = 0;
   
    unsigned long frag_left = infile_frag_count; 
    unsigned long write_ind = 0;
    unsigned long frag_this_iter = 0;
    unsigned long write_lim = 0;
    unsigned long total_frag_written = 0;

    while(sam_read1(fp_in, hdr_in, read1) >= 0) {
        // Increment for read1

        if (sam_read1(fp_in, hdr_in, read2) < 0) {
            std::string lstr = "Failed to read second read.\n";
            throw std::runtime_error(lstr);
        } else {
            // Increment for read2
        }

        if (total_frag_ind % interval == 0) {


            // Do the assignment
            if (frag_left >= interval) {
                frag_this_iter = interval;
                frag_left = frag_left - interval;
            } else {
                frag_this_iter = frag_left;
                frag_left = 0;
            }
            seq_vec.assign(frag_this_iter, 0);
            write_set.clear();
            for (unsigned long j = 0; j < frag_this_iter; j++) {
                unsigned long l_seq_val = total_frag_ind + j;
                seq_vec[j] = l_seq_val;       
            }

            // Do the shuffle
            shuffle(seq_vec.begin(), seq_vec.end(), r_engine);
            write_ind = 0;     
            double write_lim_f = (frag_this_iter * sample_p) /100.0;
            write_lim = (unsigned long) write_lim_f;  
            for (unsigned int j = 0; j < write_lim ; j++) {
                write_set.insert(seq_vec[j]);
            }

            
        }

        // Write if appropriate
        if (write_ind < write_lim) {
            if (write_set.find(total_frag_ind) != write_set.end()) {

                
                if (sam_write1(fp_out, hdr_out, read1) < 0) {
                    std::cout << "Problem with sam_write1 for read1" << "\n";
                }

                if (sam_write1(fp_out, hdr_out, read2) < 0) {
                    std::cout << "Problem with sam_write1 for read2" << "\n";
                }

                write_ind++;
                total_frag_written++;
            }
        }

        total_frag_ind++; 
        if (total_frag_ind == read_limit) {
            break;
        }
    }

    std::cout << "Total frag written: " << total_frag_written << "\n";

    bam_destroy1(read1);
    bam_destroy1(read2);
}

void sample_bamc::free_vars() {

    bam_hdr_destroy(hdr_in);
    bam_hdr_destroy(hdr_out);
    sam_close(fp_in);
    sam_close(fp_out);

}

int main(int argc, char** argv) {
    sample_bamc sbc;
    bool all_set = true;

    try {
        all_set = sbc.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    } catch(...) {
        return 0;
    } 

    if (!all_set) {
        sbc.print_help();
        return 0;
    }

    try {
        sbc.initialize();
        sbc.main_func();
        sbc.free_vars();
    } catch(const std::runtime_error& e) {
        std::cerr << "error: "  << e.what() << "\n";
        return 1;
    }

    return 0;
}


