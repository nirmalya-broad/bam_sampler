#include <iostream>
#include <string>
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

    private:
        po::options_description desc;
        std::string infile_str;
        std::string outfile_str;
        long top_seed;
        double sample_p;
        htsFile *fp_in = NULL;
        htsFile *fp_out = NULL;
        bam_hdr_t *hdr_in = NULL;
        bam_hdr_t *hdr_out = NULL;

};

void sample_bamc::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: sample_bamc --infile <sam/bam> --outfile <sam/bam>"
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

void sample_bamc::initialize() {

    bam_handle_init(infile_str, NULL, fp_in, hdr_in, "r");
    bam_handle_init(outfile_str, hdr_in, fp_out, hdr_out, "w");

}

void sample_bamc::main_func() {

    // Open the input sam/bam file; the assumption is that the file is 
    // sorted with respect to query name. So pair of reads would be one after
    // another. 

    // 1. Open the input file in the read mode

    // 2. Open the output file in the write mode
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


