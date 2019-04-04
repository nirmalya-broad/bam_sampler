#ifndef _BAM_READER_HPP
#define _BAM_READER_HPP

#include <cstring>
#include <cmath>
#include <regex>
#include <htslib/sam.h>
#include <regex.h>

class bam_reader {
    public:

    bam_reader() = delete;

    bam_reader(std::string& infile_str) {
        const char* format = NULL;
        if (has_suffix(infile_str, "sam")) {
            format = "r";
        } else if (has_suffix(infile_str, "bam")) {
            format = "rb";
        } else {
            std::string lstr = "File with illegal suffix: " + infile_str + "\n";
            throw std::runtime_error(lstr);
        }

        const char* infile_cstr = infile_str.c_str();
        if (!(lfp = sam_open(infile_cstr, format))) {
            throw std::runtime_error("Error in sam_open");
        } else {
            printf("sam_open: %s\n", infile_cstr);
        }

        lhdr = sam_hdr_read(lfp);

    }

    bam_hdr_t* get_sam_header() {
        return lhdr;
    }
    
    bam_hdr_t* get_sam_header_dup() {
        return bam_hdr_dup(lhdr);
    }

    bool has_suffix(const std::string &str, const std::string &suf)
    {
        return str.size() >= suf.size() &&
           str.compare(str.size() - suf.size(), suf.size(), suf) == 0;
    }


    ~bam_reader() {
        bam_hdr_destroy(lhdr);   
        sam_close(lfp); // clean up 
    }

    int sam_read(bam1_t* lread) {
        int ret_val = sam_read1(lfp, lhdr, lread);
        return ret_val;
    }

    private:
    std::string infile_str;    
    htsFile *lfp = NULL;
    bam_hdr_t *lhdr = NULL;


};
#endif


