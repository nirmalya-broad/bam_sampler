#ifndef __BAM_STAT_HPP
#define __BAM_STAT_HPP

#include <htslib/sam.h>


class bam_stat {

    public:

    bam_stat() = delete;
   
    bam_stat(std::string& infile_str) {
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
        }

        lhdr = sam_hdr_read(lfp);
        
        read1 = bam_init1();
        read2 = bam_init1();
    }

    bool has_suffix(const std::string &str, const std::string &suf) {
        return str.size() >= suf.size() &&
           str.compare(str.size() - suf.size(), suf.size(), suf) == 0;
    }

    unsigned long get_frag_count() {
        if (frag_count == -1) {
            frag_count = exe_frag_count();
        }
        return (frag_count);
    }

    unsigned long exe_frag_count() {

        unsigned long frag_count = 0;
        while(sam_read1(lfp, lhdr, read1) >= 0) {
            // Increment for read1

            if (sam_read1(lfp, lhdr, read2) < 0) {

                std::string lstr = "Failed to obtain second read from "
                    " get_read_count.\n";
                throw std::runtime_error(lstr);
            } else {
                // Increment for read2
            }

            frag_count++;
        }

        return frag_count;
    }

    ~bam_stat() {
        bam_hdr_destroy(lhdr);
        sam_close(lfp); // clean up
        bam_destroy1(read1);
        bam_destroy1(read2);
    }

    private:
        std::string infile_str;
        htsFile *lfp = NULL;
        bam_hdr_t *lhdr = NULL;
        long frag_count = -1;
        bam1_t* read1 = NULL;
        bam1_t* read2 = NULL;

};

#endif
