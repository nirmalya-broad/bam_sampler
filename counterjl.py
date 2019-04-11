import re
import os
from subprocess import call
from alignerutils import bam_to_sam
from alignerutils import sort_by_qname

class CounterJL:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd

    def get_samfile_paired(self, sorted_bam):
        """ This function checks the existence of sam file, and generates
        one from the bam file. Since it is single ended data, we don't care 
        about the order of the reads.        
        """
        cldict = self.cldict
        samtools_path = cldict.samtools
        if not os.path.isfile(sorted_bam):
            lstr = "Bam file not found: " + sorted_bam
            raise StandardError(lstr)
        else:
            print("sorted_bam: " + sorted_bam)
            unsorted_bam = re.sub(r'_pe.bam$', '_u.bam', sorted_bam)
            print("unsorted_bam: " + unsorted_bam)
            outsamfile = re.sub(r'_u.bam$', '.sam', unsorted_bam)
            print("outsamfile: " + outsamfile)
            if os.path.isfile(outsamfile):
                print("Returning " + outsamfile)
                return outsamfile
            else:
                if not os.path.isfile(unsorted_bam):
                    sort_by_qname(samtools_path, sorted_bam, unsorted_bam)
                if not os.path.isfile(outsamfile):
                    bam_to_sam(samtools_path, unsorted_bam, outsamfile)
                return outsamfile

    def get_samfile_single(self, sorted_bam):
        """ This function checks the existence of sam file, and generates
        one from the bam file. Since it is single ended data, we don't care 
        about the order of the reads.        
        """
        cldict = self.cldict
        if not os.path.isfile(sorted_bam):
            lstr = "Bam file not found: " + sorted_bam
            raise StandardError(lstr)
        else:
            outsamfile = re.sub("_se.bam$", '.sam', sorted_bam)
            if os.path.isfile(outsamfile):
                return outsamfile
            else:
                samtools_path = cldict.samtools
                bam_to_sam(samtools_path, sorted_bam, outsamfile)
                return outsamfile
       
    def clean_sam_single(self, outsamfile):
        print("Removing: " + outsamfile)
        os.remove(outsamfile)

    def count_paired(self, sample_id, ref_acc, sorted_bam, outdir):
        cldict = self.cldict
        JLCounter = cldict.JLCounter
        ldelim = cldict.ldelim
        Data_dir = cldict.Data_dir
        Script_dir = cldict.basepath
        Patho_dir = cldict.Patho_dir
        sampd = self.sampd
        project_id = sampd.project_id
        count_strand_rev = cldict.count_strand_rev
        #paired_only = cldict.paired_only
        patho_temp_bamdir = outdir + ldelim + "temp_bamdir"
        patho_gff = Data_dir + ldelim + ref_acc + "_ALL.gff"
        shell_script_dir = Script_dir + ldelim + "shell_scripts"
        countfile_str = outdir + ldelim + sample_id + "_" + ref_acc + ".counts"
        outsamfile = self.get_samfile_paired(sorted_bam)
        count_cmd = "sh " + JLCounter + " " + outsamfile + " " + patho_gff + \
            " " + shell_script_dir + " " + patho_temp_bamdir + " " + \
            countfile_str
        if count_strand_rev == "Y":
            count_cmd += " -STRAND_REV Y"
        #if paired_only:
        #    count_cmd += " -paired_only"
        print("Starting JL counter for RtS: " + count_cmd)
        call(count_cmd.split())
        self.clean_sam_single(outsamfile)

    def count_single(self, sample_id, ref_acc, sorted_bam, outdir, strand_rev):
        cldict = self.cldict
        JLCounter = cldict.JLCounter
        ldelim = cldict.ldelim
        Data_dir = cldict.Data_dir
        Script_dir = cldict.basepath
        Patho_dir = cldict.Patho_dir
        sampd = self.sampd
        project_id = sampd.project_id
        patho_temp_bamdir = outdir + ldelim + "temp_bamdir"
        patho_gff = Data_dir + ldelim + ref_acc + "_ALL.gff"
        shell_script_dir = Script_dir + ldelim + "shell_scripts"
        countfile_str = outdir + ldelim + sample_id + "_" + ref_acc + ".counts"
        outsamfile = self.get_samfile_single(sorted_bam)
        
        l_strand_rev = strand_rev
        count_cmd = "sh " + JLCounter + " " + outsamfile + " " + patho_gff + \
            " " + shell_script_dir + " " + patho_temp_bamdir + " " + \
            countfile_str + " -STRAND_REV " + l_strand_rev
        LC_method_val = cldict.LC_method_val
        lc_lower = LC_method_val.lower()
        if lc_lower == 'allseq':
            print("Starting JL counter for AllSeq: " + count_cmd)
        elif lc_lower == 'rts' or lc_lower == 'rts-ts':
            print("Starting JL counter for RtS: " + count_cmd)
        print("count_cmd: " + count_cmd)
        call(count_cmd.split())
        self.clean_sam_single(outsamfile)
            
