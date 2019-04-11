#!/usr/bin/env python

import os
import argparse
from alignerutils import sort_by_qname
from alignerutils import sort_bam
from alignerutils import bam_to_sam
from counterjl import CounterJL

from subprocess import call

class SBCore:

    def __init__(self, args):
        self.infile = args.infile
        self.bamdir = args.bamdir
        self.suffix = args.suffix
        self.sample_p = float(args.sample_p)
        self.top_seed = int(args.top_seed)
        self.patho_id = args.patho_id

        bamdir = self.bamdir
        self.nodupdir = self.bamdir  + "/nodupdir"
        self.datadir = args.datadir


        self.remove_dup_script = "/broad/IDP-Dx_work/nirmalya/research/read_counter/remove_dup"
        self.samtools_path = "/broad/IDP-Dx_work/nirmalya/local/bin/samtools"
        self.picard_bindir = "/broad/IDP-Dx_work/nirmalya/tools/picard/latest"
        self.sample_bam_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/sample_bam"
        self.JLCounter = "/broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/SAM_to_counts2.sh"
        self.shell_script_dir = "/broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/"

    def mark_dup_reads(self, input_bam):
        # Get the output directory which has been created before
        # Do not remove the duplicates, just let picard mark them using read is PCR or optical 
        # duplicate (0x400). We shall remove them later in pair
        ldelim = "/"
        picard_bindir = self.picard_bindir
        picard_jar = picard_bindir + "/picard.jar"

        bamdir = os.path.dirname(input_bam)
        no_dup_dir = bamdir + ldelim + "nodupdir"
        bam_file = os.path.basename(input_bam)
        output_bam_ori = no_dup_dir + ldelim + bam_file
        # Output bam is the dup marked one
        output_bam = output_bam_ori.replace("_pe.bam", "_dm.bam")

        metrics_txt = output_bam.replace('.bam', '_rm_dup_met.txt')

        # Then rename the old bam file

        mem_str = "-Xmx6G"
        dup_str = "MarkDuplicates"
        picard_cmd = "java " + mem_str + " -jar " + picard_jar + " " + dup_str + " I=" + input_bam +\
            " O=" + output_bam + " M=" + metrics_txt + " REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT"
        print("start: picard_cmd: " + picard_cmd)
        call(picard_cmd.split())
        print("end: picard_cmd:")
        return (output_bam)

    def remove_dup_reads(self, dm_sorted_bam):
        samtools_path = self.samtools_path
        remove_dup_script = self.remove_dup_script

        # 1. sort wrt queryname
        dm_unsorted_bam = dm_sorted_bam.replace("_dm.bam", "_u_dm.bam")
        if not os.path.isfile(dm_unsorted_bam):
            sort_by_qname(samtools_path, dm_sorted_bam, dm_unsorted_bam)
            print("sort by qname: " + dm_sorted_bam + " to " + dm_unsorted_bam)

        # 2. remove the pair of marked reads 0x400
        nodup_unsorted_bam = dm_unsorted_bam.replace("_u_dm.bam", "_u.bam")
        dup_bam = dm_unsorted_bam.replace("_u_dm.bam", "_u_dup.bam")
        # Call the c++ program remove_dup
        remove_dup_cmd = remove_dup_script + " -i " + dm_unsorted_bam +\
            " -o " + nodup_unsorted_bam + " -d " + dup_bam
        print("remove_dup_cmd starts: " + remove_dup_cmd)
        call(remove_dup_cmd.split())
        print("remove_dup_cmd ends")
        
        self.nodup_unsorted_bam = nodup_unsorted_bam

        # 3. sort the no_dup_unsorted_bam bam file
        nodup_sorted_bam = nodup_unsorted_bam.replace("_u.bam", "_pe.bam")
        if not os.path.isfile(nodup_sorted_bam):
            sort_bam(samtools_path, nodup_unsorted_bam, nodup_sorted_bam)
            print("sort by coordinate: " + nodup_unsorted_bam + " to " +\
                nodup_sorted_bam)
        os.remove(dm_unsorted_bam)
        print("Removed dm_unsorted_bam: " + dm_unsorted_bam)
        os.remove(dm_sorted_bam)
        print("Remove dm_sorted_bam: " + dm_sorted_bam)
        return nodup_sorted_bam

        

    def exe_sample_bam(self):
        infile_str = self.infile
        bamdir = self.bamdir
        lsuffix = self.suffix
        ltop_seed = self.top_seed
        lsample_p = self.sample_p
        infile_bname = os.path.basename(infile_str)
        out_suffix = "_" + lsuffix + "_u.bam"
        self.out_suffix = out_suffix
        outfile_bname = infile_bname.replace(".bam", out_suffix)
        out_file = bamdir + "/" + outfile_bname
        sample_bam_path = self.sample_bam_path
        sample_bam_cmd = sample_bam_path + " -i " + infile_str + " -o " + out_file + " -t " + str(ltop_seed) + " -s " + str(lsample_p)
        call(sample_bam_cmd.split())
        self.sampled_bam_u = out_file
        print("sample_bam_cmd: " + sample_bam_cmd)

    def exe_sort_bam(self, sampled_bam_u):
        sampled_bam_pe = sampled_bam_u.replace('_u.bam', '_pe.bam')
        print("Sorting by coordinate: " + sampled_bam_u + " to " + sampled_bam_pe)
        samtools_path = self.samtools_path
        sort_bam(samtools_path, sampled_bam_u, sampled_bam_pe)
        return sampled_bam_pe
        
    def exe_pcr_collapse(self):
        bamdir = self.bamdir
        nodupdir = self.nodupdir
        sampled_bam_u = self.sampled_bam_u
        # Create a new bam file in the no_dup_dir after removal of the 
        # pcr duplicates from pcr; ideally that folder should also 
        # contain the metrics related to dup removal
        sampled_bam_pe = self.exe_sort_bam(sampled_bam_u)
        self.with_dup_unsorted_bam = sampled_bam_u
        dup_marked_sorted_bam = self.mark_dup_reads(sampled_bam_pe)
        nodup_sorted_bam = self.remove_dup_reads(dup_marked_sorted_bam)
        self.nodup_sorted_bam = nodup_sorted_bam

    def count_paired(self, unsorted_bam):
        JLCounter = self.JLCounter
        ldelim = "/"
        ref_acc = self.patho_id
        shell_script_dir = self.shell_script_dir
        suffix = self.suffix
        samtools_path = self.samtools_path
         
        # get unsorted_sam
        unsorted_sam = unsorted_bam.replace('.bam', '.sam')
        bam_to_sam(samtools_path, unsorted_bam, unsorted_sam)
        outdir = os.path.dirname(unsorted_sam)
        Data_dir = self.datadir       
        patho_temp_bamdir = outdir + ldelim + "temp_bamdir"
        patho_gff = Data_dir + ldelim + ref_acc + "_ALL.gff"
        countfile_str = outdir + ldelim + suffix + "_" + ref_acc + ".counts"
        count_cmd = "sh " + JLCounter + " " + unsorted_sam + " " + patho_gff + \
            " " + shell_script_dir + " " + patho_temp_bamdir + " " + \
            countfile_str
        print("Starting JL counter for RtS: " + count_cmd)
        call(count_cmd.split())
        os.remove(unsorted_sam)


    def mainFunc(self):
        # 1. Do the actual sampling given the top_seed and percentage
        self.exe_sample_bam() 
        # Do the pcr collapse
        self.exe_pcr_collapse()

        with_dup_unsorted_bam = self.with_dup_unsorted_bam
        print("with_dup_unsorted_bam: " + with_dup_unsorted_bam)
        self.count_paired(with_dup_unsorted_bam)

        nodup_unsorted_bam = self.nodup_unsorted_bam
        print("nodup_unsorted_bam: " + nodup_unsorted_bam)
        self.count_paired(nodup_unsorted_bam)
        
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for a sample bam core instance", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--infile", dest = "infile", type = str, required = True, help = "Input file bam")    
    parser.add_argument("--bamdir", dest ="bamdir", type = str, required = True, help = "Output directory")
    parser.add_argument("--datadir", dest ="datadir", type = str, required = True, help = "Data directory")
    parser.add_argument("--suffix", dest = "suffix", type = str, required = True, help = "Suffix directory")
    parser.add_argument("--sample_p", dest = "sample_p", type = str, required = True, help = "Suffix directory")
    parser.add_argument("--top_seed", dest = "top_seed", type = str, required = True, help = "Top seed for the sampling")
    parser.add_argument("--patho_id", dest = "patho_id", type = str, required = True, help = "NCBI ref id")

    args = parser.parse_args()
    sbco = SBCore(args)
    print("About to start sbco mainFunc")
    sbco.mainFunc()
    
    
