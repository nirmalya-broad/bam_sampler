#!/usr/bin/env python

import os
import argparse
import subprocess
import ntpath

from subprocess import call 


class BamSamplerMain:
    def __init__(self, args):
        self.infile = args.infile
        self.outdir = args.outdir
        self.run_type = args.run_type
        self.steps_num = args.steps_num
        self.depth_p = args.depth_p
        self.repeat_num = args.repeat_num
        self.main_seed = args.main_seed

        self.samtools_path = "/broad/IDP-Dx_work/nirmalya/local/bin/samtools"
        self.remove_bu_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/remove_bu"        
        self.frag_count_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/frag_count"
        self.seed_gen_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/seed_gen"

        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        temp_dir = outdir + "/temp_dir"
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        self.temp_dir = temp_dir

        self.delete_mainbam = False

    def get_frag_count(self, lfile):
        frag_count_path = self.frag_count_path

        frag_count_cmd = frag_count_path + " -i " + lfile
        print("frag_count_cmd: " + frag_count_cmd)
        proc = subprocess.Popen(frag_count_cmd.split(), stdout=subprocess.PIPE)
        for line1 in proc.stdout.readlines():
            line2 = line1.rstrip()
            if '@frag_count' in line2:
                parts = line2.split()
                lcount = int(parts[1])
                return lcount
        return -1

    def is_query_sorted(self, lfile):
        # Check if the bam file is sorted by coordinates
        samtools_path = self.samtools_path

        samtools_cmd = samtools_path + " view -H " + lfile
        proc = subprocess.Popen(samtools_cmd.split(), stdout=subprocess.PIPE)
        for line1 in proc.stdout.readlines():
            line2 = line1.rstrip()
            if '@HD' in line2:
                parts = line2.split()
                for lpart in parts:
                    if 'SO:queryname' in lpart:
                        return True
        return False
        
    def sort_by_queryname(self, infile, outdir):
        infile_name = ntpath.basename(infile)
        infile_name_u = None
        if infile_name.endswith('_pe.bam'):
            infile_name_u = infile_name.replace('_pe.bam', '_u.bam')
        else:
            infile_name_u = infile_name.replace('.bam', '_u.bam')
        outfile = outdir + "/" + infile_name_u
        samtools_path = self.samtools_path
        self.sort_by_qname_internal(samtools_path, infile, outfile)
        return outfile

    def symlink_outdir(self, infile, outdir):
        infile_name = ntpath.basename(infile)

        outfile = outdir + "/" + infile_name
        os.symlink(infile, outfile)
        
    def sort_by_qname_internal(self, samtools_path, inbam, outbam):
        samcmd = samtools_path + ' sort -n -o ' + outbam + ' ' + inbam
        print("sort_qname_cmd: " + samcmd)
        #call(samcmd.split())

    def get_mapped_qsorted(self, infile_qsorted):
        remove_bu_path = self.remove_bu_path 
        infile_dirname = os.path.dirname(infile_qsorted)
        infile_basename = os.path.basename(infile_qsorted)
        mapped_name = infile_basename.replace('.bam', '_m.bam')
        unmapped_name = infile_basename.replace('.bam', '_um.bam')
        mapped_qsorted = infile_dirname + "/" + mapped_name
        unmapped_qsorted = infile_dirname + "/" + unmapped_name
        remove_bu_cmd = remove_bu_path + " -i " + infile_qsorted + " -o " +\
            mapped_qsorted + " -u " + unmapped_qsorted
        print("remove_bu_cmd started: " + remove_bu_cmd)
        #call(remove_bu_cmd.split()) 
        print("remove_bu_cmd end")
        self.infile_qsorted_mapped = mapped_qsorted
        return mapped_qsorted
        
    def exe_seed_gen(self):
        seed_gen_path = self.seed_gen_path
        infile_qsorted_mapped = self.infile_qsorted_mapped
        main_seed = self.main_seed
        run_type = self.run_type
        steps_num = self.steps_num
        depth_p = self.depth_p
        repeat_num = self.repeat_num

        seed_table = infile_qsorted_mapped.replace(".bam", "_seed.txt")
        seed_gen_cmd = seed_gen_path + " -i " + infile_qsorted_mapped +\
             " -o " + seed_table + " -m " + str(main_seed) +\
             " -r " + run_type + " -s " + str(steps_num) +\
             " -d " + str(depth_p) + " --repeat_num " + str(repeat_num)
        print("seed_gen_cmd: " + seed_gen_cmd)
        call(seed_gen_cmd.split())
        print("End of seed_gen")

    def mainFunc(self):
        # Check if the input sam/bam file is sorted by query/read name
        # If the sam/bam file is not sorted by query, then sort by query
        
        infile = self.infile
        outdir = self.outdir
        temp_dir = self.temp_dir

        if self.is_query_sorted(infile):
            print("The infile is sorted by query.")
            self.infile_qsorted = self.symlink_outdir(infile, temp_dir)
        else:
            print("The infile is not sorted by query.")
            print("infile: " + infile)
            print("temp_dir: " + temp_dir)
            self.infile_qsorted = self.sort_by_queryname(infile, temp_dir)

        # get the reads that are mapped
        infile_qsorted = self.infile_qsorted
        infile_qsorted_mapped = self.get_mapped_qsorted(infile_qsorted)

        # Get the number of fragments in the infile_qsorted_mapped
        print("Getting count of: " + infile_qsorted_mapped)
        lfrag_count = self.get_frag_count(infile_qsorted_mapped)
        print("frag_count of infile_qsorted_mapped: " + str(lfrag_count))

        # Run the seed generator from infile_qsorted_mapped
        self.exe_seed_gen()

        # Run each jobs subjobs for sampling followed by counting
        # Also inside each job do PCR collapse and then counting
        self.exe_sample_bam_core()

        # Finally prepare the count table.
        # Do count table; one for before the PCR collapse and another is
        # after the PCR collapse.



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for bam sampler.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--infile", dest = "infile", type = str, required = True, help = "Location of the input sam/bam file.")
    parser.add_argument("--outdir", dest = "outdir", type = str, required = True, help = "Location of the output directory.")
    parser.add_argument("--run_type", dest = "run_type", type = str, required = True, help = "Type of the run, steps/depth.")
    parser.add_argument("--steps_num", dest = "steps_num", type = int, default = -1, help = "Number of steps when run_type is steps")
    parser.add_argument("--depth_p", dest = "depth_p", type = int, default = -1, help = "Percentage of depth when run_type is depth")
    parser.add_argument("--repeat_num", dest = "repeat_num", type = int, default = 1, help = "Number of repetition." )
    parser.add_argument("--main_seed", dest = "main_seed", type = int, default = 12345, help = "Value of main seed")
    
    args = parser.parse_args()
    bsmo = BamSamplerMain(args)
    print("About to start bsmo mainFunc")
    bsmo.mainFunc()
    


