#!/usr/bin/env python

import os
import argparse
import subprocess
import ntpath
import shutil

from subprocess import call 
from shutil import copyfile

from sbcore import SBCore
from Parse_featureCounts import get_all_gene_counts
from Parse_featureCounts import get_all_metrics_counts


class BamSamplerMain:
    def __init__(self, args):
        self.infile = args.infile
        self.outdir = args.outdir
        self.run_type = args.run_type
        self.steps_num = args.steps_num
        self.depth_p = args.depth_p
        self.repeat_num = args.repeat_num
        self.main_seed = args.main_seed
        self.use_qsub = args.use_qsub
        self.add5 = args.add5
        self.add3  = args.add3
        self.patho_id = args.patho_id
        self.project_id = args.project_id
        self.do_ref = args.do_ref
        self.use_qsort = args.use_qsort
        self.do_remove_bu = args.do_remove_bu
        self.start_read = args.start_read
        self.do_metrics = args.do_metrics

        basepath = os.path.dirname(os.path.realpath(__file__))
        self.Script_dir = basepath
        self.SBCore = self.Script_dir + "/" + "sbcore.py"


        self.samtools_path = "/broad/IDP-Dx_work/nirmalya/local/bin/samtools"
        self.remove_bu_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/remove_bu"        
        self.frag_count_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/frag_count"
        self.seed_gen_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/seed_gen"
        self.sample_bam_path = "/broad/IDP-Dx_work/nirmalya/research/bam_sampler/sample_bam"
        self.UGER_cbp = "/broad/IDP-Dx_work/nirmalya/tools/ugetools/UGE_SUBMISSIONS/UGER_cmd_batch_processor.py"
        self.gff_parse = "/broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/gff_parse3.sh"    
        self.patho_dbpath = "/broad/IDP-Dx_storage/NCBI_files2/"    

        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        initdir = outdir + "/initdir"
        if not os.path.exists(initdir):
            os.makedirs(initdir)
        self.initdir = initdir

        datadir = outdir + "/datadir"
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        self.datadir = datadir

        bamdir = outdir + "/bamdir"
        if not os.path.exists(bamdir):
            os.makedirs(bamdir)
        self.bamdir = bamdir

        logdir = outdir + "/logdir"
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        self.logdir = logdir

        summary_dir = outdir + "/summary_dir"
        if not os.path.exists(summary_dir):
            os.makedirs(summary_dir)
        self.summary_dir = summary_dir

        UGER_cbp_dir = outdir + "/UGER_cbp"

        if self.use_qsub and os.path.exists(UGER_cbp_dir):
            shutil.rmtree(UGER_cbp_dir)
        self.UGER_cbp_dir = UGER_cbp_dir

        temp_bamdir = bamdir + "/temp_bamdir"
        if not os.path.exists(temp_bamdir):
            os.makedirs(temp_bamdir)

        nodupdir = bamdir + "/nodupdir"
        if not os.path.exists(nodupdir):
            os.makedirs(nodupdir)

        nodupdir_temp_bamdir = nodupdir + "/temp_bamdir"
        if not os.path.exists(nodupdir_temp_bamdir):
            os.makedirs(nodupdir_temp_bamdir)

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
        use_qsort = self.use_qsort
        infile_name = ntpath.basename(infile)
        infile_name_u = None
        if infile_name.endswith('_pe.bam'):
            infile_name_u = infile_name.replace('_pe.bam', '_u.bam')
        else:
            infile_name_u = infile_name.replace('.bam', '_u.bam')
        outfile = outdir + "/" + infile_name_u
        samtools_path = self.samtools_path
        if use_qsort:
            self.sort_by_qname_internal(samtools_path, infile, outfile)
        return outfile

    def symlink_outdir(self, infile, outdir):
        infile_name = ntpath.basename(infile)

        outfile = outdir + "/" + infile_name
        os.symlink(infile, outfile)
        
    def sort_by_qname_internal(self, samtools_path, inbam, outbam):
        samcmd = samtools_path + ' sort -n -o ' + outbam + ' ' + inbam
        print("sort_qname_cmd: " + samcmd)
        call(samcmd.split())

    def get_mapped_qsorted(self, infile_qsorted):
        do_remove_bu = self.do_remove_bu
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
        if do_remove_bu:
            call(remove_bu_cmd.split()) 
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
        start_read = self.start_read
        infile_frag_count = self.infile_frag_count

        print("stard_read: " + str(start_read))
        print("infile_frag_count: " + str(infile_frag_count))
        start_read_p = (start_read * 100.0) / infile_frag_count
        print("start_read_p: " + str(start_read_p))

        seed_table = infile_qsorted_mapped.replace(".bam", "_seed.txt")
        seed_gen_cmd = seed_gen_path + " -i " + infile_qsorted_mapped +\
             " -o " + seed_table + " -m " + str(main_seed) +\
             " -r " + run_type + " -s " + str(steps_num) +\
             " -d " + str(depth_p) + " --repeat_num " + str(repeat_num) +\
             " --start_read_p " + str(start_read_p)

        print("seed_gen_cmd: " + seed_gen_cmd)
        call(seed_gen_cmd.split())
        print("End of seed_gen")
        self.seed_table = seed_table
        return seed_table

    def get_sample_bam_core_cmd(self, lsuffix, ltop_seed, lsample_p):
        infile_str = self.infile_qsorted_mapped
        bamdir = self.bamdir
        patho_id = self.patho_id
        datadir = self.datadir
        SBCore = self.SBCore
        SB_cmd = SBCore + " --infile " + infile_str + " --bamdir " + bamdir +\
            " --suffix " + lsuffix + " --sample_p " + lsample_p +\
            " --top_seed " + ltop_seed + " --patho_id " + patho_id +\
            " --datadir " + datadir
        return SB_cmd
        
        
        
    def exe_sample_bam_core(self, seed_table):

        logdir = self.logdir
        use_qsub = self.use_qsub
        sbcore_job_path = logdir + "/" + "sbcore_joblist_txt"
        num_cores = 1
        lmemory = 8
        UGER_cbp = self.UGER_cbp
        UGER_cbp_dir = self.UGER_cbp_dir
        lsuf_lst = list()

        jfile = open(sbcore_job_path, "w")
        with open(seed_table) as stab:
            for line1 in stab:
                line2 = line1.rstrip()
                parts = line2.split()
                lsample = parts[0]
                lrepeat = parts[1]
                lsample_p = parts[2]
                ltop_seed = parts[3]
                lsuffix = lsample + "_" + lrepeat
                lsuf_lst.append(lsuffix)
                print(lsuffix)
                main_seed = self.main_seed
                cmd_str = self.get_sample_bam_core_cmd(lsuffix, ltop_seed,\
                    lsample_p)
                loutfile = logdir + "/" + lsuffix + "_out.txt"
                lerrfile = logdir + "/" + lsuffix + "_err.txt"
                cmd_str2 = cmd_str + " 1> " + loutfile  + " 2> " + lerrfile
                print(cmd_str2) 
                jfile.write(cmd_str2 + "\n")
        jfile.close()

        self.lsuf_lst = lsuf_lst

        joblist_cmd = UGER_cbp + " --cmds_file " + sbcore_job_path + \
                                " --batch_size 1" + \
                                " --num_cores " + str(num_cores) + \
                                " --memory " + str(lmemory) + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"        
        if use_qsub:
            call(joblist_cmd.split())

       
    def exe_gff_parser(self):
        gff_parse = self.gff_parse
        patho_dbpath = self.patho_dbpath
        datadir = self.datadir
        add3 = self.add3
        add5 = self.add5
        patho_id = self.patho_id

        gff_cmd = "sh " + gff_parse + " " + patho_dbpath + " " + datadir + " " + datadir + " " + str(add5) + " " + str(add3) + " " + patho_id
        print("gff_cmd: " + gff_cmd)
        call(gff_cmd.split())

   
    def write_read_count_table(self, lsample_lst, has_header, project_dir, subdir = ''):
        ldelim = "/"
        lproject_id = self.project_id
        ref_acc = self.patho_id
        
        loutdir = self.summary_dir
        if subdir:
            loutdir += ldelim + subdir

        if not os.path.exists(loutdir):
            os.makedirs(loutdir)

        loutfile = loutdir + ldelim +  lproject_id + "_" + ref_acc + "_counts.tsv"
        if subdir:
            project_dir += ldelim + subdir

        get_all_gene_counts(project_dir, loutfile, lsample_lst, ref_acc, has_header = has_header)

    def process_seed_info(self):
        summary_dir = self.summary_dir
        seed_table = self.seed_table

        main_seed_file = summary_dir + "/main_seed.txt"
        with open(main_seed_file, "w") as outf:
            outf.write("main seed:\t" + str(self.main_seed) + "\n")
 
        project_id = self.project_id 
        ref_acc = self.patho_id
        dest_leaf = project_id + "_" + ref_acc + "_seed.txt" 
        seed_table_leaf = os.path.basename(seed_table)
        seed_table_dest = summary_dir + "/" + dest_leaf
        copyfile(seed_table, seed_table_dest)    

 
    def mainFunc(self):
        # Check if the input sam/bam file is sorted by query/read name
        # If the sam/bam file is not sorted by query, then sort by query
        
        infile = self.infile
        outdir = self.outdir
        initdir = self.initdir
        bamdir = self.bamdir
        do_ref = self.do_ref
        use_qsort = self.use_qsort
        do_metrics = self.do_metrics

        if self.is_query_sorted(infile):
            print("The infile is sorted by query.")
            self.infile_qsorted = self.symlink_outdir(infile, initdir)
        else:
            print("The infile is not sorted by query.")
            print("infile: " + infile)
            print("initdir: " + initdir)
            self.infile_qsorted = self.sort_by_queryname(infile, initdir)

        # get the reads that are mapped
        infile_qsorted = self.infile_qsorted
        infile_qsorted_mapped = self.get_mapped_qsorted(infile_qsorted)

        # Get the number of fragments in the infile_qsorted_mapped
        print("Getting count of: " + infile_qsorted_mapped)
        lfrag_count = self.get_frag_count(infile_qsorted_mapped)
        print("frag_count of infile_qsorted_mapped: " + str(lfrag_count))
        self.infile_frag_count = lfrag_count

        # Run the seed generator from infile_qsorted_mapped
        seed_table = self.exe_seed_gen()
        
        # Run gff_parser
        if do_ref:
            self.exe_gff_parser()

        # Run each jobs subjobs for sampling followed by counting
        # Also inside each job do PCR collapse and then counting
        self.exe_sample_bam_core(seed_table)

        # Finally prepare the count table.
        # Do count table; one for before the PCR collapse and another is
        # after the PCR collapse.
        if do_metrics:
            lsuf_lst = self.lsuf_lst
            has_header = False
            self.write_read_count_table(lsuf_lst, has_header, bamdir, subdir = '')
            self.write_read_count_table(lsuf_lst, has_header, bamdir, subdir = 'nodupdir')

            self.process_seed_info()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for bam sampler.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--infile", dest = "infile", type = str, required = True, help = "Location of the input sam/bam file.")
    parser.add_argument("--outdir", dest = "outdir", type = str, required = True, help = "Location of the output directory.")
    parser.add_argument("--project_id", dest = "project_id", type = str, required = True, help = "Project id string")
    parser.add_argument("--run_type", dest = "run_type", type = str, required = True, help = "Type of the run, steps/depth.")
    parser.add_argument("--steps_num", dest = "steps_num", type = int, default = -1, help = "Number of steps when run_type is steps")
    parser.add_argument("--depth_p", dest = "depth_p", type = float, default = -1, help = "Percentage of depth when run_type is depth")
    parser.add_argument("--repeat_num", dest = "repeat_num", type = int, default = 1, help = "Number of repetition." )
    parser.add_argument("--main_seed", dest = "main_seed", type = int, default = 12345, help = "Value of main seed")
    parser.add_argument("--patho_id", dest = "patho_id", type = str, required = True, help = "NCBI ref id for the species")
    parser.add_argument('--ADD5', dest = 'add5', type = int, default = 0, help = 'ADD5 for gff parser')
    parser.add_argument('--ADD3', dest = 'add3', type = int, default = 0, help = 'ADD3 for gff parser')
    parser.add_argument('--start_read', dest = 'start_read', type = int, default = 100, help = 'Read number of starting breakpoint')
    parser.add_argument('--no_ref', dest = 'do_ref', action = 'store_false', default = True, help = 'Does not generate patho ref.')

    parser.add_argument('--no_qsub', dest = 'use_qsub', action = 'store_false', default = True, help = 'Does not submit qsub jobs.' )
    parser.add_argument('--no_qsort', dest = 'use_qsort', action = 'store_false', default = True, help = 'Does not do the infile qsort.' )
    parser.add_argument('--no_remove_bu', dest = 'do_remove_bu', action = 'store_false', default = True, help = 'Does not remove unmapped reads from bam.' )
    parser.add_argument('--no_metrics', dest = 'do_metrics', action = 'store_false', default = True, help = 'Does not do the metrics.')
    
    args = parser.parse_args()
    bsmo = BamSamplerMain(args)
    print("About to start bsmo mainFunc")
    bsmo.mainFunc()
    


