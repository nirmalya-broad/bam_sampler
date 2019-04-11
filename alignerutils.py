import os
import re
import os.path
from miscutils import exe_command
from subprocess import call
from miscutils import copyLargeFile

def sam_to_bam_sorted(outsamfile, samtools_path, tdf_path, bamdir, suf):
	outbamfile = outsamfile.replace('.sam', '_u.bam')
	sam_to_bam(samtools_path, outsamfile, outbamfile)
	suf_str = "_" + suf + ".bam"
	outsorted = outbamfile.replace('_u.bam', suf_str)
	sort_bam(samtools_path, outbamfile, outsorted)
	copy_bam(samtools_path, bamdir, outsorted)
	#os.remove(outbamfile)
	return outsorted

def sam_to_bam(samtools_path, outsamfile, outbamfile):
	bamcmd = samtools_path + ' view -b -S ' + outsamfile
	exe_command(bamcmd, outbamfile)

def bam_to_sam(samtools_path, outbamfile, outsamfile):
    samcmd = samtools_path + ' view -h ' + outbamfile
    print("Converting bam to sam: " + outbamfile)
    print("samcmd: " + samcmd)
    exe_command(samcmd, outsamfile)

def sort_by_qname(samtools_path, inbam, outbam):
    samcmd = samtools_path + ' sort -n -o ' + outbam + ' ' + inbam
    print("sort_qname_cmd: " + samcmd)
    call(samcmd.split())

def sort_bam(samtools_path, outbamfile, outsorted):
	sortedcmd = samtools_path + ' sort -o ' + outsorted + " " + outbamfile
	call(sortedcmd.split())

def copy_bam(samtools_str, bam_path, outsorted):
	ldelim = "/"
	# Copy the bam file to bam_path
	if not os.path.exists(bam_path):
		os.makedirs(bam_path)

	out_bam_base = os.path.basename(os.path.normpath(outsorted))
	copy_bam_path = bam_path + ldelim + out_bam_base
	copyLargeFile(outsorted, copy_bam_path)
	bai_cmd = samtools_str + " index " + copy_bam_path
	call(bai_cmd.split())
	print("copy_bam_path: " + copy_bam_path)

def do_tdf(tdf_str, copy_bam_path):
    tdf_cmd = "java -jar " + tdf_str + " " + copy_bam_path
    call(tdf_cmd.split())
    print("tdf_cmd: " + tdf_cmd)



