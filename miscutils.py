import os
import sys
import shutil
import gzip
import os.path
from subprocess import call

def touch(path):
    with open(path, 'a'):
        os.utime(path, None)

def open_any(lfile):
    lfile1 = None

    lfile_gz = lfile + ".gz"
    if os.path.exists(lfile):
        lfile1 = lfile
    elif os.path.exists(lfile_gz):
        lfile1 = lfile_gz
    else:
        lstr = "Warning! File does not exists: " + lfile
        print(lstr)
        #raise StandardError(lstr);
        return None

    print("Opening file: " + lfile1)
    lfile_r = os.path.realpath(lfile1)
    print("Opening file realpath: " + lfile_r)
    if lfile_r.endswith("gz"):
        return gzip.open(lfile_r, 'rb')
    else:
        return open(lfile_r, "rb")

def file_concat(source_lst, dest_file):
    # Concatenate files in source_lst to des_file in order
    copy_len = 1024*1024*10

    #for lpath in source_lst:
    #    touch(lpath)

    with open(dest_file,'wb') as wfd:
        for f in source_lst:
            fd = open_any(f)
            if fd is not None:
                shutil.copyfileobj(fd, wfd, copy_len)

def exe_command(cmd_str, outfile):
    cmd_lst = cmd_str.split()
    f_out = open(outfile, 'w')
    res = call(cmd_lst, stdout = f_out)
    f_out.close()
    return

def exe_command_stderr(cmd_str, outfile):
    cmd_lst = cmd_str.split()
    f_out = open(outfile, 'w')
    res = call(cmd_lst, stderr = f_out)
    f_out.close()
    return

def copyLargeFile(src, dest, buffer_size=16000):
    print("src: " + src + " dest: " + dest)
    with open(src, 'rb') as fsrc:
        with open(dest, 'wb') as fdest:
            shutil.copyfileobj(fsrc, fdest, buffer_size)

def get_pct(numer, denom):
    if denom == 0:
        return 0.0
    else:
        lval = (numer * 100.0) / denom
        return lval

