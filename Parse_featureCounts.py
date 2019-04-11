#!/usr/bin/env python
import pandas as pd
import csv
import ntpath
import re
import os
import os.path


def get_tab_counts(infile_str, tab_tag, series_name, numskip = 0, has_header = True):
    lfile = open(infile_str, 'r')
    tab1 = csv.reader(lfile, delimiter = '\t')
    for _ in xrange(numskip):
        next(tab1)
    if has_header:
        header = tab1.next()
    tab2 = [row for row in tab1]
    tab3 = pd.DataFrame(tab2, columns = tab2[0])
    geneid_series = tab3.iloc[:,0]
    count_series = tab3.iloc[:,-1]
    geneid_series.name = tab_tag
    count_series.name = series_name
    return geneid_series, count_series

def get_all_metrics_counts(result_dir, outfile, prefix_lst, ref_acc = "", \
        ldelim = '/'):
    tab_tag = 'metrics_type'
    header_lst = [tab_tag] + prefix_lst
    df = pd.DataFrame(columns=prefix_lst)
    geneid_series_ori = pd.Series()
    count_series_ori = pd.Series()
    for lfile_str in prefix_lst:
        lfile_path = ''
        if ref_acc:
            lfile_path = result_dir + ldelim + lfile_str + "_" + ref_acc + \
                ".metrics"
        else:
            lfile_path = result_dir + ldelim + lfile_str + ".metrics"
        geneid_series, count_series = get_tab_counts(lfile_path, tab_tag, lfile_str)  
        if geneid_series_ori.empty:
            geneid_series_ori = geneid_series
            count_series_ori = count_series
        if not geneid_series_ori.equals(geneid_series):
            raise ValueError('These two geneid_series are not equal: ' + \
                count_series_ori.name + ' and ' + count_series.name)
        df[count_series.name] = count_series
  
    if not os.path.exists(os.path.dirname(outfile)):
        try:
            os.makedirs(os.path.dirname(outfile))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    df_t = df.transpose()
    df_t.insert(0, 'sample_id', prefix_lst)
    header_lst = ['sample_id'] + list(geneid_series_ori.values)
    df_t.header = header_lst

    df_t.columns = header_lst
    df_t.to_csv(outfile, sep = '\t', index = False, header = True)

def get_all_gene_counts(result_dir, outfile, prefix_lst, ref_acc = "", \
        ldelim = '/', has_header = True):
    tab_tag = 'Geneid'
    header_lst = [tab_tag] + prefix_lst
    df = pd.DataFrame(columns=header_lst)
    geneid_series_ori = pd.Series()
    count_series_ori = pd.Series()
    for lfile_str in prefix_lst:
        lfile_path = ''
        if ref_acc:
            lfile_path = result_dir + ldelim + lfile_str + "_" + ref_acc + \
                ".counts"
        else:
            lfile_path = result_dir + ldelim + lfile_str + ".counts"
        geneid_series, count_series = get_tab_counts(lfile_path, tab_tag, lfile_str, has_header = has_header)
        if geneid_series_ori.empty:
            geneid_series_ori = geneid_series
            count_series_ori = count_series
        if not geneid_series_ori.equals(geneid_series):
            raise ValueError('These two geneid_series are not equal: ' + count_series_ori.name + ' and ' + count_series.name)
        df[count_series.name] = count_series
    df[geneid_series_ori.name] = geneid_series_ori
    if not os.path.exists(os.path.dirname(outfile)):
        try:
            os.makedirs(os.path.dirname(outfile))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    df.to_csv(outfile, sep = '\t', index = False)
