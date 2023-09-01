#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd
import os
import glob
import shutil

def read_metrics(file):
    checkm_df = pd.read_csv(file, sep='\t', index_col=0)
    checkm_flag = np.all([checkm_df['Completeness'] > 50, 
                          checkm_df['Contamination'] < 10, 
                          checkm_df['Completeness']-5*checkm_df['Contamination'] > 50
                         ], axis=0)
    return checkm_df, checkm_flag


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter bins based on checkM stats')
    parser.add_argument('checkm1', help='CheckM1 stats')
    parser.add_argument('checkm2', help='CheckM2 stats')
    parser.add_argument('input', help='Input bin folder')
    parser.add_argument('output', help='Output folder')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        checkm1 = args.checkm1
        checkm2 = args.checkm2
        src_dir = args.input
        trg_dir = args.output
        
        if not os.path.exists(trg_dir):
            os.mkdir(trg_dir)
        else:
            for f in glob.glob(trg_dir+'*.fa'):
                try:
                    os.remove(f)
                except OSError as e:
                    print("Error: %s : %s" % (f, e.strerror))
        
        checkm1_df, checkm1_flag = read_metrics(checkm1)
        checkm1_df['pass.checkm1'] = checkm1_flag
        checkm2_df, checkm2_flag = read_metrics(checkm2)
        checkm2_df['pass.checkm2'] = checkm2_flag
        stats_df = pd.concat([checkm1_df, checkm2_df], axis=1)
        qualify = stats_df.loc[np.any([stats_df['pass.checkm1'] == True, 
                                       stats_df['pass.checkm2'] == True], axis=0)]
        for bin_name in qualify.index:
            src = src_dir+'/'+bin_name+'.fa'
            shutil.copy(src, trg_dir+'/'+bin_name+'.fa')

