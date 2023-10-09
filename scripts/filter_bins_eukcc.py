#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd
import os
import glob
import shutil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter fungal bins based on EukCC stats')
    parser.add_argument('eukcc', help='EukCC stats')
    parser.add_argument('input', help='Input bin folder')
    parser.add_argument('output', help='Output folder')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        eukcc = args.eukcc
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
        
        eukcc_df = pd.read_csv(eukcc, sep='\t', index_col=0)
        eukcc_flag = np.all([eukcc_df['completeness'] > 50, 
                             eukcc_df['contamination'] < 10], axis=0)
        qualify = eukcc_df.loc[eukcc_flag]
        for bin_name in qualify.index:
            src = src_dir+'/'+bin_name
            shutil.copy(src, trg_dir+'/'+bin_name)

