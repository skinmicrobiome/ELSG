#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd
import os
import glob
import shutil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter bins based on GUNC stats')
    parser.add_argument('gunc', help='CheckM1 stats')
    parser.add_argument('input', help='Input bin folder')
    parser.add_argument('output', help='Output folder')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        gunc = args.gunc
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
        
        gunc_df = pd.read_csv(gunc, sep='\t', index_col=0)
        chimera_flag = np.any([np.all([gunc_df['contamination_portion'] > 0.05, 
                                       gunc_df['clade_separation_score'] > 0.45, 
                                       gunc_df['reference_representation_score'] > 0.5], axis=0), 
                               gunc_df['pass.GUNC'].isnull()], axis=0)
        gunc_df['pass.GUNC.lenient'] = ~chimera_flag
        chimeric_rate = sum(chimera_flag)/len(gunc_df)
        print('Chimerism rate:', chimeric_rate)
        qualify = gunc_df.loc[gunc_df['pass.GUNC.lenient'] == True]
        for bin_name in qualify.index:
            src = src_dir+'/'+bin_name+'.fa'
            shutil.copy(src, trg_dir+'/'+bin_name+'.fa')

