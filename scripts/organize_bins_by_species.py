#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd
import os
import glob
import shutil
import fnmatch
from collections import Counter

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Organize bins based on GTDB classification')
    parser.add_argument('folder', help='Folder containing bins')
    parser.add_argument('gtdb', help='Path to GTDB-tk results')
    parser.add_argument('metric', help='Bin metric file')
    parser.add_argument('prokka', help='Prokka output folder')
    parser.add_argument('output', help='Output folder')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        
        # Load bins
        genome_folder = args.folder
        genomes = sorted(fnmatch.filter(os.listdir(genome_folder), '*.fa'))
        genomes = [i.replace('.fa','') for i in genomes]
        
        # Read GTDB classifications
        gtdb_file = args.gtdb
        gtdb_df = pd.read_csv(gtdb_file, sep='\t', index_col=0)
        gtdb_df['phylum'] = [i.split(';p__')[1].split(';')[0] if ';p' in i else i for i in gtdb_df['classification']]
        gtdb_df['species'] = [i.split(';s__')[1] if ';s' in i else i for i in gtdb_df['classification']]
        gtdb_df = gtdb_df.loc[[g for g in genomes if g in gtdb_df.index]]
        
        # Read bin metrics
        checkm_file = args.metric
        checkm_df = pd.read_csv(checkm_file, sep='\t', index_col=0)
        checkm_df = checkm_df.loc[genomes]
        
        # Select qualified bins
        info_df = pd.concat([checkm_df, gtdb_df], axis=1)
        info_df = info_df.loc[np.all([info_df['Completeness'] > 90, info_df['Contamination'] < 5], axis=0)]
        species_counts = dict(Counter(info_df['species']))
        info_df['species_count'] = [species_counts[i] for i in info_df['species']]
        info_df = info_df.loc[np.all([info_df['species_count']>=2, info_df['species']!=''], axis=0)]
        qualify_species = np.unique(info_df['species'])
        print('Qualified species:', len(qualify_species))
        
        # Create folder and save files for each species
        target_dir = args.output
        if not os.path.exists(target_dir):
            os.mkdir(target_dir)
        src_dir = args.prokka
        for sp in qualify_species:
            sp_folder = 'qualifiedSpecies_'+sp.replace(' ','_')
            if not os.path.exists(target_dir+'/'+sp_folder):
                os.mkdir(target_dir+'/'+sp_folder)
            else:
                for f in glob.glob(target_dir+'/'+sp_folder+'/*.gff'):
                    try:
                        os.remove(f)
                    except OSError as e:
                        print("Error: %s : %s" % (f, e.strerror))
                for f in glob.glob(target_dir+'/'+sp_folder+'/*.faa'):
                    try:
                        os.remove(f)
                    except OSError as e:
                        print("Error: %s : %s" % (f, e.strerror))
            genome_list = info_df.loc[info_df['species']==sp].index
            print(sp_folder, len(genome_list))
            for genome in genome_list:
                shutil.copy(src_dir+'/'+genome+'/'+genome+'.gff', target_dir+'/'+sp_folder+'/'+genome+'.gff')
                shutil.copy(src_dir+'/'+genome+'/'+genome+'.faa', target_dir+'/'+sp_folder+'/'+genome+'.faa')

