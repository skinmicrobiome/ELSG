#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd
import os
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select eukaryotic viral contigs')
    parser.add_argument('contigs', help='Contig file')
    parser.add_argument('blast_file', help='BALSTn results file')
    parser.add_argument('output', help='Output viral contig file')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        contig_file = args.contigs
        blast_file = args.blast_file
        output = args.output
        
        virus_list = ['betaherpesvirus 5', 'betaherpesvirus 6', 'gammaherpesvirus', 
                      'torque teno virus', 
                      'polyomavirus 5', #alpha polyomavirus
                      'polyomavirus 3', 'polyomavirus 4', #beta polyomavirus
                      'polyomavirus 10', 'polyomavirus 11', #delta polyomavirus
                      'alphapapillomavirus',
                      'betapapillomavirus', 'gammapapillomavirus', #HPV
                      'human mastadenovirus', 'molluscum contagiosum virus', 
                      'pepper chlorotic spot',
                      'tobacco vein'
                     ]
        blastn_df = pd.read_csv(blast_file, sep='\t', header=None)
        blastn_top_df = blastn_df.drop_duplicates(keep='first', subset=[0,1])
        blastn_top_df.columns = ['query ID','ref ID','pident','length','mismatch',
                                 'gapopen','qstart','qend','sstart','send','evalue',
                                 'bitscore','qlen','slen','staxids','sscinames','stitle']
        find_list = []
        for target_virus in virus_list:
            virus_flag = [(target_virus.lower() in i.lower()) for i in blastn_top_df['stitle']]
            viral_df = blastn_top_df.loc[np.all([virus_flag, 
                                                 blastn_top_df['pident'] > 90, 
                                                 blastn_top_df['length'] > 1000,
                                                 blastn_top_df['length']/blastn_top_df['qlen'] > 0.7,
                                            ], axis=0)]
            if len(viral_df) > 0:
                for row in viral_df.iterrows():
                    find_list.append((target_virus, 
                                      row[1]['stitle'], row[1]['query ID'], row[1]['pident'], 
                                      row[1]['length'], row[1]['qlen'], row[1]['slen']))
        find_df = pd.DataFrame(find_list)
        find_df.columns = ['target virus','stitle','query ID','pident','length',
                           'qlen','slen']
        find_df = find_df.drop_duplicates(keep='first', subset=['query ID','target virus'])
        contig_list = list(find_df['query ID'])
        find_df = find_df.set_index('query ID')
        seq_lst = []
        for seq_record in SeqIO.parse(contig_file, format = "fasta"):
            if seq_record.id in contig_list:
                seq_record.description += ' '+find_df.loc[seq_record.id, 'target virus'].replace(' ', '_')
                seq_lst.append(seq_record)
        SeqIO.write(seq_lst, output, "fasta")

