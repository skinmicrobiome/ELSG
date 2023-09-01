import os
import fnmatch

parent_folder = '/data/shenz4/Vitality/data/MAG/'
run_list = sorted(fnmatch.filter(os.listdir(parent_folder), 'V*'))
run_list = [r for r in run_list if 'St' not in r]
for run in run_list[:]:
    try:
        prok_run_folder = parent_folder+run+'/03_binning_analyses/refined_bins/pass_GUNC_virusRemoved/'
        prok_bin_files = fnmatch.filter(os.listdir(prok_run_folder), '*.fa')
        prok_bin_files = [prok_run_folder+f for f in prok_bin_files]
    except:
        prok_bin_files = []
    try:
        euk_run_folder = parent_folder+run+'/03_binning_analyses/eukcc_results/predereplicated_virusRemoved/'
        euk_bin_files = fnmatch.filter(os.listdir(euk_run_folder), '*.fa')
        euk_bin_files = [euk_run_folder+f for f in euk_bin_files]
    except:
        euk_bin_files = []
    bin_files = prok_bin_files+euk_bin_files
    if len(bin_files) == 0:
        print('Skipping:', run)
        continue
    n50_list = []
    for bin_file in prok_bin_files:
        contig_lengths = []
        for seq_record in SeqIO.parse(bin_file, format = "fasta"):
            contig_lengths.append(len(seq_record.seq))
        contig_lengths = sorted(contig_lengths)[::-1]
        half_length = sum(contig_lengths)/2
        cur_length = 0
        for l in contig_lengths:
            cur_length += l
            if cur_length > half_length:
                n50_list.append((run+'_prok_'+bin_file.split('/')[-1], l, len(contig_lengths)))
                break
    for bin_file in euk_bin_files:
        contig_lengths = []
        for seq_record in SeqIO.parse(bin_file, format = "fasta"):
            contig_lengths.append(len(seq_record.seq))
        contig_lengths = sorted(contig_lengths)[::-1]
        half_length = sum(contig_lengths)/2
        cur_length = 0
        for l in contig_lengths:
            cur_length += l
            if cur_length > half_length:
                n50_list.append((run+'_euk_'+bin_file.split('/')[-1], l, len(contig_lengths)))
                break
    n50_df = pd.DataFrame(n50_list)
    n50_df.columns = ['Bin', 'N50', 'Contigs']
    n50_df.to_csv(parent_folder+run+'/03_binning_analyses/N50.txt', sep='\t', index=False)
