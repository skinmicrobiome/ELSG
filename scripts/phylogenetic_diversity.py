import dendropy
import pandas as pd

tree_file = 'ELSG.dereplicated.iqtree.treefile'
tree = dendropy.Tree.get_from_path(tree_file, "newick")
clade_lst = []
for n in tree.nodes():
    if n.taxon is None:
        continue
    else:
        clade_lst.append(n.taxon.label)
gtdb_file = 'gtdbtk.bac120.summary.tsv'
gtdb_df = pd.read_csv(gtdb_file, sep='\t', index_col=0)
gtdb_df['phylum'] = [i.split(';p__')[1].split(';')[0] if ';p' in i else i for i in gtdb_df['classification']]
gtdb_df['genus'] = [i.split(';g__')[1].split(';')[0] if ';g' in i else i for i in gtdb_df['classification']]
gtdb_df['species'] = [i.split(';s__')[1] if ';s' in i else i for i in gtdb_df['classification']]
novelty_list = []
for i in gtdb_df.index:
    if i in novel_set:
        if i in novel_set2:
            if gtdb_df.loc[i,'species']=='':
                novelty_list.append('3')
            else:
                novelty_list.append('2')
        else:
            novelty_list.append('1')
    else:
        novelty_list.append('0')
gtdb_df['novel'] = novelty_list
gtdb_df.index = [i.replace('+',' ').replace('_',' ') for i in gtdb_df.index]
gtdb_df = gtdb_df.loc[clade_lst]
print(gtdb_df.shape)
print(Counter(gtdb_df['phylum']))

pd_lst = []
for p in np.unique(gtdb_df['phylum']):
    all_genomes = gtdb_df.loc[gtdb_df['phylum']==p].index.values
    p_tree = tree.extract_tree_with_taxa_labels(labels=all_genomes)
    cur_lst = [p, p_tree.length()]
    for novel_lvl in [0,1,2]:
        if novel_lvl == 0:
            known_genomes = gtdb_df.loc[np.all([gtdb_df['phylum']==p, 
                                                gtdb_df['novel']==str(novel_lvl)], axis=0)].index.values
        else:
            known_genomes = gtdb_df.loc[np.all([gtdb_df['phylum']==p, 
                                                np.any([gtdb_df['novel']==str(i) for i in range(0,novel_lvl+1)], axis=0)], axis=0)].index.values
        if len(known_genomes) > 0:
            known_tree = tree.extract_tree_with_taxa_labels(labels=known_genomes)
            known_length = known_tree.length()
        else:
            known_length = 0
        cur_lst += [known_length]
    pd_lst.append(cur_lst)
cur_lst = ['Total', tree.length()]
for novel_lvl in [0,1,2]:
    if novel_lvl == 0:
        known_genomes = gtdb_df.loc[gtdb_df['novel']==str(novel_lvl)].index.values
    else:
        known_genomes = gtdb_df.loc[np.any([gtdb_df['novel']==str(i) for i in range(0,novel_lvl+1)], axis=0)].index.values
    known_tree = tree.extract_tree_with_taxa_labels(labels=known_genomes)
    known_length = known_tree.length()
    cur_lst += [known_length]
pd_lst.append(cur_lst)

pd_df = pd.DataFrame(pd_lst, columns=['Phylum', 'Total', 'Known', 'Known+Novel1', 'Known+Novel1,2'])
pd_df = pd_df.set_index('Phylum')
pd_df['Improvement'] = (pd_df['Total']-pd_df['Known'])/pd_df['Total']*100
pd_df['Novel1'] = (pd_df['Known+Novel1']-pd_df['Known'])/pd_df['Total']*100
pd_df['Novel2'] = (pd_df['Known+Novel1,2']-pd_df['Known+Novel1'])/pd_df['Total']*100
pd_df['Novel3'] = (pd_df['Total']-pd_df['Known+Novel1,2'])/pd_df['Total']*100
pd_df['Known_prop'] = pd_df['Known']/pd_df['Total']*100
