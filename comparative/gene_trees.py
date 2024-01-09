import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
anno_f = "/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/1444genomes_kegg_details.tab"
details_df = pd.read_csv(anno_f,sep='\t',index_col=0)

locus2ko = {locus: ko
            for ko, _d in details_df.to_dict().items()
            for genome, locus_list in _d.items()
            for locus in str(locus_list).split(',')}

l2seq = {}
for gid in tqdm(details_df.index):
    faa = f"/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files/{gid}.faa"
    _l2seq = {_.id:_ for _ in SeqIO.parse(faa,'fasta') if _.id in locus2ko}
    l2seq.update(_l2seq)

c = pd.read_excel('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/comparisons/kegg/AOB_signatures.xlsx')
target_ko = c.loc[c['Gain/Loss']=='Gain'].sort_values('corrected p-value')['KO'][:30]
for ko in target_ko:
    if not ko.startswith('K'):continue
    seqs = [l2seq[_] for _,k in locus2ko.items() if k==ko]
    odir = f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_trees/{ko}.faa'
    with open(odir,'w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
        
from glob import glob        
cmds = []
for faa in glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_trees/*.faa'):
    aln = faa.replace('.faa','.aln')
    tre = faa.replace('.faa','.newick')
    cmd = f"mafft --auto --anysymbol {faa} > {aln}; FastTree {aln} > {tre}"
    cmds.append(cmd)
from bin.multiple_sbatch import sbatch_all 
sbatch_all(cmds,thread_per_tasks=5,batch_size=1,prefix_name='gt')

p = {g.split('_')[1].replace('.','v'):g for g in details_df.index}
from ete3 import Tree
for tpath in glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_trees/*.newick'):
    tre = Tree(tpath)
    for l in tre.get_leaves():
        l.name = p[l.name.split('_')[0]]
    tre.write(outfile=tpath.replace('.newick','.renamed.nw'))
        
    