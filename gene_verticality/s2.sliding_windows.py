from os.path import *
import os
from glob import glob
from ete3 import Tree
from collections import defaultdict
from os.path import exists
from Bio import SeqIO

from pack_up_pipelines.dating.formal.final_params import get_gene_name,strategies_params,genes,gene_names

os.chdir('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces')


gene2num = {}
gene2dl = {}
for iqtree in glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/*prot/withEuk/deltaLL/*.iqtree') + glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/cog25_prot/deltaLL/*.iqtree'):
    gene = iqtree.split('/')[-1].replace('.iqtree','')
    rows = open(iqtree).read().strip().split('\n')
    idx = [idx for idx,v in enumerate(rows) if 'deltaL  bp-RELL' in v][0]
    r1, r2 = rows[idx+2],rows[idx+3]
    r1 = [_ for _ in r1.strip().split(' ') if _]
    r2 = [_ for _ in r2.strip().split(' ') if _]
    if r2[2]=='0':
        gene2dl[gene] = float(r1[2])
    else:
        gene2dl[gene] = float(r2[2])
    gene2num[gene] = len(Tree(iqtree.replace('.iqtree','.treefile')).get_leaf_names())
# smallest deltaL indicate fewer topological inconsistence

setname2sorted_genes = defaultdict(list)
for gene,deltall in gene2dl.items():
    setname = get_gene_name(gene)
    setname2sorted_genes[setname].append(gene)
setname2sorted_genes = {k: sorted(v,key=lambda x: gene2dl[x],reverse=True) for k,v in setname2sorted_genes.items()}

import pandas as pd
_df = pd.read_csv('/home-db/pub/protein_db/CDD/cdd.info.dat/cddid_all.tbl',sep='\t',header=None,index_col=0)
cog25toname = _df.loc[[int(_) for _ in setname2sorted_genes['COG25']],:][1].to_dict()


# be careful which part of genes are removed. those with higher/lower values.
gene2tax = {}
for setname,g in tqdm(zip(gene_names,genes)):
    for _g in g:
        _indir=strategies_params[setname][2]
        gene2tax[_g] = [_.id for _ in SeqIO.parse(join(_indir,_g+'.aln'),'fasta')]
            
import re
import numpy as np
from pack_up_pipelines.dating.smaller_genomes.quick_run import gen_pf_list
from tqdm import tqdm
cmds = []
bdir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/sliding_windows'
len_intervals = 10
nstep = 1
for setname,g in tqdm(zip(gene_names,genes)):
    if setname not in [#'COG25',
                       #'P39',
                       'G18',
                       #'M24'
                       ]:continue
    sortedg = sorted(g,key=lambda x:gene2dl[x])[::-1]  # descending
    base = 0
    while 1:
        used_genes = sortedg[base:base+len_intervals]
        if len(used_genes) != len_intervals:
            break
        mean_deltaLL = np.mean([gene2dl[gene] for gene in used_genes])
        odir = f"{bdir}/{setname}_{len_intervals}/move{base}_{mean_deltaLL:.2f}"
        if not exists(odir):
            os.makedirs(odir)
        with open(f'{odir}/R{setname}.list','w') as f1:
            f1.write('\n'.join(used_genes))
        full_tax = [_.split('\t')[0] for _ in open(strategies_params[setname][0]).read().strip().split('\n')]
        if len(set([_ 
                    for k,v in gene2tax.items() 
                    for _ in v if k in used_genes]))<len(full_tax):
            base += nstep
            continue
        ### change the partition scheme here
        phy_file = f'{odir}/R{setname}_pf/2pf.phy'
        fodir = dirname(dirname(phy_file))
        if not exists(phy_file):
            gen_pf_list(f'{odir}/R{setname}.list',
                    indir=strategies_params[setname][2] ,
                    pf_range=[2])
        ### change the calibration it used here
        if setname == 'COG25':
            caltree = strategies_params[setname][1]
            cal_set = re.findall("_([ABP][1-9]+).",caltree)[0]
        else:
            caltree = strategies_params[setname][4]
            cal_set = re.findall("_([ABP][1-9]+[ABP][1-9]+).",caltree)[0]
        _omcmc = f"{fodir}/{cal_set}_MCMC_updated"
        cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {_omcmc} -rg '1 30 1' -sf 30 -c 2 -p 1 "
        if not exists(_omcmc):
            cmds.append(cmd)
        base += nstep
from bin.multiple_sbatch import sbatch_all
sbatch_all(cmds,thread_per_tasks=1,prefix_name='sw',fixed_cluster='others')






