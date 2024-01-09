import os
os.chdir("/mnt/home-backup/thliao/AOB/analysis/20210713_repeat")
import pandas as pd 
from tqdm import tqdm
from glob import glob
from ete3 import Tree
from api_tools import to_binary_shape

anno_f = "/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/near_genomes_kegg_binary.tab"
anno_df = pd.read_csv(anno_f,sep='\t',index_col=0)
def simp(gid):
    return gid.split('_')[-1].split('.')[0]
#! get identified XMO genomes from the contanated gene tree
# sampling analyzed genomes
import re
gene_tre = Tree("/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/xmoCBA_remove_archaea/xmoCBA_renamed.newick",3)
def get_g(nodes):
    _tmp = []
    for n in [n for n in gene_tre.traverse() if n.name in nodes]:
        _tmp.extend( [re.findall("GC[FA]_\d*\.\d",l)[0] 
                                 for l in n.get_leaf_names()
                                 if l.startswith('GC')])
    return _tmp
suffix2full = {_.split('_')[-1].split('.')[0]:_ for _ in dm.index}
beta_amo = get_g(['I141_S100'])
gamma_amo = get_g(['I53_S100'])
nitro_amo = get_g(['I142_S100'])
pmo_lineage = get_g(['I52_S100']) # gamma
pxm_lineage = get_g(["I82_S100"])
pmo_lineage2 = get_g(['I38_S100'])
from sklearn.neighbors import NearestNeighbors
from collections import Counter
dm = pd.read_csv('./ref_phy/mash_dist/pairwise_merged_db.tab',sep='\t',index_col=0)
neigh = NearestNeighbors(metric='precomputed')
neigh.fit(dm)
def get_near_genomes(seeds,threshold = 300):
    seeds = [_ for _ in seeds if _ in dm.index ]
    fold = threshold//len(seeds)
    N = neigh.kneighbors(dm.loc[seeds,:],fold*len(seeds),return_distance = False)
    v = sorted(Counter(N.flatten()).items(),key=lambda x:x[1])[::-1][:fold*len(seeds)]
    v = [dm.columns[_[0]] for _ in v]
    return v
near_genomes = []
for l in [beta_amo,
            gamma_amo,
            nitro_amo,
            pmo_lineage,
            pxm_lineage,
            pmo_lineage2]:
    tmp = get_near_genomes(l)
    near_genomes += tmp
near_genomes = list(near_genomes)

text = to_binary_shape({k:['retained'] for k in near_genomes},
                       info2style={'retained':{'color':'#ff0000',
                                               'shape':"5"}},
                       unfilled_other=True,
                       info_name='neighbouring genomes'
                       )
with open('./comparisons/near_genomes.txt', 'w') as f1:
    f1.write(text)