"""
This script include two kinds of genome sampling
First for genome sampling for the reference genome sets (pre-download sampling)
1. sample from summary files
2. build phylogenomic tree
3. visulization

Seoncd for genome sampling for dating analysis
1. select seeds based on three criterias
    a. maximize kmer based similiarity
    b. completeness
2. select backbone genomes
    same as above

Third
The last part would generate a species tree for the reconciliation
"""

from api_tools import *
from bin import ngd
import pandas as pd
from visualization.sunburst_taxon import *
from bin.ngd import *
from ete3 import Tree
from collections import defaultdict
import pandas as pd
from os.path import *
import os
from subprocess import check_call
from tqdm import tqdm
from api_tools.itol_func import *
from collections import Counter
from api_tools.IO_for import read_summary
info_df = pd.read_excel('./xmoA_containing_genomes.xlsx', index_col=0)
sub_info_df = info_df.loc[info_df['superkingdom'] != 'Archaea', :]

tax_tab = "/home-user/thliao/.cache/ncbi-genome-download/taxonomy.tab"
tax_df = pd.read_csv(tax_tab, sep='\t', index_col=0)
tax_df.index = [_.split('_')[-1] for _ in tax_df.index]
genome2tax = tax_df.to_dict(orient='index')

m1 = read_summary('/home-user/thliao/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt')
sub_tax_df = tax_df.reindex([_.split('_')[-1].split('.')[0] for _ in m1.index])
# focal phyla except APF (Proteobacteria, Firmuctes, Actinobacteria)

# sub_tax_df = tax_df.reindex([_.split('_')[-1].split('.')[0] for _ in info_df.index])
# sub_tax_df.groupby('phylum').size().sort_values()
# Acidobacteria                       1
# Nanoarchaeota                       1
# Armatimonadetes                     1
# Candidatus Binatota                 1
# Candidatus Dadabacteria             1
# Candidatus Korarchaeota             1
# Candidatus Marinimicrobia           1
# candidate division CSSED10-310      1
# Candidatus Rokubacteria             2
# Chloroflexi                         4
# Crenarchaeota                       4
# candidate division NC10             6
# Verrucomicrobia                    18
# Actinobacteria                     23
# Nitrospirae                        60
# Thaumarchaeota                    246
# Proteobacteria                    510

n_gids = sub_tax_df.index[sub_tax_df['class'].isin(
    ['Betaproteobacteria', 'Gammaproteobacteria'])]

analyzed_genomes = []
analyzed_genomes.extend(sub_tax_df.loc[sub_tax_df['phylum'] == 'Nitrospirae'])  # 1083
print(len(analyzed_genomes))
analyzed_genomes.extend(sub_tax_df.loc[sub_tax_df['phylum'] == 'candidate division NC10'])  # 72
print(len(analyzed_genomes))
# ! First genome sampling
target_class = ['Betaproteobacteria',
                'Alphaproteobacteria',
                'Gammaproteobacteria']
def statistic_number(tax, threshold=30):
    """
    The threshold indicates
    if the number of genomes within this genus is larger than 50
    """
    sub_df = tax_df.loc[tax_df['class'] == tax]
    sub_df = sub_df.reindex([_ for _ in sub_df.index if _.startswith('GCA')])
    # level2names = {l: sub_df[l].unique()
    #                for l in sub_df.columns}
    level2name2counts = {l: sub_df.groupby(l).size().to_dict()
                         for l in sub_df.columns}
    # if the genus level contain over 50 genomes, it will collapsed with most completed one
    genus_over_threshold = [n
                            for n, c in level2name2counts["genus"].items()
                            if c >= threshold]
    remaining_ids = sub_df.loc[~sub_df['genus'].isin(
        genus_over_threshold)].index
    genus2repr = {}
    for g in genus_over_threshold:
        sub_df2 = sub_df.loc[sub_df['genus'] == g]
        sub_m_df = m1.reindex(sub_df2.index).sort_values('assembly_level')
        gids = list(sub_m_df.index[:2])
        genus2repr[g] = gids
    final_ids = list(remaining_ids) + \
        list(set([_ for v in genus2repr.values() for _ in v]))
    print(
        f"{tax}: original number is {sub_df.shape[0]}, remaining numer is {len(final_ids)}")
    return final_ids


p2ids = {}
all_ids = []
for p in target_class:
    _ids = statistic_number(p)
    _ids = set(_ids)
    p2ids[p] = [_ for _ in m1.index if _.split('.')[0] in _ids]
    all_ids.extend(p2ids[p])
# Betaproteobacteria: original number is 22258, remaining numer is 2544
# Alphaproteobacteria: original number is 27204, remaining numer is 6334
# Gammaproteobacteria: original number is 769332, remaining numer is 5879

# refine with prior knowledge of AOB
# Nitrococcaceae family (Gammaproteobacteria) and the Nitrosomonadaceae family (Betaproteobacteria)
gamma_well_known = sub_tax_df.index[(sub_tax_df['class'] == 'Gammaproteobacteria') & (sub_tax_df['family'] == 'Chromatiaceae')]
beta_well_known = sub_tax_df.index[(sub_tax_df['class'] == 'Betaproteobacteria') & (sub_tax_df['family'] == 'Nitrosomonadaceae')]
print(len(gamma_well_known),len(beta_well_known))
# 246 231

# build phylogeny in each phylum
tax2gids = {}
tax2gids['Nitrospirae'] = sub_tax_df.index[sub_tax_df['phylum'] == 'Nitrospirae']  
tax2gids['NC10'] = sub_tax_df.index[sub_tax_df['phylum'] == 'candidate division NC10']  
tax2gids.update(p2ids)

def build_phy(tax, gids):
    odir = f"./ref_phy/{tax}"
    if not exists(odir):
        os.system(f"mkdir -p {odir}")
    with open(f'{odir}/gids.list', 'w') as f1:
        f1.write('\n'.join(gids))
    bac120_odir = f"./ref_phy/{tax}/bac120_anno/"
    cmd = f"""python3 ~/script/evol_tk/dating_workflow/step_script/extract_bac120.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/bac120_annotate -o {bac120_odir} -evalue 1e-50 -gl {odir}/gids.list;
python3 /home-user/thliao/script/evol_tk/dating_workflow/step_script/get_genome_list_with_percentage.py -i {bac120_odir} -o {odir}/over20p_gids.list -s faa -num_p 20 ;
python3 ~/bin/batch_run/batch_mafft.py -i {bac120_odir} -s faa -o {bac120_odir}/aln -gl {odir}/over20p_gids.list -m mafft -f """
    check_call(cmd, shell=1)
    # cmd = f"""python3 ~/bin/batch_run/batch_any.py -i  {bac120_odir}/aln -s aln -o {bac120_odir}/refined -ns clipkit -cmd 'clipkit -m smart-gap  {infile} -o {ofile} -l ' -f     """
    # check_call(cmd, shell=1)
    # cmd = f"""python3 /home-user/thliao/script/evol_tk/dating_workflow/bin/concat_aln.py -i {bac120_odir}/refined -o {odir}/trees/concat/bac120_over20p.clipkit -s clipkit -gl {odir}/over20p_gids.list -ct partition -no_graph -simple    """
    # check_call(cmd, shell=1)

    check_call(cmd, shell=1)
    cmd = f"""python3 ~/bin/batch_run/batch_trimal.py -i {bac120_odir}/aln -o {bac120_odir}/refined -ns trimal;
    """
    check_call(cmd, shell=1)
    cmd = f"""python3 /home-user/thliao/script/evol_tk/dating_workflow/bin/concat_aln.py -i {bac120_odir}/refined -o {odir}/trees/concat/bac120_over20p.trimal -s trimal -gl {odir}/over20p_gids.list -ct partition -no_graph -simple    """
    check_call(cmd, shell=1)
    if not exists(f"{odir}/trees/iqtree"):
        os.makedirs(f"{odir}/trees/iqtree")
    if not exists(f"{odir}/trees/fasttree"):
        os.makedirs(f"{odir}/trees/fasttree")

    cmds = [#f"iqtree -nt 30 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/trees/iqtree/bac120_over20p_clipkit -s {odir}/trees/concat/bac120_over20p.clipkit",
            f"iqtree -nt 30 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/trees/iqtree/bac120_over20p_trimal -s {odir}/trees/concat/bac120_over20p.trimal",
            #f"iqtree -nt 30 -m TESTONLY -madd LG+C20+G,LG+C30+G,LG+C40+G,LG+C50+G,LG+C60+G -mset WAG,LG,JTT -mrate E,I,G,I+G  -wbtl -bb 1000 -redo -pre {odir}/trees/iqtree/bac120_over20p_test -s {odir}/trees/concat/bac120_over20p.trimal",
            ]
    return cmds
# iqtree -s /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/bac120.trimal -m LG+C60+F+G -ft /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/tree/bac120.contree
# iqtree -nt 30 -m TESTONLY -madd LG+C20+G,LG+C30+G,LG+C40+G,LG+C50+G,LG+C60+G -mset WAG,LG,JTT -mrate E,I,G,I+G  -wbtl -bb 1000 -redo -pre /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/tree/bac120_Complex -s /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/bac120.trimal

# Nitrospirae 1083
# NC10 72
# Betaproteobacteria 2518
# Alphaproteobacteria 6300
# Gammaproteobacteria 5860

cmds = []
for tax in ['Nitrospirae',
            'Betaproteobacteria',
            'Gammaproteobacteria']:
    gids = open(
        "/home-user/thliao/data/cyano/ref_tree_new/used_genomes.list").read().split('\n') + tax2gids[tax]
    gids = list(set([_ for _ in gids if _]))
    cmds.append(build_phy(tax, gids))

# !following annotation
colormap = {'Nitrospirae': '#edc31d',
            'Betaproteobacteria': '#956cb4',
            'Gammaproteobacteria': '#78fce0',
            'Cyanobacteria': '#F44336',
            "Alphaproteobacteria": "#8c613c"}
gids = [_ for _ in open('./ref_phy/all/gids.list').read().strip().split('\n') ]
gids += ['GCA_000014865.1',
 'GCA_002109495.1',
 'GCA_003015145.1',
 'GCA_008189685.1',
 'GCA_000469665.2',
 'GCA_000264455.2',
 'GCA_900107585.1',
 'GCA_006385135.1',
 'GCA_002924445.1',
 'GCA_007197755.1']
# ten alphaproteobacteria
gid2l = {}
for gi in gids:
    if genome2tax[gi.split('.')[0].split('_')[-1]]['phylum'] == 'Nitrospirae':
        gid2l[gi] = 'Nitrospirae'
    elif genome2tax[gi.split('.')[0].split('_')[-1]]['phylum'] == 'Proteobacteria':
        gid2l[gi] = genome2tax[gi.split('.')[0].split('_')[-1]]['class']

gid2l.update({k: 'Cyanobacteria'
              for k in open("/home-user/thliao/data/cyano/ref_tree_new/used_genomes.list").read().split('\n')})
text = to_color_strip(gid2l, colormap, info_name='lineage names')
with open('./ref_phy/all/lineage_colorstrip.txt', 'w') as f1:
    f1.write(text)

import re
gene_tre = Tree("/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB/NoArc/xmoCBA_renamed.newick",3)
genome_in_amolineage = []
for n in [n for n in gene_tre.traverse() if n.name in ['I112_S100','I53_S100']]:
    genome_in_amolineage.extend([re.findall("GC[FA]_\d*\.\d",l)[0]
                                 for l in n.get_leaf_names()
                                 if l.startswith('GC')])
gene2locus_id_te = open(
    "/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/related_faa/annotations/gene2locus_id").read().split('\n')
genome2genes = defaultdict(set)
for row in gene2locus_id_te:
    if row:
        for prefix in ['GCA_','GCF_']:
            genome = prefix+row.split('\t')[1].split('_')[0].replace('v', '.')
            gene = row.split('\t')[0]
            if gene.startswith('amo') and genome in genome_in_amolineage:
                genome2genes[genome].add(gene)
            elif not gene.startswith('amo'):
                genome2genes[genome].add(gene)

info2style = {'amoA': {'color': '#ff0000', 'info': 'amoA','shape':'2'},
              'amoB': {'color': '#ff0000', 'info': 'amoB','shape':'2'},
              'amoC': {'color': '#ff0000', 'info': 'amoC','shape':'2'},
              'cycA': {'color': '#00ff00', 'info': 'cycA','shape':'2'},
              'cycB': {'color': '#00ff00', 'info': 'cycB','shape':'2'},
              'hao': {'color': '#f87a0d', 'info': 'hao','shape':'2'},
              'nxrB': {'color': '#0000ff', 'info': 'nxrB','shape':'2'},
              'nxrA': {'color': '#0000ff', 'info': 'nxrA','shape':'2'},}
text = to_binary_shape(genome2genes,
                       info2style=info2style,
                       unfilled_other=True,
                       manual_v=['amoA', 'amoB', 'amoC', 'hao',
                                 'cycA', 'cycB', 'nxrA', 'nxrB'],
                       # extra_replace=
                       )
with open('./ref_phy/all/genes_binary.txt', 'w') as f1:
    f1.write(text)


#! auto collapse (collaps internal nodes which don't contain any xmoA containing genomes)
home = os.getenv("HOME")
indir = f'{home}/template_txt/'
collapse_template = join(indir, "collapse.txt")
collapse_template_txt = open(collapse_template).read()
def get_collapsed_n(tre, necessary_l):
    collapsed_n = []
    for n in tqdm(tre.traverse()):
        if n.is_root():
            continue
        if (not n.is_leaf()) and (n not in collapsed_n):
            leaves = n.get_leaf_names()
            if set(leaves).intersection(set(necessary_l)):
                continue
            else:
                collapsed_n.extend(list(n.traverse()))
    set_names = set([_.name for _ in collapsed_n])
    purfied_collapsed_n = collapsed_n[::]
    for n in collapsed_n:
        if n.up.name in set_names:
            purfied_collapsed_n.remove(n)
    return purfied_collapsed_n
tre = Tree("./ref_phy/all/trees/fasttree/bac120_over20p_renamed.newick", 3)
necessary_l = set(genome2genes)
collapsed_n = get_collapsed_n(tre, necessary_l)
collapsed_leaves = set([_.name
                        for v in collapsed_n
                        for _ in v.traverse()
                        if _.is_leaf()])
with open("./ref_phy/all/collapsed.txt", 'w') as f1:
    f1.write(collapse_template_txt + '\n' +
             '\n'.join([_.name for _ in collapsed_n]))

cn2l = {}
for cn in collapsed_n:
    sub_leaves = cn.get_leaf_names()
    l = [gid2l.get(_,'') for _ in sub_leaves]
    l_count = Counter(l)
    cn2l[cn.name] = list(sorted(l_count,
                                key=lambda x:l_count[x]))[-1]  # most common one
gid2l.update(cn2l)
text = to_color_strip(gid2l, colormap, info_name='lineage names')
with open('./ref_phy/all/lineage_with_collapsed_nodes_colorstrip.txt', 'w') as f1:
    f1.write(text)


cn2new_name = {}
for cn in collapsed_n:
    num_g = len(cn.get_leaf_names())
    cn2new_name[cn.name] = f"Condensed Node ({num_g})"
n2n = {k:k for k in tre.get_leaf_names()}
cn2new_name.update(n2n)
text = to_label(cn2new_name)
with open("./ref_phy/all/collapsed_id2name.txt", 'w') as f1:
    f1.write(text)


# ! kmers based downsampling
gids = [_ for _ in open(
    "./ref_phy/all/over20p_gids.list").read().split('\n') if _]
info_df = pd.read_excel('xmoA_containing_genomes.xlsx',index_col=0)
gids += list(info_df.index)
gids = list(set(gids))
gids.remove('GCF_000764305.1')

collect_faas = []
for name in set(gids):
    faa = f"/mnt/maple/thliao/data/NCBI/modified_data/prokka_o/{name}/{name}.faa"
    if not exists(faa):
        print(name)
    collect_faas.append(faa)
with open('./ref_phy/quality/micomplete_o/input.tab', 'w') as f1:
    f1.write('\n'.join([faa+'\tfaa' for faa in collect_faas]))
cmd = "miComplete ./ref_phy/quality/micomplete_o/input.tab --weights Bact105 --hmms Bact105  --threads 30 > ./ref_phy/quality/micomplete_o/output.completeness"
check_call(cmd, shell=1)


cmds = []
for fna in tqdm(collect_fnas):
    name = basename(fna).replace('.fna', '')
    cmd = f"mash sketch -o ./ref_phy/mash_dist/db/{name} {fna} -s 10000"
    if not exists(f"./ref_phy/mash_dist/db/{name}.msh"):
        cmds.append(cmd)
check_call("mash paste ./ref_phy/mash_dist/merged_db ./ref_phy/mash_dist/db/*.msh", shell=1)
check_call("mash dist -p40 ./ref_phy/mash_dist/merged_db.msh ./ref_phy/mash_dist/merged_db.msh > ./ref_phy/mash_dist/pairwise_merged_db.dist ", shell=1)
g2g_mash_v = defaultdict(dict)
for row in tqdm(open("./ref_phy/mash_dist/pairwise_merged_db.dist"),total=54937456):
    s1, s2, mash_v, p, hashed = row.split('\t')
    g1 = s1.split('/')[-1].replace('.fna', '')
    g2 = s2.split('/')[-1].replace('.fna', '')
    g2g_mash_v[g1][g2] = float(mash_v)
    g2g_mash_v[g2][g1] = float(mash_v)
df = pd.DataFrame.from_dict(g2g_mash_v)
gids = list(df.index)
df = df.reindex(gids, columns=gids)
df.to_csv(f"./ref_phy/mash_dist/pairwise_merged_db.tab", sep='\t', index=1)



dm_df = pd.read_csv(f"./ref_phy/mash_dist/pairwise_merged_db.tab", sep='\t',index_col=0)
dm = DistanceMatrix.read(f"./ref_phy/mash_dist/pairwise_merged_db.tab")
tre = Tree("./ref_phy/all/trees/fasttree/bac120_over20p_renamed.newick", 3)
## prototype selection
import re
gene_tre = Tree("/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/xmoCBA_remove_archaea/xmoCBA_renamed.newick",3)
def get_g(nodes):
    _tmp = []
    for n in [n for n in gene_tre.traverse() if n.name in nodes]:
        _tmp.extend( [re.findall("GC[FA]_\d*\.\d",l)[0]
                                 for l in n.get_leaf_names()
                                 if l.startswith('GC')])
    return _tmp

beta_amo = get_g(['I141_S100'])
gamma_amo = get_g(['I53_S100'])
nitro_amo = get_g(['I142_S100'])
pmo_lineage = get_g(['I52_S100']) # gamma
pxm_lineage = get_g(["I82_S100"])
pmo_lineage2 = get_g(['I38_S100'])
from ..third_party.prototypeSelection import prototype_selection_destructive_maxdist

seeds = []
for l in [beta_amo,
            gamma_amo,
            nitro_amo,
            pmo_lineage,
            pxm_lineage,
            pmo_lineage2]:
    seeds += list(set([_ for _ in l if _ in dm.ids]))
prototypes = prototype_selection_destructive_maxdist(dm, 500,seedset=seeds)
gsets = prototypes+ open(
        "/home-user/thliao/data/cyano/ref_tree_new/used_genomes.list").read().split('\n')
gsets = list(set([_ for _ in gsets if _]))

gsets += ['GCA_000014865.1',
 'GCA_002109495.1',
 'GCA_003015145.1',
 'GCA_008189685.1',
 'GCA_000469665.2',
 'GCA_000264455.2',
 'GCA_900107585.1',
 'GCA_006385135.1',
 'GCA_002924445.1',
 'GCA_007197755.1']
gsets.remove('GCA_002682985.1')
# ten alphaproteobacteria
cmds = build_phy('neighbouring_genomes',gsets)
# sampled2 contains ten more alphaprotoebacteria than sampled_1

text = to_binary_shape({k:['retained'] for k in gsets},
                       info2style={'retained':{'color':'#ff0000',
                                               'shape':"5"}},
                       unfilled_other=True,
                       info_name='retained genomes')
with open('./ref_phy/neighbouring_genomes/selected_genomes.txt', 'w') as f1:
    f1.write(text)


# for reconciliation
operon_df = pd.read_csv('./gene_annotations/xmoCBA/operon_info.tab',sep='\t',index_col=0)
cmds = build_phy('recon_species',list(set(operon_df.index)))
# sampled2 contains ten more alphaprotoebacteria than sampled_1

text = to_binary_shape({k:['retained'] for k in prototypes},
                       info2style={'retained':{'color':'#ff0000',
                                               'shape':"5"}},
                       unfilled_other=True,
                       info_name='retained genomes'
                       )
with open('./ref_phy/all/selected_genomes.txt', 'w') as f1:
    f1.write(text)