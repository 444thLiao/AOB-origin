from os.path import *
import os
from glob import glob
from ete3 import Tree
from collections import defaultdict
from os.path import exists
from Bio import SeqIO

from pack_up_pipelines.dating.formal.final_params import get_gene_name,strategies_params,genes,gene_names

os.chdir('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces')
# calculating delta-LL
cmds = []
for contree in glob('phy_files/*prot/withEuk/*.treefile') + glob('phy_files/cog25_prot/*.treefile'):
    gene = contree.split('/')[-1].replace('.treefile','')
    gn = get_gene_name(gene)
    if gn != 'P39':continue
    aln = contree.replace('.treefile','.aln')
    genomes = list([_.id for _ in SeqIO.parse(aln,'fasta')])
    t = ''
    st = strategies_params[gn][3]
    st = Tree(st,1)
    st.prune(genomes)
    t+=st.write()+'\n'
    t+=open(contree).read()+'\n'
    tlist = contree.replace('.treefile','.list')
    with open(tlist,'w') as f1:
        f1.write(t)
    odir = dirname(contree)+'/deltaLL'
    if not exists(odir):
        os.makedirs(odir)
    cmd = f"iqtree -nt 10 -m TESTMERGE -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -pre {odir}/{gene} -s {aln} -n 0 -zb 10000 -au -z {tlist}"
    #if not exists(f"{odir}/{gene}.iqtree"):
    cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=4,prefix_name='gt')


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

### batch draw phylogenetic tree
from api_tools.IO_for.itol import get_itoltree
from glob import glob
import pandas as pd
from tqdm import tqdm
from visualization.validated_tips.merged_multipl_pdf import merged_pdfs_with_texts

setname2sorted_genes = defaultdict(list)
for gene,deltall in gene2dl.items():
    setname = get_gene_name(gene)
    setname2sorted_genes[setname].append(gene)
setname2sorted_genes = {k: sorted(v,key=lambda x: gene2dl[x],reverse=True) for k,v in setname2sorted_genes.items()}


_df = pd.read_csv('/home-db/pub/protein_db/CDD/cdd.info.dat/cddid_all.tbl',sep='\t',header=None,index_col=0)
cog25toname = _df.loc[[int(_) for _ in setname2sorted_genes['COG25']],:][1].to_dict()

wrong_genes = ['ENGA','']


idir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb'
for setname,sorted_genes in setname2sorted_genes.items():
    if setname !='G18':continue
    #genes = sorted_genes[:3] + sorted_genes[-2:]
    #genes = sorted_genes[:2] + sorted_genes[-2:]
    genes =  sorted_genes
    for g in genes:
        if setname == 'COG25':
            name = cog25toname[int(g)]
        else:
            name = g
        tpath = glob(f'phy_files/*prot/withEuk/{g}.treefile')
        if not tpath:
            tpath = glob(f'phy_files/cog25_prot/{g}.treefile')
        tpath = tpath[0]
        tre = Tree(tpath)
        #tre.set_outgroup(tre.get_midpoint_outgroup())
        tre.write(outfile='./tmp2.tre')
        #os.system(f"FastRoot.py -i ./tmp2.tre -o ./tmpMV.tre -m MV")
        os.system(f"~/software/mad/mad  ./tmp2.tre; mv ./tmp2.tre.rooted /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb/MAD_rooted_deltaLL/{name}.tre")
        get_itoltree(f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb/MAD_rooted_deltaLL/{name}.tre',
                    outfile=f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb/MAD_rooted_deltaLL2/{name} ({gene2dl[g]}).pdf",
                    anno_files=[f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/ALE_archive/names.txt',
                                f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/ALE_archive/lineage_colorstrip.txt'],
                    keys={'ignore_branch_length':'1',
                        'line_width':"10"})
# Get the largest 3 and bottom 1 genes for each set of genes
for setname,sorted_genes in setname2sorted_genes.items():
    genes = sorted_genes[:2] + sorted_genes[-2:]
    files = []
    for g in genes:
        if setname == 'COG25':
            name = cog25toname[int(g)]
        else:
            name = g
        # the leftmost one has the samller likelihood
        files.append(f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb/MAD_rooted_deltaLL/{name} ({gene2dl[g]}).pdf")
    merged_pdfs_with_texts(files,
                           ofile=f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb/MAD_rooted_deltaLL//{setname}.pdf' ,shape=4)  


# draw a heatmap with plotly
import pandas as pd
df = pd.read_csv('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/cog25_prot/deltaLL/summary.csv',index_col=0)

### draw species tree
idir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb'
for setname in setname2sorted_genes:
    tpath = strategies_params[setname][3]
    tre = Tree(tpath,8)
    tre.write(outfile='./tmp2.tre')
    get_itoltree('./tmp2.tre',
                outfile=f"{idir}/{setname} ST.pdf",
                anno_files=[f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/ALE_archive/names.txt',
                            f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/ALE_archive/lineage_colorstrip.txt'],
                keys={'ignore_branch_length':'1',
                    'line_width':"10"})
merged_pdfs_with_texts(glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb/* ST.pdf'),
                       ofile=f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ipynb/MERGED.pdf' ,shape=1)  

import plotly.graph_objects as go
def get_heatmap(df,ofile):
    fig = go.Figure(data=go.Heatmap(
                   z=df.values,
                   x=df.columns,
                   y=df.index,
                   colorscale='Viridis'))
    fig.update_layout(
        title='delta-LL',
        xaxis_nticks=36)
    fig.write_image(ofile)



###! stepwisely removing genes that

# be careful which part of genes are removed. those with higher/lower values.
import re
from pack_up_pipelines.dating.smaller_genomes.quick_run import gen_pf_list
cmds = []
bdir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise'
for setname,g in zip(gene_names,
                     genes):
    if setname not in ['COG25',
                       'P39',
                       'G18',
                       'M24']:continue
    sortedg = sorted(g,key=lambda x:gene2dl[x])[::-1]  # ascending
    #print(setname,gene2dl[sortedg[0]],gene2dl[sortedg[-1]])
    for i in range(0,40):
        num_removed_g = i*2
        remaining_genes = sortedg[num_removed_g:]
        if len(remaining_genes) <= 3:break
        odir = f"{bdir}/{setname}/r{num_removed_g}"
        if not exists(odir):
            os.makedirs(odir)
        with open(f'{odir}/R{setname}.list','w') as f1:
            f1.write('\n'.join(remaining_genes))
        phy_file = f'{odir}/R{setname}_pf/3pf.phy'
        fodir = dirname(dirname(phy_file))
        if not exists(phy_file):
            gen_pf_list(f'{odir}/R{setname}.list',
                    indir=strategies_params[setname][2] ,
                    pf_range=[3])
        if setname in ['COG25']:
            caltree = strategies_params[setname][1]
            cal_set = re.findall("_([ABP][1-9]+).",caltree)[0]
            _omcmc = f"{fodir}/{cal_set}_MCMC_updated"
            cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {_omcmc} -rg '1 30 1' -sf 30 -c 2 -p 1 "
            if not exists(_omcmc):
                cmds.append(cmd)
        else:
            caltree = strategies_params[setname][4]
            cal_set = re.findall("_([ABP][1-9]+[ABP][1-9]+).",caltree)[0]
            calname = caltree.split('/')[-1].replace('.newick','')
            _omcmc = f"{fodir}/{cal_set}_MCMC_updated"
            cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {_omcmc} -rg '1 30 1' -sf 30 -c 2 -p 1 "
            if not exists(_omcmc):
                cmds.append(cmd)
from bin.multiple_sbatch import sbatch_all
sbatch_all(cmds,thread_per_tasks=1,prefix_name='3pf',fixed_cluster='others')




cmds = []
bdir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise'
for setname,g in zip(gene_names,genes):
    sortedg = sorted(g,key=lambda x:gene2dl[x])[::-1]  # ascending
    # for i in [0,1,2,3,4]:
    i = 2
    num_removed_g = i*2
    remaining_genes = sortedg[num_removed_g:]
    odir = f"{bdir}/{setname}/r{num_removed_g}"
    for p in [1,2,3,4,5]:
        phy_file = f'{odir}/R{setname}_pf/{p}pf.phy'
        if not exists(phy_file):
            gen_pf_list(f'{odir}/R{setname}.list',
                indir=strategies_params[setname][2] ,
                pf_range=[p])
    caltree = strategies_params[setname][1]
    f_odir = f"{bdir}/mcmc_o/B15_R{setname}.1pf.r{num_removed_g}"
    cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {f_odir} -rg '1 30 1' -sf 30 -c 2 -p 1 "
    if not exists(f_odir):
        cmds.append(cmd)
    caltree = strategies_params[setname][4]
    if not caltree:
        continue
    calname = caltree.split('/')[-1].replace('.newick','')
    f_odir = f"{bdir}/mcmc_o/{calname}_R{setname}.1pf.r{num_removed_g}"
    cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {f_odir} -rg '1 30 1' -sf 30 -c 2 -p 1 "
    if not exists(f_odir):
        cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=1,prefix_name='1pf')



cmd = f"python3 ~/script/evol_tk/dating_workflow/figtree2itol.py -i2 /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise/P39/r24/B5P3_MCMC_updated/mcmc_for/FigTree.tre -i /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/tree/bac120_Plastid_revised.newick -o /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise/P39/r24/B5P3_MCMC_updated/mcmc_for/fig3.newick -f 8"

