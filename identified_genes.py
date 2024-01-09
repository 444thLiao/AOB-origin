from collections import defaultdict
import pandas as pd
from glob import glob
from os.path import *
import os
from subprocess import check_call
from tqdm import tqdm
from api_tools.itol_func import *
from Bio import SeqIO
from collections import Counter

def run_blast(out_dir,in_proteins_dir=None,genome_ids=None,redo=False,
              dbfile = ''):
    
    if not exists(out_dir):
        os.system(f"mkdir -p {out_dir}")
    infile = dbfile #"/mnt/home-backup/thliao/AOB/whole_tree_20200412/rebuild_gene_tree/nxrA/prepared.faa"

    # redo = False
    cmds = []
    faa_list = []
    if in_proteins_dir is not None:
        faa_list = glob(join(in_proteins_dir, '*.faa'))
    elif genome_ids is not None:
        faa_list = [f"/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files/{gid}.faa" for gid in genome_ids]
    else:
        raise IOError
    for faa in tqdm(faa_list):
        name = basename(faa).replace('.faa', '')
        ofile = f"{out_dir}/{name}.tab"
        # cmd = f"blastp -db {infile} -query {faa} -outfmt 6 -max_target_seqs 10000000 -evalue 1e-20 -num_threads 20 > {ofile}"
        cmd = f"diamond blastp --quiet -q {faa} -o {ofile} -d {infile} -p 20 -k 0 -f 6"
        if not exists(ofile) or redo:
            check_call(cmd, shell=1)
    for _ in tqdm(cmds):
        check_call(_, shell=1)
        
def parse_blast(indir=None,
                outdir=None,
                infiles=[],
                locus2gene={},
                seq_file='',
                locus2genome_mapping=None):
    if outdir is not None and not exists(outdir):
        os.makedirs(outdir)
    if not infiles:
        infiles = glob(join(indir, '*.tab'))
        
    def get_genome2genes(sub_df,locus2gene,locus_any):
        
        _sub_df = sub_df.loc[sub_df[1].isin(locus2gene)]
        l2g = {}
        for l,idx_list in _sub_df.groupby(0).groups.items():
            genes = [locus2gene[_] for _ in sub_df.loc[idx_list,1]]
            gene = sorted(genes,key=lambda x:Counter(genes)[x])[-1]
        # genes = [locus2gene[_] for _ in _sub_df[1]]
            l2g[l] = gene
        # l = _sub_df[0]
        # l2g = dict(zip(l,genes))
        # locus2gene_sub.update()
        # sub_df.loc[:,'gene'] = [locus2gene.get(_,'') for _ in sub_df[1]]
        # genes = set(genes)
        genes = list(set(l2g.values()))
        
        if locus2genome_mapping is None:
            genome = "GCA_"+locus_any.split('_')[0].replace('v','.')
        else:
            genome = locus2genome_mapping[locus_any.rpartition('_')[0]]
        return {genome: genes},l2g
    removed_false_positive_ref_faa = ['Y14338.1']
    
    infile = seq_file #"/mnt/home-backup/thliao/AOB/whole_tree_20200412/rebuild_gene_tree/nxrA/prepared.faa"
    curated_gene2length = {_.id: len(str(_.seq))
                        for _ in SeqIO.parse(infile, format='fasta')}
    genome2genes = {}
    locus2gene_focal_phyla = {}

        
    for tab in tqdm(infiles):
        if os.path.getsize(tab) == 0:
            continue
        df = pd.read_csv(tab,sep='\t',header=None)
        df = df.loc[~df[1].isin(removed_false_positive_ref_faa),:]
        if df.shape[0]==0:
            continue
        df.loc[:,'total_length'] = [curated_gene2length[_] for _ in df[1]]
        df.loc[:,'qcov'] = abs(df[7]-df[6])/df['total_length']
        df.loc[:,'gene'] = [locus2gene[_] for _ in df[1]]
        sub_df = df.loc[(df[10]<=1e-20) & (df[2] >=30) & (df['qcov'] >0.4) ]
        sub_df = sub_df.groupby(0).head(3)
        g = df.iloc[0,0]
        _dict,l2g_d = get_genome2genes(sub_df,locus2gene,g)
        genome2genes.update(_dict)
        locus2gene_focal_phyla.update(l2g_d)
        
    # write out the results incase repeat run above codes
    with open(join(outdir,'genome2genes.txt'),'w') as f1:
        for genome,genes in genome2genes.items():
            f1.write('\t'.join([genome] + list(genes)) + '\n')
    with open(join(outdir,'locus2genes.list'),'w') as f1:
        for locus,gene in locus2gene_focal_phyla.items():
            f1.write('\t'.join([locus,gene])+ '\n')
    return genome2genes,locus2gene_focal_phyla

# info2style = dict(cnxrA=dict(color="#3d5afe"),
#                 pnxrA=dict(color="#3d5afe"), 
#                 narG1=dict(color="#f87a0d"), 
#                 narG2=dict(color="#f87a0d"),
#                 )
# text = to_binary_shape({k:[v] for k,v in locus2gene_focal_phyla.items()},
#                     info2style=info2style,
#                     unfilled_other=True)
# with open(join(outdir,'locus2gene.txt'),'w') as f1:
#     f1.write(text)


# # write out annotation
# _info2style = dict(nxrA=dict(color="#3d5afe"), 
#                 narG=dict(color="#f87a0d"), 
#                 both=dict(color="#00796B"))
# #text = to_color_range(_g2g, info2color=info2style)
# text = to_binary_shape(genome2genes,info2style=_info2style,unfilled_other=True)
# with open(join(outdir,'target_genes_binary.txt'),'w') as f1:
#     f1.write(text)


# info2style = dict(cnxrA=dict(color="#3d5afe"),
#                 pnxrA=dict(color="#3d5afe"), 
#                 narG1=dict(color="#f87a0d"), 
#                 narG2=dict(color="#f87a0d"),
#                 )
# text = to_binary_shape(genome2genes_finer,info2style=info2style,unfilled_other=True,info_name='finer nxrA')
# with open(join(outdir,'target_genes_binary_finer.txt'),'w') as f1:
#     f1.write(text)
    

# all_locus = {_:_ for _ in locus2gene_focal_phyla_finer}
# text = to_label(all_locus)
# with open(join(outdir,'reset_names.txt'),'w') as f1:
#     f1.write(text)
        




base_odir = "/mnt/home-backup/thliao/AOB/analysis/update_calibrations"

all_genomes = open('/mnt/home-backup/thliao/AOB/whole_tree_20200412/all_1471.txt').read().strip('\n').split('\n')
base_in_db_dir = '/mnt/home-backup/thliao/AOB/analysis/query_db/fa/'

# make db
target_genes  = ["amoA", "amoB", "amoC", "cycA", "cycB", "hao", "nxrA", "nxrB"]

ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'
ref_df = pd.read_excel(ref_file, index_col=None)
ref_df = ref_df.loc[ref_df.loc[:, 'note'] != 'removed', :]
gb = ref_df.groupby('outgroup/ref for which KO')

ko2gene_name = {'K00370': 'nxrA',
                'K00371': 'nxrB',
                'K10535': 'hao',
                'K10944': 'amoA',
                'K10945': 'amoB',
                'K10946': 'amoC'}
for ko,idx in gb.groups.items():
    sub_df = ref_df.loc[idx,:]
    genes = list(sub_df['AA accession'])
    with open(f'/mnt/home-backup/thliao/AOB/analysis/query_db/fa/{ko2gene_name[ko]}_db.faa','w') as f1:
        for g,seq in zip(list(sub_df['AA accession']),
                         list(sub_df['seq'])):
            f1.write(f'>{g}\n{seq}\n')
locus2gene = dict(zip([_.strip() for _ in ref_df['AA accession']],
                      ref_df['gene name']))
with open('/mnt/home-backup/thliao/AOB/analysis/query_db/fa/locus2gene.list','w') as f1:
    f1.write('\n'.join([f"{k}\t{v}" for k,v in locus2gene.items()]))

for infa in glob('/mnt/home-backup/thliao/AOB/analysis/query_db/fa/*_db.faa'):
    if not exists(infa+'.dmnd'):
        cmd = f"diamond makedb --in {infa} -d {infa}"
        check_call(cmd,shell=1)

# run blast
for infa in sorted(glob('/mnt/home-backup/thliao/AOB/analysis/query_db/fa/*_db.faa')):
    gene = basename(infa).split('_')[0]
    run_blast(out_dir=join(base_odir,'gene_annotations',gene),
              genome_ids=all_genomes,
              dbfile=infa,
              )

        
# parse blast
gene2results = {}
for infa in sorted(glob('/mnt/home-backup/thliao/AOB/analysis/query_db/fa/*_db.faa')):
    gene = basename(infa).split('_')[0]
    genome2genes,locus2gene_focal_phyla = parse_blast(indir=join(base_odir,'gene_annotations',gene),
            seq_file = infa,
            locus2gene=locus2gene,
            outdir=join(base_odir,'gene_annotations',gene+'_o'),
            )
    gene2results[gene] = genome2genes,locus2gene_focal_phyla
    

info2style = {'amoA': {'color': '#ff0000',
                       'info': 'amoA'},
              'amoB': {'color': '#ff0000',
                       'info': 'amoB'},
              'amoC': {'color': '#ff0000',
                       'info': 'amoC'},
              'hao': {'color': '#b68100',
                      'info': 'hao'},
              'nxrA': {'color': '#4b85c1',
                       'info': 'nxrA'},
              'nxrB': {'color': '#4b85c1',
                       'info': 'nxrB'},
              }

genome2genes = defaultdict(list)
for gene,(g2g,_2) in gene2results.items():
    genomes = [k for k,v in g2g.items() if gene in v]
    for g in genomes:
        genome2genes[g].append(gene)
genome2genes = dict(genome2genes)

requested_genes = ['amoA','amoB','amoC','hao', 'nxrA','nxrB']
text = to_binary_shape(genome2genes, 
                       info2style, 
                       manual_v=requested_genes,
                       info_name="nitrification",
                       unfilled_other=True)
with open(join(base_odir,"genome2genes.txt"), "w") as f1:
    f1.write(text)





locus2gene = {}
for of in glob('./gene_annotations/*_o/locus2genes.list'):
    _dict = {row.split('\t')[0]:row.split('\t')[1] 
             for row in open(of).read().strip('\n').split('\n')}
    locus2gene.update(_dict)


homologues = {'amoA':['pmoA',
                      'pxmA',
                      'emoA',
                      'amoA'
                      ],
              'amoB':['pmoB',
                      'pxmB',
                      'emoB',
                      'amoB'
                      ],
              'amoC':['pmoC',
                      'pxmC',
                      'emoC',
                      'amoC'
                      ],
              }


color_scheme_g = {"amo": {"color": "#e53238", "shape": "2"},
                  "bmo": {"color": "#00C853", "shape": "2"},
                  "pxm": {"color": "#FF6F00", "shape": "2"},
                  "pmo": {"color": "#0064d3", "shape": "2"},
                  "emo": {"color": "#D500F9", "shape": "2"},

                  "nxrA": {"color": "#408bd2", "shape": "2"},
                  "unknown": {"color": "#d509e9", "shape": "2"},
                  "narG": {"color": "#f87a0d", "shape": "2"},
                  "dms_tor": {"color": "#f70000", "shape": "2"},
                  }

for name, genes in homologues.items():
    sub_l2gene = {k:v for k,v in locus2gene.items() if v in genes}
    genome2genes = defaultdict(list)
    for l,gene in sub_l2gene.items():
        genome = f"GCA_"+l.split('_')[0].replace('v','.')
        genome2genes[genome].append(gene)
    text = to_binary_shape(genome2genes,
                           {v:color_scheme_g[v[:-1]] for v in set(sub_l2gene.values())},
                           unfilled_other=True,
                           info_name=f'{name} homologues')
    
    with open(f'./gene_annotations/itol_txt/genome2genes_{name}.txt','w') as f1:   
        f1.write(text)
        
all_locus2seq = {}
all_genome = set([f"GCA_{l.split('_')[0].replace('v','.')}" for l in locus2gene])
for genome in tqdm(all_genome):
    faa = f"/home-user/thliao/data/NCBI/modified_data/direct_protein_files/{genome}.faa"
    records = list(SeqIO.parse(faa,'fasta'))
    
    locus = genome.split('_')[-1].replace('.','v')
    all_locus = [l for l,g in locus2gene.items() if l.startswith(locus)]
    
    all_locus2seq.update({r.id:r for r in records if r.id in all_locus})

all_genes = set(locus2gene.values())
for name,genes in homologues.items():
    sub_l2g = {l:g for l,g in locus2gene.items() if g in genes}
    seqs = [all_locus2seq[_] for _ in sub_l2g]
    with open(f'./gene_annotations/genes/{name}_homologues.faa','w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
    