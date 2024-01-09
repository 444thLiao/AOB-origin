
"""
Using the latest genome set updated at 20220818

further retrieve the main taxonomic groups
gids = tax2gids['Nitrospirae']+gamma_well_known+beta_well_known+tax2gids['NC10']

1. Chromatiaceae family (Gammaproteobacteria) 
2. the Nitrosomonadaceae family (Betaproteobacteria)
3. NC10
4. Nitrospirae

{'Betaproteobacteria': 88,
'Nitrospirae': 73,
'NC10': 5,
'Gammaproteobacteria': 2}
Only add above genomes, no one is belong to gamma-AOB
suggesting the rareness of gamma-AOB
         
"""
from bin.multiple_sbatch import sbatch_all
from api_tools import *
from bin import ngd
import pandas as pd
from visualization.sunburst_taxon import *
from bin.ngd import *
from collections import defaultdict
import pandas as pd
from os.path import *
import os
from tqdm import tqdm
from api_tools.itol_func import *
from api_tools.IO_for import read_summary
from collections import Counter

info_df = pd.read_excel('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/xmoA_containing_genomes.xlsx', index_col=0)
tax_tab = "/home-user/thliao/.cache/ncbi-genome-download/taxonomy.tab"
# metadata = "/home-user/thliao/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt"
# metadata2 = "/home-user/thliao/.cache/ncbi-genome-download/genbank_archaea_assembly_summary.txt"
# m2 = read_summary(metadata2)
# m1 = read_summary(metadata)

tax_df = pd.read_csv(tax_tab, sep='\t', index_col=0)
genome2tax = tax_df.to_dict(orient='index')
for _,v in list(genome2tax.items()):
    if _.replace('GCF','GCA') not in genome2tax:
        genome2tax[_.replace('GCF','GCA')] = v
# focal phyla except APF (Proteobacteria, Firmuctes, Actinobacteria)
sub_info_df = info_df.loc[info_df['superkingdom'] != 'Archaea', :]

# sub_info_df.groupby('phylum').size().sort_values()
# Acidobacteria                       1
# Arthropoda                          1
# Candidatus Binatota                 1
# Candidatus Dadabacteria             1
# Candidatus Marinimicrobia           1
# candidate division CSSED10-310      1
# Candidatus Rokubacteria             2
# Chloroflexi                         4
# candidate division NC10             7
# Verrucomicrobia                    35
# Actinobacteria                     41
# Nitrospirae                        63
# Proteobacteria                    698
chosen_phyla = list(sub_info_df['phylum'].unique())
chosen_phyla = [_ for _ in chosen_phyla if str(_) != 'nan']
genome_id = set([genome
                 for genome, _d in genome2tax.items()
                 if _d['phylum'] in chosen_phyla and genome.startswith('GC')])
print(len(set([_.split('_')[-1] for _ in genome_id])))
# 962808  --> 1123417


Counter([genome2tax['GCA_' + _]['phylum']
         for _ in set([_.split('_')[-1] for _ in genome_id])])
# total number of genomes within related phyla
# Counter({'Proteobacteria': 1059271,
#          'Actinobacteria': 46103,
#          'Chloroflexi': 5927,
#          'Verrucomicrobia': 5848,
#          'Acidobacteria': 3475,
#          'Nitrospirae': 1507,
#          'Candidatus Marinimicrobia': 897,
#          'Candidatus Dadabacteria': 102,
#          'Candidatus Binatota': 46,
#          'candidate division NC10': 86,
#          'Candidatus Rokubacteria': 149,
#          'candidate division CSSED10-310': 6})

n_gids = sub_info_df.index[sub_info_df['class'].isin(
    ['Betaproteobacteria', 'Gammaproteobacteria'])]

# missing_gids = []
# for _ in m1.index :
#     if _.split('.')[0] not in genome2tax:
#         missing_gids.append(_)
# print(len(missing_gids))
analyzed_genomes = []
analyzed_genomes.extend([_ for _ in m1.index 
                         if genome2tax.get(_.split('.')[0].replace('GCF','GCA'),{}).get('phylum','') == 'Nitrospirae'])  # 1491
print(len(analyzed_genomes))
analyzed_genomes.extend([_ for _ in m1.index 
                         if genome2tax.get(_.split('.')[0].replace('GCF','GCA'),{}).get('phylum','') == 'candidate division NC10'])  # 83
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
# Betaproteobacteria: original number is 26502, remaining numer is 4772
# Alphaproteobacteria: original number is 34823, remaining numer is 9677
# Gammaproteobacteria: original number is 896467, remaining numer is 9148

# refine with prior knowledge of AOB
# Nitrococcaceae family (Gammaproteobacteria) and the Nitrosomonadaceae family (Betaproteobacteria)

# 
m1 = m1.loc[[_ for _ in m1.index if _.split('.')[0] in genome2tax]]
gamma_well_known = [_ 
                    for _ in m1.index 
                    if (genome2tax[_.split('.')[0]]['class'] == 'Gammaproteobacteria') and (genome2tax[_.split('.')[0]]['family'] == 'Chromatiaceae')]
beta_well_known = [_ 
                   for _ in m1.index 
                   if (genome2tax[_.split('.')[0]]['class'] == 'Betaproteobacteria') and (genome2tax[_.split('.')[0]]['family'] == 'Nitrosomonadaceae')]
print(len(gamma_well_known),len(beta_well_known))
# 246 231 (20210714)
# 290 287 (20230408)
remaining_ids = []
#for gid in tqdm(all_ids):
for gid in gamma_well_known+beta_well_known+analyzed_genomes:   
    if exists(f"/home-user/thliao/data/NCBI/genbank/bacteria/{gid}"):
        continue
    remaining_ids.append(gid)
print(len(remaining_ids))  
ngd.main(odir='/home-user/thliao/data/NCBI',
         ids_list=remaining_ids)

existing_gids = open('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/20221013_all_related.gids').read().split('\n')
new_gids = existing_gids+gamma_well_known+beta_well_known+analyzed_genomes
new_gids = list(set(new_gids))

tax_color_scheme = {'Archaea': '#358f0f',
                    'Nitrospirae': '#edc31d',
                    'Alphaproteobacteria': '#8c613c',
                    'Betaproteobacteria': '#956cb4',
                    'Gammaproteobacteria': '#78fce0',
                    "Other-proteobacteria": "#F57F17",
                    'Chloroflexi': '#e41a1c',
                    'Actinobacteria': '#11FF11',
                    'Planctomycetes': '#FF66bb',
                    "Verrucomicrobia": "#2e91e5",
                    "NC10": "#74a45b",
                    }


cmds = []
for gene in ["amoA", "amoB", "amoC",
             ]:
    for g in new_gids:
        infaa = f"/home-user/thliao/data/NCBI/modified_data/direct_protein_files/{g}.faa"
        if not exists(infaa) :
            continue
        ofile = f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes/{g}_{gene}"
        if exists(ofile):
            continue
        cmd = f"blastp -query {infaa} -db /home-user/thliao/project/nitrogen_cycle/curated_genes/{gene}.faa -out {ofile} -num_threads 5 -outfmt 6 -max_target_seqs 10000000 -evalue 1e-20"
        cmds.append(cmd)
sbatch_all(cmds,batch_size=50,thread_per_tasks=5,prefix_name='blastp',fixed_cluster='others')

t = {'PF13435.9':'cycA','PF03264.17':'cycB'}
cmds = [] 
for g in new_gids:
    infaa = f"/home-user/thliao/data/NCBI/modified_data/direct_protein_files/{g}.faa"
    for _ in t:
        ofile = f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes/cyc/{g}_{_}.out"
        cmd = f"hmmsearch --cut_ga  --domtblout {ofile} --noali --notextw --cpu 10 /home-user/thliao/db/protein_db/Pfam.v33.1/separated/{_}.hmm {infaa} "
        if not exists(ofile):
            cmds.append(cmd)
sbatch_all(cmds,batch_size=50,thread_per_tasks=5,prefix_name='hmm',fixed_cluster='others')

from api_tools.tk import read_hmmsearch
cyc_locus2gene = {}
cyc_genome2gene = defaultdict(list)
for tab in tqdm(glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes/cyc/*.out')):
    genome,pf = tab.split('/')[-1].replace('.out','').rsplit('_',1)
    try:
        f_df = read_hmmsearch(tab)
        cyc_genome2gene[genome].append(list(f_df['target name']))
        for _ in list(f_df['target name']):
            cyc_locus2gene[_] = [t[pf],genome]
    except:
        continue
with open('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes/cyc/cyc_annotations','w') as f1:
    for k,(gene,genome) in cyc_locus2gene.items():
        f1.write(f"{k}\t{gene}\t{genome}\n")





removed_genes = open('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/removed_gene_ids').read().strip().split('\n')
related_genomes = set()
gene2locus_id = defaultdict(list)
for otab in glob("/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes/*"):
    if '_' not in otab.split('/'):
        continue
    genome,gene = otab.split('/')[-1].rsplit('_',1)
    for row in open(otab):
        if float(row.strip('\n').split("\t")[-2]) <= 1e-20:
            # if row.split('\t')[1] in removed_genes:
            #     # remove amo from heterotrophic nitrifier 
            #     continue
            related_genomes.add(genome)
            locus_id = row.split('\t')[0]
            gene2locus_id[gene].append(locus_id)
gene2locus_id = {k:set(v) for k,v in gene2locus_id.items()}

locus2seq = {}
for genome in related_genomes:
    faa = f"/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files/{genome}.faa"
    locus2seq.update({_.id:_ for _ in SeqIO.parse(faa,'fasta')})
locus2new = {}
locus2gene = {}
for gene in ["xmoA", "xmoB", "xmoC"]:
    ref_xmo = list(SeqIO.parse(f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/xmoCBA_remove_archaea/{gene}.faa','fasta'))
    l = gene2locus_id['a'+gene[1:]]
    for _ in l:
        locus2new[_]=['new']
    seqs = [locus2seq[_] for _ in l]
    seqs += list(ref_xmo)
    seqs = list({_.id:_ for _ in seqs}.values())
    seqs = [_ for _ in seqs if _.id.split('_amo')[0] not in removed_genes]
    for _ in seqs:
        if _.id[-1] in 'ABC':continue
        locus2gene[_.id] = gene
    with open(f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/{gene}.faa','w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
related_genomes = set([f"GCA_"+ _.split('_')[0].replace('v','.') for _ in locus2gene])

l2seq = {}
for gene in ["xmoA", "xmoB", "xmoC"]:
    infaa = f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/{gene}.faa'
    l2seq.update({_.id:_ for _ in SeqIO.parse(infaa,'fasta')})
with open('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/xmo.faa','w') as f1:
    SeqIO.write(l2seq.values(),f1,'fasta-2line')

   
os.system(f"/home-user/thliao/software/kofamscan/exec_annotation --cpu 10 -p /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/xmo.hal -k /mnt/home-db/pub/protein_db/kegg/v20230301/ko_list --tmp-dir ./tmpkofam -o /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/xmo.kofamout -f mapper-one-line --no-report-unannotated /mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/xmo.faa")
xmo2locus = defaultdict(list)
locus2xmo = {}
for row in open('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/xmo.kofamout'):   
    locus2xmo[row.split('\t')[0]] = row.split('\t')[1].strip('\n')
    xmo2locus[row.split('\t')[1].strip('\n')].append(row.split('\t')[0])

for gene,ko in [("xmoA",'K10944'), ("xmoB",'K10945'), ("xmoC",'K10946')]:
    infaa = f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/{gene}_afterko.faa'
    with open(infaa,'w') as f1:
        SeqIO.write([l2seq[_] for _ in xmo2locus[ko]],f1,'fasta-2line')
    aln = infaa.replace('_afterko.faa','.aln')
    os.system(f"mafft --auto {infaa} > {aln}")
    os.system(f"FastTree {aln} > {aln.replace('.aln','.newick')}")

_genome2tax = {}
gene2tax = {}
for l in locus2xmo:
    sl = l.split('v')[0]
    genome = f"GCA_{sl}"
    if genome not in genome2tax:continue
    row = genome2tax[genome]
    if row['phylum'] == 'Proteobacteria' and (str(row['class']).split('proteo')[0] in ['Alpha', 'Beta', 'Gamma']):
        gene2tax[l] = row['class']
    elif row['phylum'] == 'Proteobacteria':
        gene2tax[l] = 'Other-proteobacteria'
    elif row['phylum'] == "candidate division NC10":
        gene2tax[l] = "NC10"
    elif row['phylum'] in ['Crenarchaeota', 'Thaumarchaeota']:
        gene2tax[l] = 'Archaea'
    elif row['phylum'] in tax_color_scheme:
        gene2tax[l] = row['phylum']
    if l in gene2tax:
        _genome2tax[genome] = gene2tax[l]
text = to_color_strip(gene2tax,
                      {k:v for k,v in tax_color_scheme.items() if k in list(gene2tax.values())},
                      dataset_name='taxonomic',other_params={'width':'50','strip_width':'50'})
with open('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/extra_genomes_xmoCAB2023/gene_taxcolorstrip.txt','w') as f1:
    f1.write(text)
