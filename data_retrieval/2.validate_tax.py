"""
This script include following processes
1. map clustered protein to taxonomic linages
2. download xmoA-containing genomes
3. blastp search the presence of cycAB, nxrAB and amoCBA
"""

from collections import defaultdict
from re import I
from subprocess import check_call
from Bio.SeqIO import parse
import pandas as pd
import os
from for_software.for_cd_hit.parse import parse_clstr
from collections import defaultdict
os.chdir("/mnt/home-backup/thliao/AOB/analysis/20230331_update")

cmd = f"cd-hit -i ./unique_prot_xmoA.faa -o ./prot_cdhit90_xmoA"
cmd = f"/home-user/thliao/software/kofamscan/exec_annotation -o -k /mnt/home-db/pub/protein_db/kegg/v20230301/ko_list -p /mnt/home-db/pub/protein_db/kegg/v20230301/profiles --cpu 30 -f mapper-one-line  --no-report-unannotated ./prot_cdhit90_xmoA"
cluster2seqs,cluster2repr = parse_clstr("/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/prot_cdhit90_xmoA.clstr")
repr2seqs = {}
for k,v in cluster2repr.items():
    repr2seqs[v] = cluster2seqs[k]
    
kofam_validated_seqs = [_.split('\t')[0] 
                        for _ in open(f'/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/cdhit90_xmoA.kofamout').readlines()
                        if _ and _.strip('\n').split('\t')[1] == 'K10944']
print(len(set(kofam_validated_seqs)))  # 887   --> 929 (2023)
# brocast the kofam annotation following the cd-hit
kofam_validated_seqs = set(kofam_validated_seqs)
xmoA_seqs_ = [_ for k,v in repr2seqs.items() for _ in v if k in kofam_validated_seqs]

seqs = list(SeqIO.parse('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/unique_prot_xmoA.faa','fasta'))
seq2seqs = {_.id: [k.split(' ')[0] for k in _.description.split('\x01')] for _ in seqs}
xmoA_seqs = [v for _ in xmoA_seqs_ for v in seq2seqs[_]]
print(len(xmoA_seqs))  # 109092    ---> 112385
pid2assembly = pd.read_csv('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/xmoA_basic.out/pid2assembly.tab',sep='\t',index_col=0)
xmoA_assembly = pid2assembly.loc[pid2assembly.index.isin(xmoA_seqs_)]

used_genomes = [_.replace('GCF','GCA') for _ in xmoA_assembly['assembly_ID']
                if str(_)!='nan']
used_genomes = list(set(used_genomes))
len(used_genomes)  # 887 --> 822
tax_df = pd.read_csv('/home-user/thliao/.cache/ncbi-genome-download/taxonomy.tab',sep='\t',index_col=0)
tax_df = tax_df.reindex([_.split('.')[0] for _ in used_genomes])
tax_df.index = used_genomes
missing_tax = tax_df.index[tax_df['species'].isna()]

convertor = NCBI_convertor(missing_tax,db='assembly')
pid2tid = convertor.get_taxon()
pid2taxon = convertor.get_taxons_from_tid()
sub_df = pd.DataFrame.from_dict(pid2taxon).T
sub_df = sub_df.reindex(missing_tax)
final_tax_df = pd.concat([tax_df,sub_df],axis=0)
final_tax_df.to_csv('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/assembly2taxon.tsv',sep='\t',index=1)

####previous one
operon_info = pd.read_excel('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/xmoCBA/operon_info.xlsx',index_col=0)
all_genomes = operon_info['genome'].unique()
a = set(all_genomes).difference(set(used_genomes))
b = set(used_genomes).difference(set(all_genomes))
# a is because the previous one has download allgenomes. some genomes doesn't have protein deposited in the nr data
# b is the updated one
# ! download
from bin.ngd import *
missing_gids = [_ for _ in used_genomes if not exists(f'/home-user/thliao/data/NCBI/modified_data/prokka_o/{_}')]
ngd.download(**{"assembly_accessions": ','.join(missing_gids),
                        "dry_run": False,
                        "use_cache":True, # to avoid it automatic download/update the summary file
                        "section": 'genbank',
                        "parallel": 2,
                        "output": db_dir,
                        "groups": 'bacteria',  # if not assign this, it will take long time to iterate all groups
                        "file_formats": 'genbank'})


# info_df = pd.DataFrame.from_dict(aid2taxon,orient='index')
# for aid,xmoA_seqs in assembly2xmoA_seqs.items():
#     info_df.loc[aid,'xmoA locus'] = ';'.join(sorted(xmoA_seqs))
# info_df = info_df.reindex(columns=['superkingdom', 
#                                    'phylum', 
#                                    'class',
#                                    'order', 
#                                    'family',
#                                    'genus',
#                                    'species', 'xmoA locus'])
# info_df = info_df.sort_values(list(info_df.columns))
# info_df.to_excel('./xmoA_containing_genomes.xlsx',index=1)

# info_df = pd.read_excel('./xmoA_containing_genomes.xlsx',index_col=0)
## validate the relationship between locus and genome
### in sometimes, the identical proteins from different genomes with different locus id will wrongly result in different genomes share same locus id. 


## ! validate the presence of cycAB
# from os.path import basename,dirname
# os.chdir('/mnt/home-backup/thliao/AOB/analysis/20210713_repeat')
# info_df = pd.read_excel('./xmoA_containing_genomes.xlsx',index_col=0)
# #bac_info_df = info_df.loc[info_df['superkingdom']=='Bacteria',:]

# genome_dir = './related_faa'
# if not exists(genome_dir):
#     os.system(f"mkdir {genome_dir}")
    
# missing_aids = []    
# for aid in info_df.index:
#     db_faa = f"/home-user/thliao/data/NCBI/modified_data/direct_protein_files/{aid}.faa"
#     if exists(db_faa) and not exists(join(genome_dir,basename(db_faa))):
#         os.system(f"ln -sf `realpath {db_faa}` {genome_dir}/")
#     elif exists(join(genome_dir,basename(db_faa))):
#         pass
#     else:
#         missing_aids.append(aid)
 
        
        
genome_dir = '/mnt/ivy/thliao/project/AOB/analysis/20230331_update/related_faa'
cmd = f"cat {genome_dir}/*.faa > ./merged.faa"
if not exists('./merged.faa'):
    check_call(cmd,shell=1)

tmp_dir = './related_faa/annotations'
collected_gs = '/home-user/thliao/project/nitrogen_cycle/curated_genes/'
# genome2collect_genes = defaultdict(list)
gene2locus_id = defaultdict(list)
genes = ["amoA", "amoB", "amoC", "cycA", "cycB", "hao",'nxrA','nxrB']
do_genes = glob(join(collected_gs, '*.faa'))
do_genes = [_ for _ in do_genes if basename(_).replace('.faa', '') in genes]
cmds = []
for db in do_genes:
    gene_name = basename(db).replace('.faa', '').strip()
    otab = join(tmp_dir, f'{gene_name}.tab')
    if not exists(dirname(otab)):
        os.makedirs(dirname(otab))
    if (not exists(otab)):
        cmd = f"blastp -query ./merged.faa -db {db} -out {otab} -num_threads 2 -outfmt 6 -max_target_seqs 10000000 -evalue 1e-3 "
        cmds.append(cmd)
    for row in open(otab):
        if float(row.strip('\n').split("\t")[-2]) <= 1e-20:
            if row.split('\t')[0] == 'Y14338.1':
                # remove amo from heterotrophic nitrifier 
                continue
            genome_name = "GCA_" + row.split("\t")[0].split('_')[0].replace('v','.')
            locus_id = row.split("\t")[0]
            # genome2collect_genes[genome_name].append(gene_name)
            gene2locus_id[gene_name].append(locus_id)
with open(f'{tmp_dir}/gene2locus_id','w') as f1:
    for gene,locus_list in gene2locus_id.items():
        for l in set(locus_list):
            f1.write(f"{gene}\t{l}\n")

