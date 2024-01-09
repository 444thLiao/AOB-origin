import os
os.chdir("/mnt/home-backup/thliao/AOB/analysis/20210713_repeat")
import pandas as pd 
from tqdm import tqdm
from glob import glob
from ete3 import Tree
from os.path import *
from collections import Counter,defaultdict

near_genomes = open("./comparisons/near_genomes.txt").read().strip().split('\n')
near_genomes = [_.split('\t')[0] for _ in near_genomes if len(_.split('\t'))==2]
from bin.multiple_sbatch import sbatch_all
# cmds = []
# for gid in near_genomes:
#     tab = f'/mnt/maple/thliao/data/NCBI/modified_data/annotations/kegg/{gid}.tab'
#     if not exists(tab):continue
#     cmd = f"python /home-user/thliao/script/evol_tk/bin/kofamanno_s2.py -i {tab} -o /home-user/thliao/data/NCBI/modified_data/annotations/kofamscan_output "
#     if exists(f"/home-user/thliao/data/NCBI/modified_data/annotations/kofamscan_output/{gid}.kofamout"):continue
#     cmds.append(cmd)
# sbatch_all(cmds,batch_size=100,fixed_cluster='others')

# kegg_all_hmm = "/mnt/home-backup/thliao/kofam/20190810/kegg.HMM"   # version is 20190810
# cmds = []
# for gid in gids:
#     ofile = f'/mnt/maple/thliao/data/NCBI/modified_data/annotations/kegg/{gid}.tab'
#     faa = f"/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files/{gid}.faa"
#     if not exists(faa):
#         print(gid)
#         continue
#     if not exists(ofile):
#         cmd = f"hmmsearch --tblout {ofile} --acc --noali --notextw --cpu 10 {kegg_all_hmm} {faa} "
#         cmds.append(cmd)
# sbatch_all(cmds,thread_per_tasks=10,batch_size=40,fixed_cluster='others')

# other genes
from bin.multiple_sbatch import sbatch_all
db  = './comparisons/custom/ncyA_prot_ref.faa'
cmds = []
for genome in near_genomes:
    faa = f"/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files/{genome}.faa"
    cmd = f"blastp -query {faa} -db {db}  -evalue 1e-50 -num_threads 5 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' > ./comparisons/custom/nycA/{genome}.tab"
    if not exists(f"./comparisons/custom/nycA/{genome}.tab"):
        cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=5,batch_size=50)


## 
from api_tools.IO_for.read import parse_hmmscan
from ete3 import Tree
gids = Tree('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/xmoA_genome_tree/phylo_20220926/related_faa.tre.treefile').get_leaf_names()
near_genomes = open("./comparisons/near_genomes.txt").read().strip().split('\n')
near_genomes = [_.split('\t')[0] for _ in near_genomes if len(_.split('\t'))==2]
near_genomes += gids
gids = [_ for _ in set(near_genomes) if _.startswith('GC')]
# ref = pd.read_csv('gene_annotations/kegg/near_genomes_kegg_binary.tab',sep='\t',index_col=0)
# genomes = list(ref.index)
# missing_genomes = []
genome2ko2locus =  defaultdict(lambda:defaultdict(list))
odir = "/mnt/maple/thliao/data/NCBI/modified_data/annotations/kofamscan_output/"
for genome in tqdm(gids):
    ofile = join(odir,f'{genome}.kofamout')
    if not exists(ofile):
        ofile = join(odir,f"{genome.replace('GCF','GCA')}.tab")
        if not exists(ofile):
            print(genome)
            continue
    ko2locus = defaultdict(list)
    for row in open(ofile).read().strip().split('\n'):
        rows = row.split('\t')
        locus = rows[0]
        # kos = ';'.join(rows[1:])
        for ko in rows[1:]:
            ko2locus[ko].append(locus)
    genome2ko2locus[genome] = ko2locus
details_df = pd.DataFrame.from_dict({genome: {k:','.join(x) for k,x in vdict.items()} for genome,vdict in genome2ko2locus.items()},
                                    orient='index')
anno_f = "/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/near_genomes_kegg_details.tab"
details_df.to_csv(anno_f,sep='\t',index=1)


genome2ko2num = {genome: {k:len(set(x)) for k,x in vdict.items()} for genome,vdict in genome2ko2locus.items()}
num_df = pd.DataFrame.from_dict(genome2ko2num,orient='index')
num_df = num_df.fillna(0)
bin_df = num_df.where(num_df==0,1)
anno_f = "/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/near_genomes_kegg_binary.tab"
bin_df.to_csv(anno_f,sep='\t',index=1)
num_df.to_csv('/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/near_genomes_kegg_numeric.tab',sep='\t',index=1)
other_df = pd.read_csv('./comparisons/custom/ncyA_anno.tab',sep='\t',index_col=0)




"""


/home-user/thliao/software/kofamscan/exec_annotation -o ~/test2.kegganno.out -f mapper-one-line --no-report-unannotated ./test2.faa --cpu 10 -p /mnt/home-db/pub/protein_db/kegg/v20230301/profiles

python3 /home-user/thliao/bin/batch_run/batch_any.py -i /home-user/jjtao/db/kegg/profiles -s hmm -o LN1_tot_annotated_tab -ns tab -np 32 -cmd "hmmsearch --tblout {ofile} --acc --noali --notextw --cpu 1 {infile} LN1_tot.proteins.faa > ./log.txt"

"""
# ~ = /home-user/thliao
# ./ = 




import os
target_dirs = ['/mnt/ivy/nxiang/practice/metagenome/gene_centric_Sep4/total_reads',
               '/mnt/ivy/nxiang/practice/metagenome/gene_centric_Sep4/total_reads',]
output_dir = ''
for d in target_dirs:
    for r in os.listdir(d) :
        if r.endswith('.faa'):
            faa = r
            name = r.replace('.faa','')
            cmd = f"/home-user/thliao/software/kofamscan/exec_annotation -o {output_dir}/{name}.kofamout -f mapper-one-line --no-report-unannotated {d}/{faa} --cpu 5 -p /mnt/home-db/pub/protein_db/kegg/v20230301/profiles"
            job_file = f'{output_dir}/{name}_kofam.sh'
            with open(job_file,'w') as f1:
                f1.write('#!/bin/bash\n')
                f1.write(cmd)
            os.system(f"sbatch -n 5 {job_file}")
            



