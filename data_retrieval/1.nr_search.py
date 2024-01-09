"""
Retrieve xmoA homologues from the nr database

"""

import os
import pandas as pd
from collections import defaultdict 
from tqdm import tqdm
from Bio import SeqIO

os.chdir("/home-user/thliao/project/AOB/analysis/20230331_update")

# cmds = []
# for _ in 'ABC':
#     cmds.append(f"""blastp -query ./nr_search/ref_seqs/xmo{_}.faa -db /home-user/thliao/db/protein_db/NCBI/nr_20210713/blast_db/nr -out ./nr_search/nt_output/xmo{_}_retrieve_all.blast -max_target_seqs 10000000 -evalue 1e-20 -num_threads 9 -task blastp-fast -outfmt '6 qaccver saccver staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore' """)
cmds = []
for _ in 'ABC':
    cmds.append(f"""diamond_v2.1.6 blastp --quiet -q ./nr_search/ref_seqs/xmo{_}.faa --db /mnt/ivy/thliao/db/protein_db/NCBI/nr_20230325/diamond/nr.dmnd -o ./nr_search/diamond/xmo{_}_retrieve_all.blast -k 0 --evalue 1e-20 -p 10 --outfmt 6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp """)
    
# for xmoA 70502.36s user 5050.85s system 199% cpu

for unit in "ABC":
    faa_path = f"./nr_search/ref_seqs/xmo{unit}.faa"
    ref2length = {_.id:len(_.seq) for _ in SeqIO.parse(faa_path,'fasta')}
    # bewlow are pre-geneated using only 1e-20 threshold
    # subj2length = {p.split(' ')[0]:len(_.seq)
    #                for _ in SeqIO.parse(f"./nr_search/nt_output/unique_prot_xmo{_}.faa",'fasta')
    #                for p in _.description.split('\x01')}

    df = pd.read_csv(f"nr_search/diamond/xmo{unit}_retrieve_all.blast", sep='\t', header=None)
    print(df.shape)
    df.columns = "qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp".split(' ')
    df = df.loc[(df['evalue'] <= 1e-20),:]
    # accidentally added for xmoA
    df = df.loc[df['qseqid']!='Y14338.1',:]
    # df.loc[:,'length'] = [ref2length[q] for q in df['qaccver']]
    # df.loc[:,'cov'] = (df["qstart"] - df["qend"]).abs()/df.loc[:,'length'] *100
    # df.loc[:,'sublength'] = [subj2length[q] for q in df['saccver']]
    # df.loc[:,'subcov'] = (df["sstart"] - df["send"]).abs()/df.loc[:,'slen'] *100    
    df = df.loc[(df['scovhsp']>=60), :]
    sub_df = df.loc[:, ['sseqid', 'qstart', 'qend']]
    unique_id = sub_df['sseqid'].unique()
    unique_id = set(unique_id)
    print(len(unique_id))
    acc2seq_full = {}
    for seq in tqdm(SeqIO.parse('/mnt/ivy/thliao/db/protein_db/NCBI/nr_20230325//nr', 
                                format='fasta'),
                    total=535574845):
        if seq.id in unique_id:
            acc2seq_full[seq.id] = seq
        if len(acc2seq_full) == len(unique_id):
            break

    with open(f'./nr_search/diamond/unique_prot_xmo{unit}.list', 'w') as f1:
        f1.write('\n'.join(unique_id))

    ofaa = f'./nr_search/diamond/unique_prot_xmo{unit}.faa'
    with open(ofaa,'w') as f1:
        SeqIO.write(list(acc2seq_full.values()),f1,'fasta-2line')

# pid2
from bin.multiple_sbatch import batch_iter
from bin.ncbi_convertor.func_factory import *
from os.path import *
order_id_list = [_ for _ in open('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/unique_prot_xmoA.list').read().split('\n')]
existed_df = pd.read_csv('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/pid2taxdf.tab',sep='\t',index_col=0)
new_ids = set(order_id_list).difference(set(existed_df.index))

all_iters = batch_iter(order_id_list,500)
for i,n_iter in enumerate(all_iters):
    convertor = NCBI_convertor(n_iter,db='protein')
    ofile = f'/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/xmoA_basic.out/chunk{i}.out'
    if exists(ofile):
        continue
    pid2assembly_dict = convertor.pid2assembly()
    with open(ofile, 'w') as f1:
        print('#accession ID\tGI\tassembly_ID\tnuccore ID\tstart\tend\tstrand', file=f1)
        for pid in convertor.origin_ids:
            if convertor.GI is not None:
                GI = convertor.GI[pid]
            else:
                GI = 'NA'
            _dict = pid2assembly_dict[pid]
            for assembly_id,info in zip(_dict['assembly'],_dict['nuc_info']):
                print(f'{pid}\t{GI}\t{assembly_id}\t' + '\t'.join(map(str, info)), file=f1)
    tqdm.write('finish writing into ' + ofile)
from glob import glob
dfs = []
for i in glob('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/xmoA_basic.out/chunk*.out'):
    dfs.append(pd.read_csv(i, sep='\t'))
result_df = pd.concat(dfs,axis=0)
result_df.to_csv('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/xmoA_basic.out/pid2assembly.tab',sep='\t')

# missing_gids = set(order_id_list).difference(set(result_df['#accession ID']))

# assemblied = result_df.loc[~result_df['assembly_ID'].isna(),:]
# print(assemblied.shape)

all_iters = batch_iter(order_id_list,500)
for i,n_iter in enumerate(all_iters):
    convertor = NCBI_convertor(n_iter,db='protein')
    ofile = f'/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/xmoA_basic.out/chunk{i}.taxdf'
    if exists(ofile):
        continue
    pid2tid = convertor.get_taxon()
    pid2taxon = convertor.get_taxons_from_tid()
    pd.DataFrame.from_dict(pid2taxon).to_csv(ofile,sep='\t')
    tqdm.write('finish writing into ' + ofile)
dfs = []
for i in glob('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/xmoA_basic.out/chunk*.taxdf'):
    dfs.append(pd.read_csv(i, sep='\t',index_col=0).T)
result_df = pd.concat(dfs,axis=0)
result_df.to_csv('/mnt/ivy/thliao/project/AOB/analysis/20230331_update/nr_search/diamond/xmoA_basic.out/taxdf.tab',sep='\t')
