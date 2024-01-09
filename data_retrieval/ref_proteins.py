from Bio import SeqIO
from collections import defaultdict
import pandas as pd

def add_ref_seq(sub_df, used_records):
    used_ids = [str(_.id) for _ in used_records]
    id2info = {}
    record_need_dropped_ids = []
    for _, row in sub_df.iterrows():
        aa_id = row['AA accession'].strip()
        gene_name = row['gene name']
        seq = row['seq']
        info = row['phylum/class']
        if aa_id in used_ids:
            record_need_dropped_ids.append(aa_id)

        used_records.append(SeqIO.read(io.StringIO(
            f'>{aa_id}_{gene_name}\n{seq}'), format='fasta'))
        id2info[f'{aa_id}_{gene_name}'] = info
    final_records = [_ for _ in used_records if str(
        _.id) not in record_need_dropped_ids]
    return final_records, id2info


outgroup_gene_names = {'K10944': ['bmo'],
                       'K10945': ['bmo'],
                       'K10946': ['bmo']}
ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'
ref_df = pd.read_excel(ref_file, index_col=None)
ref_df = ref_df.loc[ref_df.loc[:, 'note'] != 'removed', :] 
final_name2grouping = defaultdict(set)
org2uniq_int = {}
uniq_int2org = {}
for ko in outgroup_gene_names.keys():
    sub_ref_df = ref_df.loc[ref_df.loc[:,
                                        'outgroup/ref for which KO'] == ko, :]
    unit = sub_ref_df.iloc[0,1][-1]  # subunit of amoCAB
    df = sub_ref_df.copy()
    orgs = [org for org in df['org name']]
    org_l1 = [org for org in orgs if orgs.count(org)==1]  # organism occur only once
    sub_df = df.loc[df['org name'].isin(org_l1)]  # data with organism occur only once
    
    for idx,row in sub_df.iterrows():
        aa_id = row['AA accession'].strip()
        gene_name = row['gene name']
        org_name = row['org name']
        # for the first one
        if not org2uniq_int:
            uniq_int = 1
            org2uniq_int[org_name] = uniq_int
            uniq_int2org[uniq_int] = org_name
        # for others
        if org_name in org2uniq_int:
            f_name = f"ref{org2uniq_int[org_name]}"
        else:
            uniq_int = max(org2uniq_int.values())+1
            org2uniq_int[org_name] = uniq_int
            uniq_int2org[uniq_int] = org_name
            f_name = f"ref{uniq_int}"
        final_name2grouping[f_name].add(f"{aa_id}_{gene_name}")
        
unit2df = {}
names_lt_2 = [uniq_int2org[int(k.replace('ref',''))] 
              for k,v in final_name2grouping.items() if len(v)>=2]
for ko in ['K10944','K10945','K10946']:
    _sub_ref_df = ref_df.loc[ref_df.loc[:,
                                        'outgroup/ref for which KO'] == ko, :]
    
    sub_ref_df = _sub_ref_df
    # sub_ref_df = _sub_ref_df.loc[_sub_ref_df.loc[:,
    #                                              'phylum/class'] != 'Thaumarchaeota', :]
    unit = sub_ref_df.iloc[0,1][-1]  # subunit of amoCAB
    sub_df = sub_ref_df.loc[sub_ref_df['org name'].isin(names_lt_2)]  # data with organism occur only once
    unit2df[unit] = sub_df
    
# unit2df.keys()
# dict_keys(['A', 'B', 'C'])
    