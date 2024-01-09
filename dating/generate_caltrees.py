import itertools
from dating_workflow.calibrations.cal_gen import get_mapping
from final_params import strategies_params

mapping_dict = get_mapping()

### validations
from ete3 import Tree
gene2tree = {}
for gene,i in strategies_params.items():
    gene2tree[gene] = Tree(i[3],8)
# project unique, naming issue
# for gene,tre in gene2tree.items():
#     #print(gene)
#     for (node,scheme), (lca,name) in mapping_dict.items():
#         if lca == 'ROOT':continue
#         try:
#             if tre.get_common_ancestor(lca.split('|')):
#                 print('')
#         except:
#             pass

### main
import pandas as pd
from os.path import join
import re

df = pd.read_excel('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal/Table S2.xlsx', index_col=0)
df.columns = [_.split(' ')[0] for _ in df.columns]   
# df = df.reindex(columns=[_[0] for _ in mapping_dict])

pattern = f"([ABP][1-9]+)"

def get_caltxt(df,calnames,gene_scheme):
    groups = re.findall(pattern,calnames)
    subdf = df.loc[groups,:]
    caltxt = ''
    for g in groups:
        for n in subdf.columns:
            _caltxt = str(subdf.loc[g,n])
            if _caltxt == 'nan':
                continue
            v = [v for k,v in mapping_dict.items() if k[0] == n]
            if len(v)!=1:
                v = [v for k,v in mapping_dict.items() if k[0] == n and k[1]==gene_scheme]
            #print(g,n,v)
            lca = v[0][0]
            if lca:
                caltxt += f"{lca}\t{_caltxt}\t{n}\n" 
    return caltxt
from api_tools.for_tree.format_tree import add_cal_api
genes = ['G18', 'M24', 'P39', 'COG25']
gene2cals = {}
gene2cals['COG25'] = [f"B{i}" for i in range(1,10)]
gene2cals['G18'] = [f"B{i}" for i in range(1,10)] + [f"B{i}A{j}" for i in range(1,10) for j in [1,2]]
gene2cals['M24'] = [f"B{i}" for i in range(1,10)] + [f"B{i}A{j}" for i in range(1,10) for j in [1,2]]
gene2cals['P39'] = [f"B{i}" for i in range(1,10)] + [f"B{i}P{j}" for i in range(1,10) for j in [1,3]]

odir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal'
for gene,tre in gene2tree.items():
    in_newick = strategies_params[gene][3]
    for cal in gene2cals[gene]:
        if gene in ['P39']:
            gene_scheme = 'P'
        elif gene in ['COG25']:
            gene_scheme = 'B'
        else:
            gene_scheme = 'A'
        calibration_txt = get_caltxt(df,cal,gene_scheme)
        with open('./tmp.txt','w') as f1:
            f1.write(calibration_txt)
        with open(f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal_txt/{cal}.txt','w') as f1:
            f1.write(calibration_txt)
        add_cal_api(in_newick,
                    out_newick=join(odir,f"{gene}_{cal}.newick"),
                    calibration_txt='./tmp.txt',
                    format=8)
odir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal/itol_txt'
for gene,tre in gene2tree.items():
    in_newick = strategies_params[gene][3]
    cal = gene2cals[gene][-1]
    if gene in ['P39']:
        gene_scheme = 'P'
    elif gene in ['COG25']:
        gene_scheme = 'B'
    else:
        gene_scheme = 'A'
    calibration_txt = get_caltxt(df,cal,gene_scheme)
    with open('./tmp.txt','w') as f1:
        f1.write(calibration_txt)
    draw_cal_itol('./tmp.txt',f"{odir}/{cal}_{gene}")
    cmd = f"python ~/bin/format_newick.py itol-cal -c ./tmp.txt -o {odir}/{cal}_{gene}"
    os.system(cmd)


### sequential 
odir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal'
c = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/cal_latest/sequential_sp27_normal_hard.txt'

t = open(c).read().strip().split('\n')[4:]
t = get_caltxt(df,'B5') + '\n'.join(t)
with open('./tmp.txt','w') as f1:
    f1.write(t)
in_newick = strategies_params['M24'][3]
add_cal_api(in_newick,
            out_newick=join(odir,f"M24_B5SEQ1.newick"),
            calibration_txt='./tmp.txt',
            format=8)

t = open(c).read().strip().split('\n')[4:]
t = get_caltxt(df,'B5') + '\n'.join(t)
with open('./tmp.txt','w') as f1:
    f1.write(t)
in_newick = strategies_params['G18'][3]
add_cal_api(in_newick,
            out_newick=join(odir,f"G18_B5SEQ1.newick"),
            calibration_txt='./tmp.txt',
            format=8)


