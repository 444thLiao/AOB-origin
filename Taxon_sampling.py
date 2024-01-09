# using the genomes downloaded from NOB project
# firstly, we annotated them
# according to the phylogeny with 1471 genomes
from Bio import SeqIO
from dating_workflow.step_script.quick_sampling import *
from dating_workflow.step_script.extract_cog25 import parse_annotation 
from os.path import *
from ete3 import Tree            
from glob import glob
from subprocess import check_call
from tqdm import tqdm
import os
os.chdir('/mnt/home-backup/thliao/AOB/analysis/update_calibrations')

genome2cog25 = parse_annotation(expanduser('~/data/NCBI/modified_data/cog25_annotate/'), 
                                top_hit=True, evalue=1e-20)
st = Tree('/mnt/home-backup/thliao/AOB/whole_tree_20200412/trees/fasttree/all_1471_new.newick',3)
target_nodes = []
sampling()


# perform generax for deduce the target nodes
base_odir = './'
locus2gene = {}
for of in glob('./gene_annotations/*_o/locus2genes.list'):
    _dict = {row.split('\t')[0]:row.split('\t')[1] 
             for row in open(of).read().strip('\n').split('\n')}
    locus2gene.update(_dict)

for faa in ["/mnt/home-backup/thliao/AOB/whole_tree_20200412/rebuild_gene_tree/cycA/from_genomes_annotated.faa",
            "/mnt/home-backup/thliao/AOB/whole_tree_20200412/rebuild_gene_tree/cycB/from_genomes_annotated.faa"]:
    records = list(SeqIO.parse(faa,'fasta'))
    name = faa.split('/')[-2]
    locus2gene.update({_.id:name for _ in records})

genome2locus = defaultdict(list)
for l in tqdm(locus2gene):
    genome = 'GCA_'+ l.split('_')[0].replace('v','.')
    genome2locus[genome].append(l)
    
id2seq = {}
for genome,locus_list in tqdm(genome2locus.items()):
    faa = f"/home-user/thliao/data/NCBI/modified_data/direct_protein_files/{genome}.faa"
    records = list(SeqIO.parse(faa,'fasta'))
    for l in locus_list:
        id2seq[l] = [r for r in records if r.id == l][0]
        
with open(join(base_odir,'gene_annotations','all_annotated_locus.faa'),'w') as f1:
    SeqIO.write(id2seq.values(),f1,'fasta-2line')



def build_genetree(ofaa,redo=False):
    gtree_f = f"{ofaa.replace('.faa','.newick')}"
    aln_f = f"{ofaa.replace('.faa','.trimal')}"
    cmd = f"""
ginsi {ofaa} > {ofaa.replace('.faa','.aln')}
trimal -automated1 -in {ofaa.replace('.faa','.aln')} -out {ofaa.replace('.faa','.trimal')}
FastTreeMP {ofaa.replace('.faa','.trimal')} > {gtree_f}
    """
    if not exists(gtree_f) or redo:
        check_call(cmd,shell=1)
    t = Tree(gtree_f,0)
    t.resolve_polytomy()
    t.write(outfile=gtree_f)
            
    return gtree_f,aln_f

id2seq = {r.id:r 
          for r in SeqIO.parse(join(base_odir,'gene_annotations','all_annotated_locus.faa'),
                               format='fasta')}

stree = "/mnt/home-backup/thliao/AOB/whole_tree_20200412/trees/fasttree/all_1471_new.newick"
leaves = Tree(stree,3).get_leaf_names()

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
              "cycA":['cycA'],
              "cycB":['cycB'],
              }

generax_odir = "/mnt/ivy/thliao/project/AOB/analysis/update_calibrations/generax/indata_full"
collect_parameters = []
redo = False
for gene in ['amoA','amoB','amoC','cycA','cycB']:
    locus = [_ 
            for _,v in locus2gene.items() 
            if v in homologues[gene] and "GCA_"+_.split('_')[0].replace('v','.') in leaves]
    l2g_list = join(generax_odir,f'{gene}_locus2genome.txt')
    if not exists(dirname(l2g_list)):
        os.system(f'mkdir -p {dirname(l2g_list)}')
    with open(l2g_list,'w') as f1:
        f1.write('\n'.join([f"{_l}\t{'GCA_'+ _l.split('_')[0].replace('v','.')}" for _l in locus]))
        
    ofaa = join(generax_odir,f"{gene}.faa")

    seqs = [id2seq[_] for _ in locus]
    for _ in seqs:
        _.name=_.description = ''
    SeqIO.write(seqs,
                open(ofaa,'w'),format='fasta-2line')
    gtree_f,aln_f = build_genetree(ofaa,redo=redo)
    
    collect_parameters.append((gene,gtree_f,aln_f,l2g_list))
        
with open(join(generax_odir,f"genrax.ctl"),'w') as f1:
    f1.write(f"""[FAMILIES] # this is a comment""" + '\n')
    for gene,sgtree,align,mapping in collect_parameters:
        f1.write(f"""
- {gene}
starting_gene_tree = {sgtree}
alignment = {align}
mapping = {mapping}
subst_model = LG+I+G4    """)
        
t = Tree(stree,3)
t.set_outgroup("I337_S1")
t.resolve_polytomy()

with open(join(generax_odir,'stree.newick'),'w') as f1:
    f1.write(t.write(format=3))
cmd = f"mpiexec -np 32 generax -s {join(generax_odir,'stree.newick')} -p {join(generax_odir,'output')} -f {join(generax_odir,f'genrax.ctl')} --seed 100 --unrooted-gene-tree"
print(cmd)




from dating_workflow.toolkit.fromGenerax import *
from dating_workflow.bin.parse_mcmc import *

gene2target_nodes = {}
for xml_file in glob("/mnt/home-backup/thliao/AOB/analysis/update_calibrations/generax/all_fixed/reconciliations/*.xml"):
    gene = basename(xml_file).split('_')[0]
    stree = join(base_odir,'generax','indata','stree.newick')
    #generax2mcmctree(xml_file,stree,gene,dating_o='./dating/reconciliation_o')
    phylum_name = xml_file.split('/')[-4]
    st =  Tree(stree,format=3)
    tmp_name = phylum_name+'_'+gene
    _p2node,_p2node_transfer_receptor = get_p2node(xml_file,stree,key=tmp_name)
    target_nodes = list(_p2node.values())[0] + list(_p2node_transfer_receptor.values())[0]
    gene2target_nodes[gene] = [_ for _ in target_nodes if _ and _!='NoName']  #NoName maybe the original root or the internal nodes added by resolve_polytomy
gene2target_nodes = {k:list(set(v)) for k,v in gene2target_nodes.items()}


# manual curated target nodes
target_nodes_info = {"I435_S1":{'name':'beta AMO','type':['AOB']},
                     "I426_S1":{'name':'gamma AMO+PMO','type':['AOB','MOB']},
                     "I554_S1":{'name':"gamma AMO",'type':['AOB']},
                     "I1030_S1":{'name':'Nitrospirae AMO','type':['AOB']},
                     "I359_S1":{'name':'NC10 PMO','type':['MOB']},
                     "I472_S1": {'name':"verruco PMO",'type':['MOB']},   
                     "I632_S1":{'name':"gamma PMO",'type':['MOB']}
                     }

target_nodes = list(target_nodes_info)
target_nodes.append('I164_S1')  # LCA of alpha-proteobacteria
target_nodes.append('I172_S1')  # LCA of proteobacteria
must_in_genomes = open("/mnt/home-backup/thliao/cyano/ref_genomes_list.txt").read().strip('\n').split('\n')
must_in_genomes.append('GCA_002842395.1')
cluster2genomes = get_cluster('/mnt/home-backup/thliao/AOB/analysis/update_calibrations/generax/indata/stree.cluster.list')
g2cluster = {v:c for c,d in cluster2genomes.items() for v in d}
retained_ids = sampling(st,
                        target_nodes,
                        must_in = must_in_genomes ,
                        node2cluster=g2cluster,
                        genome2cog25=genome2cog25,
                        max_num_up_each_level=1,
                        max_num_down=8)
print(len(retained_ids))

# manual adjust
retained_ids.remove('GCA_001788395.1')
retained_ids.remove('GCA_000212735.1')

text = to_binary_shape({g:['keep'] for g in retained_ids},
                    {"keep":{"color":"#88b719"}})
text = to_color_range({g:'keep' for g in retained_ids},
                        {"keep":"#88b719"} )
_odir = join(base_odir,'dating','sampling_ids')
if not exists(_odir):
    os.system(f"mkdir -p {_odir}")
with open(join(_odir,f'sampling.txt'),'w') as f1:
    f1.write(text)
with open(join(_odir,f'sampling.list'),'w') as f1:
    f1.write('\n'.join(retained_ids))

num_g = len(retained_ids)
final_tree = f"./trees/final/dating_generax_{num_g}g.newick"
final_phy = f"./trees/phy_files/{num_g}g_cog25_nucl.phy"
cmd = f"""
python3 ~/script/evolution_relative/dating_workflow/step_script/extract_bac120.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/bac120_annotate -o ./bac120_extract/dating_generax -evalue 1e-50 -gl {join(_odir,f'sampling.list')}

python3 ~/bin/batch_run/batch_mafft.py -i  ./bac120_extract/dating_generax -s faa -o ./bac120_extract/dating_generax -gl {join(_odir,f'sampling.list')} -m mafft -f ;
python3 ~/bin/batch_run/batch_trimal.py -i  ./bac120_extract/dating_generax -s aln -o ./bac120_extract/dating_generax
python3 /home-user/thliao/script/evolution_relative/dating_workflow/toolkit/concat_aln.py -i ./bac120_extract/dating_generax -o ./trees/concat/dating_generax.trimal -s trimal -gl {join(_odir,f'sampling.list')} -ct partition -no_graph
"""
check_call(cmd,shell=1)

cmd = f"""
mkdir -p ./trees/iqtree; iqtree -nt 30 -m LG -redo -pre ./trees/iqtree/dating_generax_guide -s ./trees/concat/dating_generax.trimal
mkdir -p ./trees/iqtree; iqtree -nt 30 -wbtl -bb 1000 -m LG+C20+F+G -redo -ft ./trees/iqtree/dating_generax_guide.treefile -s ./trees/concat/dating_generax.trimal -pre ./trees/iqtree/dating_generax
format_newick.py set1 -i ./trees/iqtree/dating_generax.contree -r "GCA_900552525.1,GCA_014377495.1" -o {final_tree} -f 0
"""
print(cmd)
check_call(cmd,shell=1)


cmd = f"""
## cog25 annotation and retrieval
python3 /home-user/thliao/script/evolution_relative/dating_workflow/step_script/extract_cog25.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/cog25_annotate -pd /mnt/home-backup/thliao/NCBI/modified_data/prokka_o -o ./cog25_single/nucl_e20 -evalue 1e-20 -gl {join(_odir,f'sampling.list')} -ot nucl 
python3 /home-user/thliao/script/evolution_relative/dating_workflow/step_script/extract_cog25.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/cog25_annotate  -o ./cog25_single/prot_e20 -evalue 1e-20 -gl {join(_odir,f'sampling.list')}

## preparation for prot
#python3 ~/bin/batch_run/batch_mafft.py -i ./cog25_single/prot_e20 -s faa -o ./cog25_single/dating_generax_prot -f -m ginsi -gl {join(_odir,f'sampling.list')};
#python3 ~/bin/batch_run/batch_trimal.py -i ./cog25_single/dating_generax_prot -o ./cog25_single/dating_generax_prot;
#python3 ~/script/evolution_relative/dating_workflow/toolkit/concat_aln.py -i ./cog25_single/dating_generax_prot -ct phy -gl {join(_odir,f'sampling.list')} -o ./trees/phy_files/dating_generax_cog25.trimal -s trimal -no_graph

## preparation for nucl
python3 ~/bin/batch_run/batch_mafft.py -i ./cog25_single/nucl_e20 -s ffn -o ./cog25_single/dating_generax_nucl -f -m mafft -gl {join(_odir,f'sampling.list')};
python3 ~/bin/batch_run/batch_trimal.py -i ./cog25_single/dating_generax_nucl -o ./cog25_single/dating_generax_nucl;
python3 ~/script/evolution_relative/dating_workflow/toolkit/concat_aln.py -i ./cog25_single/dating_generax_nucl -ct phy -gl {join(_odir,f'sampling.list')} -o {final_phy.rpartition('.')[0]}.trimal -s trimal -no_graph -p '1,2'

format_newick.py add-cal -i {final_tree} -o ./dating/cal_tree/{num_g}g_set10.newick -c ./dating/calibrations_set/cal_set10.txt -f 3

## dating !!!!! (using prot)
#python3 ~/script/evolution_relative/dating_workflow/bin/dating_pro.py -i ./trees/phy_files/dating_generax_cog25.phy -it ./dating/cal_tree/{num_g}g_set10.newick -o ./dating/dating_generax_{num_g}g/prot/template -p 1 -rg '1 30 1' -sf 5 -c 2  &
## dating !!!!! (using nucl)

python3 ~/script/evolution_relative/dating_workflow/bin/dating_pro.py -i {final_phy} -it ./dating/cal_tree/{num_g}g_set10.newick -o ./dating/dating_generax_{num_g}g/nucl/template -p 1 -rg '1 100 1' -sf 5 -c 2 -nucl 
"""
print(cmd)


lca_list = ';'.join([_.replace('|',',') for n,_ in node2LCA.items()])
name_list = ';'.join([target_nodes_info[n]['name'] for n,_ in node2LCA.items()])
cmd = f"""
python3 ~/script/evolution_relative/dating_workflow/bin/parse_mcmc.py -i ./dating/dating_generax_{num_g}g/nucl/clock2_diff_cal -ns '{lca_list}' -name '{name_list}' -o ./dating/dating_generax_{num_g}g/dated_results/nucl_clock2 -disable_plot
"""
check_call(cmd,shell=1)



# scheme2 using crown group calibration
g_l = open(join(_odir,f'sampling.list')).read().strip('\n').split('\n')
g_l.remove('GCA_900552525.1')
g_l.remove('GCA_001899315.1')
g_l.remove('GCA_014377495.1')
# above are the three genome outgroup for the oxygenic cyanobacteria
# long branch
# g_l.remove('GCA_001788395.1')

new_num_g = len(g_l)
new_list_file = join(_odir,f"CG_calir_{new_num_g}g.list")
with open(new_list_file,'w') as f1:
    f1.write('\n'.join(g_l))
print(len(g_l))


cmd = f"""

nw_prune {final_tree} -f {new_list_file} -v > ./trees/final/CG_calir_{new_num_g}g.newick

python3 ~/bin/batch_run/batch_mafft.py -i ./cog25_single/nucl_e20 -s ffn -o ./cog25_single/schemeCG_nucl -f -m mafft -gl {new_list_file};
python3 ~/bin/batch_run/batch_trimal.py -i ./cog25_single/schemeCG_nucl -o ./cog25_single/schemeCG_nucl;
python3 ~/script/evolution_relative/dating_workflow/toolkit/concat_aln.py -i ./cog25_single/schemeCG_nucl -ct phy -gl {new_list_file} -o ./trees/phy_files/schemeCG_cog25_nucl.trimal -s trimal -no_graph

format_newick.py add-cal -i ./trees/final/CG_calir_{new_num_g}g.newick -o ./dating/cal_tree/CG_{new_num_g}g_set10.newick -c ./dating/calibrations_set/CG_cal_set10.txt -f 3
python3 ~/script/evolution_relative/dating_workflow/bin/dating_pro.py -i ./trees/phy_files/schemeCG_cog25_nucl.phy -it ./dating/cal_tree/CG_{new_num_g}g_set10.newick -o ./dating/CG_{new_num_g}g/nucl/template -p 1 -rg '1 100 1' -sf 5 -c 2 -nucl
"""
print(cmd)




def get_new_node_name(ori_t,aft_t,target_nodes):
    retained_ids = aft_t.get_leaf_names()
    
    ori_name2node = {n.name:n for n in ori_t.traverse()}
    node2LCA = {}
    for n in target_nodes:
        node = ori_name2node[n]
        l1 = [l for l in node.children[0].get_leaf_names() if l in retained_ids][0]
        l2 = [l for l in node.children[1].get_leaf_names() if l in retained_ids][0]
        node2LCA[n] = f"{l1}|{l2}"
    return node2LCA

base_odir = '/mnt/home-backup/thliao/AOB/analysis/update_calibrations'

tagged_nodes = list(target_nodes_info)

stree = join(base_odir,'generax','indata','stree.newick')
st = Tree(stree,3)  
t = join(final_tree)
t = Tree(t,3)
node2LCA = get_new_node_name(st,t,tagged_nodes)
# manual adjust for some topologo inconsistence
node2LCA.update({"I426_S1":"GCA_002356115.1|GCA_001566965.1"})
from api_tools import pie_chart
text = pie_chart({k:{_t:1 for _t in _d['type']} for k,_d in target_nodes_info.items()},
                 {'AOB':'#B31B1B','MOB':"#2243B6"},
                 dataset_label='target nodes')
with open(join(base_odir,'itol','target_nodes_for_oritree.txt'),'w') as f1:
    f1.write(text)
    
text = pie_chart({v:{_t:1 for _t in target_nodes_info[k]['type']} for k,v in node2LCA.items()},
                 {'AOB':'#B31B1B','MOB':"#2243B6"},
                 dataset_label='target nodes')
with open(join(base_odir,'itol','target_nodes_for_sampled_tree.txt'),'w') as f1:
    f1.write(text)


lca_list = ';'.join([_.replace('|',',') for n,_ in node2LCA.items()])
name_list = ';'.join([target_nodes_info[n]['name'] for n,_ in node2LCA.items()])
cmd = f"""
python3 ~/script/evolution_relative/dating_workflow/bin/parse_mcmc.py -i ./dating/{num_g}g/nucl/clock2_diff_cal -ns '{lca_list}' -name '{name_list}' -o ./dating/{num_g}g/dated_results/nucl_clock2 -disable_plot
"""
check_call(cmd,shell=1)

# get only two node (beta LCA & gamma AOB+MOB)
lca_list = ';'.join(["GCA_002842395.1,GCA_003268955.1","GCA_000143085.1,GCA_001566965.1",
                     "GCA_003268955.1,GCA_001899235.1"])
name_list = ';'.join(["Beta LCA","Gamma AMO+PMO",'Beta AOB'])
cmd = f"""
python3 ~/script/evolution_relative/dating_workflow/bin/parse_mcmc.py -i ./dating/{num_g}g/nucl/clock2_diff_cal -ns '{lca_list}' -name '{name_list}' -o ./dating/{num_g}g/dated_results/nucl_clock2 -disable_plot
"""
check_call(cmd,shell=1)
##

cmd = f"""
python3 ~/script/evolution_relative/dating_workflow/figtree2itol.py -i {final_tree} -i2 ./dating/{num_g}g/nucl/clock2_diff_cal/{num_g}g_set24_run1/FigTree.tre -o ./dating/{num_g}g/dated_tree/{num_g}g_C21/{num_g}g_C21.newick; 
"""
check_call(cmd,shell=1)



