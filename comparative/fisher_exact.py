"""

"""

from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact,spearmanr
from bin.other_convertor.classify_kos import *
from os.path import *
from api_tools import to_binary_shape


def get_ko_info(ko):
    df = pd.read_csv('/mnt/home-backup/thliao/protein_db/kegg/ko_info.tab',
                     sep='\t', index_col=0, header=None)
    df.index = [_.split(':')[-1] for _ in df.index]
    koids = sep_ko(ko)

    sub_df = df.reindex(koids)
    text = ""
    for kid, row in sub_df.iterrows():
        des = row[1]
        name = des.split('; ')[0]
        text += '\t'.join([kid, name, '#ff77aa', '2', des]) + '\n'
    print(text)

def sep_ko(ko_list):
    kos = []
    if isinstance(ko_list, list):
        return ko_list
    for _ in ko_list.split('K'):
        if _.strip(string.printable[10:]).strip(' '):
            kos.append('K' + _.strip(string.printable[10:]).strip(' '))
    return kos

def fe_text(g1, g2, ko_df):
    ko2tab = {}
    ko2odd_ratio = {}
    ko2pvalue = {}
    for ko in tqdm(set(ko_df.columns)):
        ko_g1_p = sum(ko_df.loc[g1, ko] > 0)
        ko_g1_a = sum(ko_df.loc[g1, ko] == 0)
        ko_g2_p = sum(ko_df.loc[g2, ko] > 0)
        ko_g2_a = sum(ko_df.loc[g2, ko] == 0)
        tab_num = [[ko_g1_p, ko_g1_a],
                   [ko_g2_p, ko_g2_a]]
        ko2tab[ko] = tab_num
        oddratio, p = fisher_exact(tab_num, alternative="two-sided")
        ko2pvalue[ko] = p
        ko2odd_ratio[ko] = oddratio

    fe_corrected_p = multipletests([_
                                    for k, _ in ko2pvalue.items()],
                                   method='fdr_bh')
    fe_corrected_ko2p = dict(zip([k
                                  for k, _ in ko2pvalue.items()],
                                 fe_corrected_p[1]))
    sig_ko_list = {k: v
                   for k, v in fe_corrected_ko2p.items()
                   if v <= 0.05}
    return ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list

def output_sig_kos(sig_kos,X,y_bin,ofile):
    text = '\t'.join(['KO','description','Gain/Loss','corrected p-value'])+'\n'
    for ko,p in sig_kos.items():
        g1 = X.index[y_bin==1]
        g2 = X.index[y_bin==0]        
        ratio_AOB = X.loc[g1,ko].sum()/len(g1) # AOB
        ratio_others = X.loc[g2,ko].sum()/len(g2) # non-AOB
        if ratio_AOB > ratio_others:
            t = 'Gain'
        else:
            t = 'Loss'
        text += f"{ko}\t{get_ko_infos(ko).get(ko,ko)}\t{t}\t{p}\n"

    br_kos = ko_classified_br(sig_kos)
    md_kos = ko_classified_module(sig_kos)
    md2info = get_md_infos(md_kos.keys())
    info_kos = get_ko_infos(sig_kos)
    infos = get_br_info(br_kos)
    
    df = pd.read_csv(io.StringIO(text),sep='\t')
    df = df.sort_values(['Gain/Loss','KO'])
    
    for md, kos in tqdm(md_kos.items()):
        kos = [_ for _ in kos if _ in df.index]
        df.loc[kos, "module ID"] = md
        df.loc[kos, "module name"] = md2info[md]['name']
        df.loc[kos, "module all ko"] = md2info[md]['all ko']
                
    df.to_excel(ofile,index=0)

def main(annotation_df, tree, group1, group2, ofile=None,
         stop_at_annotation=False, just_g2=False,
         pass_id=False, return_df=False):
    """
    

    :param annotation_df:
    :param tree:
    :param group1: some specific group which is smaller
    :param group2: bigger group which normally comprise of group1
    :param ofile:
    :return:
    """
    if ofile is not None and not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))

    if pass_id:
        g1 = group1
        g2 = group2
    else:
        LCA = tree.get_common_ancestor(group1.split(','))
        larger_LCA = tree.get_common_ancestor(group2.split(','))

        g = larger_LCA.get_leaf_names()
        g1 = LCA.get_leaf_names()
        # target group
        if not just_g2:
            g2 = set(g).difference(set(g1))
        else:
            g2 = set(g)

    # remainning group
    ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list = fe_text(
        g1, g2, annotation_df)

    if stop_at_annotation:
        return ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list

    # summarized df
    br_kos = ko_classified_br(sig_ko_list)
    md_kos = ko_classified_module(sig_ko_list)
    md2info = get_md_infos(md_kos.keys())
    info_kos = get_ko_infos(sig_ko_list)
    infos = get_br_info(br_kos)

    if len(infos) == 0:
        raise Exception

    df = pd.concat(infos, axis=0)
    df.loc[:, 'des'] = [info_kos.get(_, '') for _ in df.index]

    new_df = pd.DataFrame(index=sig_ko_list, columns=['num_p in target', 'num_a in target', "ratio in target (%)",
                                                      'num_p in remained', 'num_a in remained', "ratio in remained (%)",
                                                      ])
    for ko in new_df.index:
        row = ko2tab[ko][0][::]
        row.append(round(row[0] / (row[0] + row[1]) * 100, 4))
        row += ko2tab[ko][1]
        row.append(round(row[3] / (row[4] + row[3]) * 100, 4))
        new_df.loc[ko, :] = row
    new_df = new_df.reindex(index=df.index)

    final_df = pd.concat([new_df, df], axis=1)
    for md, kos in tqdm(md_kos.items()):
        kos = [_ for _ in kos if _ in final_df.index]
        try:
            final_df.loc[kos, "module ID"] = md
            final_df.loc[kos, "module name"] = md2info[md]['name']
            final_df.loc[kos, "module all ko"] = md2info[md]['all ko']
        except:
            pass
            print(kos)
    final_df = final_df.sort_values(
        ['num_p in target', 'top', 'A', 'B', 'C', 'des'])

    final_df.loc[:, 'raw pvalue'] = [ko2pvalue[_] for _ in final_df.index]
    final_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_]
                                           for _ in final_df.index]

    if return_df:
        return final_df
    else:
        final_df.to_excel(ofile, index_label='K number')


from ete3 import Tree
anno_f = "/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/1444genomes_kegg_binary.tab"
anno_df = pd.read_csv(anno_f,sep='\t',index_col=0)
all_gids = list(anno_df.index)
# print(len(set(all_gids)))

stree = "/home-backup/thliao/AOB/analysis/20210713_repeat/ref_phy/sampled3/trees/iqtree/bac120_over20p_renamed.newick"
st = Tree(stree, format=3)
name2node = {n.name:n for n in st.traverse()}
all_genomes = st.get_leaf_names()
assert set(all_gids) == set(all_genomes)

_draw_df = pd.read_csv("/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/info_df.tab",sep='\t',index_col=0)
sub_df = _draw_df.reindex(all_genomes)

g1 = sub_df.index[(sub_df['interested lineages']=='Nitrospirae') & (sub_df['operon pro'] !='Nitro-AMO')]
g2 = sub_df.index[sub_df['operon pro'] =='Nitro-AMO'  ]
main(anno_df,
     Tree(stree, 3),
     group1=g1,
     group2=g2,
     pass_id=True,
     ofile=f'./comparisons/kegg/comammox.xlsx')

# the format shoule as follow:
# former is the bigger lineage, the later is the functional lineage
compare_groups = {"gamma-pmo2amo": ['I131_S100','I246_S100'],
                  "gamma-non2xmo": ['I89_S87','I131_S100'],
                  "Nitropsira": ["I37_S100","I92_S100"],
                  "Beta": ["I54_S100","I104_S100"],}
for name,n_list in compare_groups.items():
    n1,n2 = [name2node[_] for _ in n_list]
    main(anno_df,
     Tree(stree, 3),
     group1=set(n1.get_leaf_names()).difference(n2.get_leaf_names()),
     group2=n2.get_leaf_names(),
     pass_id=True,
     ofile=f'./comparisons/kegg/{name}.xlsx')


#! pearson correlation
from statsmodels.stats.multitest import multipletests
import numpy as np
def get_related_ko(y_bin, sub_anno_df, final_P=1e-50):
    ko2p = {}
    ko2coef = {}
    for ko, col in tqdm(sub_anno_df.iteritems()):
        test = spearmanr(y_bin, col)
        ko2p[ko] = test.pvalue
        ko2coef[ko] = test.correlation
    ko2p = {ko: v for ko, v in ko2p.items() if str(v) != 'nan'}
    fe_corrected_p = multipletests([_ for k, _ in ko2p.items()],
                                   method='fdr_bh')
    fe_corrected_ko2p = dict(zip([k
                                  for k, _ in ko2p.items()],
                                 fe_corrected_p[1]))
    sig_ko_list = {k: v
                   for k, v in fe_corrected_ko2p.items()
                   if v <= final_P}
    # print(len(sig_ko_list))
    return fe_corrected_ko2p, ko2coef

from scipy.stats import spearmanr
sub_idx = _draw_df.index[_draw_df['interested lineages'].isin(['beta','gamma'])]
y_bin = np.array([1 if _ == 'AMO' else 0 for _ in _draw_df.reindex(sub_idx)['operon'] ])
other_df = pd.read_csv('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/comparisons/custom/ncyA_anno.tab',sep='\t',index_col=0)
sub_anno_df = anno_df.loc[sub_idx,:]
sub_anno_df = pd.concat([other_df.reindex(sub_idx),sub_anno_df],axis=1).fillna(0)
ko2p,ko2coef = get_related_ko(y_bin,sub_anno_df,0.05)
ko2p = {ko:v for ko,v in ko2p.items() if abs(ko2coef[ko])>=0.4}
output_sig_kos(ko2p,sub_anno_df,y_bin,'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/comparisons/kegg/AOB_signatures.xlsx')



sub_idx = _draw_df.index[_draw_df['interested lineages']=='Nitrospirae']
y_bin = np.array([1 if _ == 'Nitro-AMO' else 0 for _ in _draw_df.reindex(sub_idx)['operon pro'] ])
other_df = pd.read_csv('./comparisons/custom/ncyA_anno.tab',sep='\t',index_col=0)
sub_anno_df = anno_df.loc[sub_idx,:]
sub_anno_df = pd.concat([other_df.reindex(sub_idx),sub_anno_df],axis=1)
ko2p,ko2coef = get_related_ko(y_bin,sub_anno_df,0.05)
ko2p = {ko:v for ko,v in ko2p.items() if abs(ko2coef[ko])>=0.4}
output_sig_kos(ko2p,ko2coef,'./comparisons/kegg/comammox_signatures.xlsx')



sub_idx = _draw_df.index[_draw_df['interested lineages']=='gamma']
y_bin = np.array([1 if 'AMO' in _ else 0 for _ in _draw_df.reindex(sub_idx)['operon pro'] ])
other_df = pd.read_csv('./comparisons/custom/ncyA_anno.tab',sep='\t',index_col=0)
sub_anno_df = anno_df.loc[sub_idx,:]
sub_anno_df = pd.concat([other_df.reindex(sub_idx),sub_anno_df],axis=1)
ko2p,ko2coef = get_related_ko(y_bin,sub_anno_df.fillna(0),0.05)
ko2p = {ko:v for ko,v in ko2p.items() if abs(ko2coef[ko])>=0.4}
output_sig_kos(ko2p,ko2coef,'./comparisons/kegg/gammaAOB_signatures.xlsx')



sub_idx = _draw_df.index[(_draw_df['interested lineages']=='gamma') & (_draw_df['operon'].isin(['AMO','PXM;PMO','PMO']))]
y_bin = np.array([1 if 'AMO' in _ else 0 for _ in _draw_df.reindex(sub_idx)['operon'] ])
other_df = pd.read_csv('./comparisons/custom/ncyA_anno.tab',sep='\t',index_col=0)
sub_anno_df = anno_df.loc[sub_idx,:]
sub_anno_df = pd.concat([other_df.reindex(sub_idx),sub_anno_df],axis=1)
ko2p,ko2coef = get_related_ko(y_bin,sub_anno_df.fillna(0),0.05)
ko2p = {ko:v for ko,v in ko2p.items() if abs(ko2coef[ko])>=0.4}
output_sig_kos(ko2p,ko2coef,'./comparisons/kegg/gammaAOBvsMOB_signatures.xlsx')


from api_tools import to_color_strip
target_ko = ['K00018']
s = anno_df.loc[:,target_ko].to_dict(orient='index')
text = to_color_strip({g:list(v)[0] 
                        for g,v in s.items() if v[target_ko[0]]==1},
                       {k:"#ff0000" for k in target_ko},
                       other_params={'MARGIN':'50',"STRIP_WIDTH":"50"})

target_ko = ['K00368','K09816','K09815','K09817','K03320']
s = anno_df.loc[:,target_ko].to_dict(orient='index')
text = to_binary_shape({g:[_ for _,a in v.items() if a!=0]
                        for g,v in s.items() },
                       {k:{'color':"#ff0000"} for k in target_ko},
                       unfilled_other=True,
                       other_params={'MARGIN':'50',"HEIGHT_FACTOR":"3",
                                     'SYMBOL_SPACING':'10'})
with open('./tmp.txt','w') as f1:f1.write(text)



