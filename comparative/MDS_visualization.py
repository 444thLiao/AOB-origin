"""
load required information
calculate the pairwise distance 


"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn import decomposition, manifold
from os.path import exists


def confidence_ellipse(x, y, n_std=1.96, size=100):
    """
    Get the covariance confidence ellipse of *x* and *y*.
    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.
    size : int
        Number of points defining the ellipse
    Returns
    -------
    String containing an SVG path for the ellipse
    
    References (H/T)
    ----------------
    https://matplotlib.org/3.1.1/gallery/statistics/confidence_ellipse.html
    https://community.plotly.com/t/arc-shape-with-path/7205/5
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    theta = np.linspace(0, 2 * np.pi, size)
    ellipse_coords = np.column_stack([ell_radius_x * np.cos(theta), ell_radius_y * np.sin(theta)])
    
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    x_scale = np.sqrt(cov[0, 0]) * n_std
    x_mean = np.mean(x)

    # calculating the stdandard deviation of y ...
    y_scale = np.sqrt(cov[1, 1]) * n_std
    y_mean = np.mean(y)
  
    translation_matrix = np.tile([x_mean, y_mean], (ellipse_coords.shape[0], 1))
    rotation_matrix = np.array([[np.cos(np.pi / 4), np.sin(np.pi / 4)],
                                [-np.sin(np.pi / 4), np.cos(np.pi / 4)]])
    scale_matrix = np.array([[x_scale, 0],
                            [0, y_scale]])
    ellipse_coords = ellipse_coords.dot(rotation_matrix).dot(scale_matrix) + translation_matrix
        
    path = f'M {ellipse_coords[0, 0]}, {ellipse_coords[0, 1]}'
    for k in range(1, len(ellipse_coords)):
        path += f'L{ellipse_coords[k, 0]}, {ellipse_coords[k, 1]}'
    path += ' Z'
    return path



anno_f = "/mnt/home-backup/thliao/AOB/analysis/20210713_repeat/gene_annotations/kegg/near_genomes_kegg_binary.tab"
anno_df = pd.read_csv(anno_f,sep='\t',index_col=0)
all_gids = list(anno_df.index)

# taxonomic info
metadata = "/home-user/thliao/.cache/ncbi-genome-download/taxonomy.tab"
tax_df = pd.read_csv(metadata,sep='\t',index_col=0)
tax_df.index = [_.split('_')[-1] for _ in tax_df.index]
sub_tax_df = tax_df.loc[~tax_df.index.duplicated(),:]
sub_tax_df = sub_tax_df.reindex([_.split('.')[0].split('_')[-1] for _ in all_gids])
sub_tax_df.index = all_gids

request_taxon = 'phylum'
gid2taxon = {}
for g, r2n in sub_tax_df.to_dict(orient='index').items():
    gid2taxon[g] = r2n.get('phylum', '')
    if r2n.get('phylum', '') == 'Proteobacteria' and r2n.get('class', ''):
        gid2taxon[g] = r2n['class']
    elif r2n.get('phylum', '') == "candidate division NC10":
        gid2taxon[g] = "NC10"
gid2taxon = {k: v for k, v in gid2taxon.items() if v}

operon_info = pd.read_excel('./gene_annotations/xmoCBA/operon_info.xlsx',index_col=0)

cb_dist = squareform(pdist(anno_df, metric="cityblock"))
kegg_cb_dist_df = pd.DataFrame(
    cb_dist, index=anno_df.index, columns=anno_df.index
)
kegg_cb_dist_df.to_csv(
    f"./comparisons/pairwise_distance/Kegg_bin_cityblock_distant.tab",
    sep="\t",
    index=1,
    index_label="Assembly IDs",
)


otab = f"./comparisons/pairwise_distance/kegg_MDS_cityblock.tab"
mds = manifold.MDS(dissimilarity="precomputed", n_jobs=-1,random_state=123)
projected_X = mds.fit_transform(kegg_cb_dist_df.values)
_draw_df = pd.DataFrame(projected_X[:, :2])
_draw_df.index = anno_df.index
_draw_df.to_csv(otab, sep="\t", index=1, index_label="Assembly IDs")


if not exists(
    f"./comparisons/pairwise_distance/Kegg_bin_cityblock_distant.tab"
):
    cb_dist = squareform(pdist(anno_df, metric="cityblock"))
    kegg_cb_dist_df = pd.DataFrame(
        cb_dist, index=anno_df.index, columns=anno_df.index
    )
    kegg_cb_dist_df.to_csv(
        f"./comparisons/pairwise_distance/Kegg_bin_cityblock_distant.tab",
        sep="\t",
        index=1,
        index_label="Assembly IDs",
    )
else:
    kegg_cb_dist_df = pd.read_csv(
        f"./comparisons/pairwise_distance/Kegg_bin_cityblock_distant.tab",
        sep="\t",
        index_col=0,
    )
    
otab = f"./comparisons/pairwise_distance/kegg_MDS_cityblock.tab"
if not exists(otab):
    mds = manifold.MDS(dissimilarity="precomputed", n_jobs=-1)
    projected_X = mds.fit_transform(kegg_cb_dist_df.values)
    _draw_df = pd.DataFrame(projected_X[:, :2])
    _draw_df.index = anno_df.index
    _draw_df.to_csv(otab, sep="\t", index=1, index_label="Assembly IDs")
else:
    _draw_df = pd.read_csv(otab, sep="\t", index_col=0)
    if 0 not in _draw_df.columns:
        new_cols = []
        for c in _draw_df.columns:
            try:
                new_cols.append(int(c))
            except:
                new_cols.append(c)
        _draw_df.columns = new_cols
operon_containing_genomes = [_ 
                             for _ in _draw_df.index 
                             if _ in list(operon_info.index)]
_draw_df.loc[:,'operon'] = [operon_info.loc[[_],'inferred xmo type'][0]
                            if _ in operon_containing_genomes else 'no' 
                            for _ in _draw_df.index]
interested_lineages = []
for gid in _draw_df.index:
    if gid2taxon[gid] == 'Alphaproteobacteria':
        interested_lineages.append('alpha')
    elif gid2taxon[gid] == 'Betaproteobacteria':
        interested_lineages.append('beta')
    elif gid2taxon[gid] == 'Gammaproteobacteria':
        interested_lineages.append('gamma')
    elif gid2taxon[gid] == 'Nitrospirae':
        interested_lineages.append('Nitrospirae')
    else:
        interested_lineages.append('others')
_draw_df.loc[:,'interested lineages'] = interested_lineages
tmp = []
for g,row in _draw_df.iterrows():
    o = row['operon'] 
    if o == 'AMO' and row['interested lineages'] == 'gamma':
        tmp.append('gamma-AMO')
    elif o == 'AMO' and row['interested lineages'] == 'beta':
        tmp.append('beta-AMO')
    elif o == 'AMO' and row['interested lineages'] == 'Nitrospirae':
        tmp.append('Nitro-AMO')
    else:
        tmp.append(o)
_draw_df.loc[:,'operon pro'] = tmp
_draw_df.loc[:,'show_text'] = [gid2taxon[gid] for gid in _draw_df.index]



#! example of draing MDS with ellipse
## more interactive scripts are provided in the `ipython_notebook/AOB viewing.ipynb` 


from IPython.display import Image
import plotly.graph_objects as go

sub_d = _draw_df.sort_values("interested lineages",ascending=True)
c_map = {"AMO":"#e53238",
        "no":"#B0BEC5",
        "EMO":"#d500f9",
        "PMO":"#0064d3",
        "PXM":"#ff6f00",
        "PXM;PMO":"#0064d3"}
fig = go.Figure()
for g,idx in sorted(sub_d.groupby('operon').groups.items(),reverse=True):
    fig.add_scatter(x=sub_d.loc[idx,0],
                   y=sub_d.loc[idx,1],
                   mode='markers',
                    name=g,
                   marker=dict(size=8,color=c_map[g]))
fig.layout.height = 500
fig.layout.width = 700
fig.layout.xaxis.title.text = 'MDS1'
fig.layout.yaxis.title.text = 'MDS2'
x = _draw_df.loc[_draw_df['interested lineages'] == 'Nitrospirae', 0]
y = _draw_df.loc[_draw_df['interested lineages'] == 'Nitrospirae', 1]
fig.add_shape(type='path',
              path=confidence_ellipse(x, y),
              line=dict(color="#edc31d",
                        width=3,
                        dash="dot", ),
              )
x = _draw_df.loc[_draw_df['interested lineages'] == 'beta', 0]
y = _draw_df.loc[_draw_df['interested lineages'] == 'beta', 1]
fig.add_shape(type='path',
              path=confidence_ellipse(x, y),
              line=dict(color="#956cb4",
                        width=3,
                        dash="dot", ),
              )

x = _draw_df.loc[_draw_df['interested lineages'] == 'gamma', 0]
y = _draw_df.loc[_draw_df['interested lineages'] == 'gamma', 1]
fig.add_shape(type='path',
              path=confidence_ellipse(x, y),
              line=dict(color="#78fce0",
                        width=3,
                        dash="dot", ),
              )

fig.layout.paper_bgcolor='#ffffff'
fig.layout.plot_bgcolor='#ffffff'
fig.show()
# img_bytes1 = fig.to_image(format="png", engine="kaleido")
# Image(img_bytes1)

fig.write_image('./fig3b.pdf', engine="kaleido",validate=True)


#! heatmap visualziation are provided in the `ipython_notebook/AOB viewing.ipynb` 