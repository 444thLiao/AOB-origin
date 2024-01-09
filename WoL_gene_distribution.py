# for woL phylogeny

import os
from collections import Counter, defaultdict
from glob import glob
from os.path import *
from subprocess import check_call

import pandas as pd
from api_tools.itol_func import *
from Bio import SeqIO
from tqdm import tqdm

from ..identified_genes import run_blast, parse_blast, info2style

os.chdir("/home-user/thliao/data/AOB/analysis/update_calibrations")

# 1. annotate the all nitrification-related genes
# infile = "/home-user/thliao/project/NOB/gene_phylogeny/prepared.faa"
in_proteins_dir = "/mnt/ivy/thliao/data/WoL_2019/weboflife/proteins/pergenomes"

locus2genome_WoL = {}
for faa in tqdm(glob(in_proteins_dir + "/*.faa")):
    genome = basename(faa).replace(".faa", "")
    for r in SeqIO.parse(faa, format="fasta"):
        locus = r.id.rpartition("_")[0]
        locus2genome_WoL[locus] = genome
base_odir = "/mnt/home-backup/thliao/AOB/analysis/update_calibrations"
base_in_db_dir = "/mnt/home-backup/thliao/AOB/analysis/query_db/fa/"

out_dir = join(base_odir, "Wol_annotations")

for infa in sorted(glob("/mnt/home-backup/thliao/AOB/analysis/query_db/fa/*_db.faa")):
    gene = basename(infa).split("_")[0]
    run_blast(out_dir=join(out_dir, gene), in_proteins_dir=in_proteins_dir, dbfile=infa)


# db locus2genes
locus2gene = {
    row.split("\t")[0]: row.split("\t")[1]
    for row in open("/mnt/home-backup/thliao/AOB/analysis/query_db/fa/locus2gene.list")
    .read()
    .strip("\n")
    .split("\n")
}

curated_gene2length = {}
target_genes_dict = {}
for infa in glob("/mnt/home-backup/thliao/AOB/analysis/query_db/fa/*_db.faa"):
    gene = basename(infa).split("_")[0]
    records = list(SeqIO.parse(infa, format="fasta"))
    curated_gene2length.update({_.id: len(_.seq) for _ in records})
    _list = [l for l, g in locus2gene.items() if g == gene]
    _list = [_.id for _ in records if _.id in _list]
    target_genes_dict[gene] = _list

gene2results = {}
for infa in sorted(glob("/mnt/home-backup/thliao/AOB/analysis/query_db/fa/*_db.faa")):
    gene = basename(infa).split("_")[0]
    genome2genes, locus2gene_focal_phyla = parse_blast(
        indir=join(base_odir, "Wol_annotations", gene),
        seq_file=infa,
        locus2gene=locus2gene,
        outdir=join(out_dir, gene + "_o"),
        locus2genome_mapping=locus2genome_WoL,
    )
    gene2results[gene] = genome2genes, locus2gene_focal_phyla


# 2. parse the annotation output or build gene phylogeny, again.
collect_df = []
genome2match = defaultdict(list)
collected_genes = defaultdict(list)

for gene, ref_id_lists in target_genes_dict.items():
    all_tabs = glob(join(out_dir, gene, "*.tab"))
    for tab in tqdm(all_tabs):
        _tmp = []
        name = basename(tab).replace(".tab", "")
        if not os.path.getsize(tab):
            continue
        df = pd.read_csv(tab, sep="\t", header=None)
        df.loc[:, "total_length"] = [curated_gene2length[_] for _ in df[1]]
        sub_df = df.loc[
            (df[10] <= 1e-20)
            & (df[2] >= 30)
            & (abs(df[7] - df[6]) / df["total_length"] > 0.65)
        ]
        sub_df = sub_df.sort_values(10).groupby(0).head(1)
        ref_matched_df = sub_df.loc[df[1].isin(ref_id_lists)]

        if ref_matched_df.shape[0]:
            genome2match[name].append(gene)
            _tmp += list(set(list(ref_matched_df.iloc[:, 0])))
        collected_genes[gene] += list(ref_matched_df.iloc[:, 0].unique())

# overwrite the presence of nxrA using curated list
gene_lists = (
    open(
        "/home-user/thliao/project/NOB/gene_phylogeny/newly_extracted_faa/manual_curated2.list"
    )
    .read()
    .split("\n")
)
locus2gene_nxr = defaultdict(list)
for row in gene_lists:
    if not row:
        continue
    rows = row.split("\t")
    # merge narG1 and narG2
    locus2gene_nxr[rows[0]] = "narG" if rows[2].startswith("narG") else rows[2]

_genome2nxr = defaultdict(list)
for l, g in locus2gene_nxr.items():
    if l.rpartition("_")[0] not in locus2genome_WoL:
        continue
    _genome2nxr[locus2genome_WoL[l.rpartition("_")[0]]].append(g)
    print(locus2genome_WoL[l.rpartition("_")[0]])

# map nxrB according to nxrA classification
# MAP_DICT = {'narG':'narH',
#             "nxrA":'nxrB'}
# nxrB_dict = {}
# nxrB_homologues = collected_genes['nxrB']
# for l in nxrB_homologues:
#     genome = locus2genome_WoL[l.rpartition('_')[0]]
#     nxrA_classification = list(set(_genome2nxr[genome]))
#     if not nxrA_classification:
#         continue
#     if len(nxrA_classification) == 1:
#         nxrB_dict[l] = MAP_DICT[nxrA_classification[0]]
#     else:
#         print(nxrA_classification)

for genome, matched_genes in genome2match.items():
    _l = matched_genes.copy()
    if "nxrA" in _l:
        _l.remove("nxrA")
        genome2match[genome] = _genome2nxr[genome] + _l


# --

## annotation those genes on the reference phylogeny
all_ids = [
    basename(_).replace(".faa", "") for _ in glob(join(in_proteins_dir, "*.faa"))
]
tax_tab = "/home-user/thliao/data/wol_db/data/taxonomy_df.tab"
tax_df = pd.read_csv(tax_tab, sep="\t", index_col=0)
tax_df = tax_df.reindex(all_ids)


# from bin.ncbi_convertor.toolkit import tax2tax_info
# genome2tid = {row.split('\t')[0]:row.split('\t')[1]
#               for row in open("/home-user/thliao/data/wol_db/data/taxonomy/ncbi/taxids.txt").read().split('\n') if row}
# genome2taxon_info = {genome:tax2tax_info(tid) for genome,tid in genome2tid.items()}
# tax_df = pd.DataFrame.from_dict(genome2taxon_info,orient='index')
with pd.ExcelWriter(
    join(out_dir, "presence_of_AMO.xlsx"), engine="xlsxwriter"
) as writer:
    p2num = tax_df.groupby("phylum").size()
    for gene in ["amoA", "amoB", "amoC", "cnxrA", "pnxrA", "narG"]:
        tax_df.loc[[g for g, v in genome2match.items() if gene in v], gene] = 1
        sub_df = tax_df.fillna(0)
        particualr_gene_p2num = (
            sub_df.loc[sub_df[gene] == 1, :].groupby("phylum").size()
        )
        particualr_gene_ratio = pd.DataFrame(
            (particualr_gene_p2num / p2num).dropna().sort_values()
        )
        particualr_gene_ratio.loc[:, "num"] = particualr_gene_p2num.reindex(
            particualr_gene_ratio.index
        )
        # print(gene,particualr_gene_ratio)
        particualr_gene_ratio.to_excel(writer, sheet_name=gene)


info2style["narG"] = dict(color="#f87a0d")
# text = to_color_range(new_genome2match, info2color=info2style,
#                       dataset_name="target genes", no_legend=False)

info2style = {
    "amoA": {"color": "#ff0000", "info": "amoA"},
    "amoB": {"color": "#ff0000", "info": "amoB"},
    "amoC": {"color": "#ff0000", "info": "amoC"},
    "hao": {"color": "#b68100", "info": "hao"},
    "cnxrA": {"color": "#4b85c1", "info": "nxrA"},
    "pnxrA": {"color": "#4b85c1", "info": "nxrB"},
    "narG": {"color": "#f87a0d", "info": "hao"},
}
text = to_binary_shape(
    genome2match,
    info2style=info2style,
    unfilled_other=True,
    manual_v=["amoA", "amoB", "amoC", "hao", "cnxrA", "pnxrA", "narG"],
)
with open(join(out_dir, "matched_genes_binary.txt"), "w") as f1:
    f1.write(text)
