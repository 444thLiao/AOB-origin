############### delta-statstic

import sys
sys.path.insert(0, "/home-user/thliao/ref_github/delta-statistic/Delta-Python")
# Imports
import numpy as np
import os
import pandas as pd

# Delta functions
from delta_functs import delta
from pastml.tree import read_tree, name_tree
from pastml.acr import acr
from pastml.annotation import preannotate_forest
from pastml import col_name2cat
from collections import defaultdict, Counter

def _validate_input(tree_nwk, data, data_sep=",", single_tree_file=False):
    """tree_nwk      = Represents the path to the Newick file containing the tree or a string with the tree itself.
    data             = Represents the path to the data file or DataFrame used for annotation with leaf states.
    data_sep         (default: ',')   = Separator used in the data file.
    single_tree_file (default: False) = Boolean value that specifies whether the input tree is provided as a single file."""

    if single_tree_file == False:
        with open(
            tree_nwk, "r"
        ) as f:  # Reads the tree from a Newick file and returns its roots
            nwks = f.read().replace("\n", "")
        roots = [read_tree(tree_nwk)]
    else:
        roots = [read_tree(tree_nwk)]  # Reads the newick tree and returns its roots

    column2annotated = (
        Counter()
    )  # Counter to keep track of the number of times each column is annotated
    column2states = defaultdict(
        set
    )  # Dictionary to store the unique states for each column

    # Read the data as a pandas DataFrame
    if type(data) is pd.DataFrame:
        df = data
    else:
        df = pd.read_csv(data, sep=data_sep, index_col=0, header=0, dtype=str)
    df.index = df.index.map(str)
    df.columns = [col_name2cat(column) for column in df.columns]
    columns = df.columns

    node_names = set.union(
        *[{n.name for n in root.traverse() if n.name} for root in roots]
    )  # Get the names of the nodes in the tree
    df_index_names = set(df.index)  # Get the index names from the DataFrame
    common_ids = list(
        node_names & df_index_names
    )  # Find the common IDs between node names and DataFrame index names

    # strip quotes if needed
    if not common_ids:
        node_names = {_.strip("'").strip('"') for _ in node_names}
        common_ids = node_names & df_index_names
        if common_ids:
            for root in roots:
                for n in root.traverse():
                    n.name = n.name.strip("'").strip('"')

    # Preannotate the forest with the DataFrame
    preannotate_forest(roots, df=df)

    # Populate the column2states dictionary with unique states for each column
    for c in df.columns:
        column2states[c] |= {_ for _ in df[c].unique() if pd.notnull(_) and _ != ""}

    num_tips = 0

    # Count the number of annotated columns for each node
    column2annotated_states = defaultdict(set)
    for root in roots:
        for n in root.traverse():
            for c in columns:
                vs = getattr(n, c, set())
                column2states[c] |= vs
                column2annotated_states[c] |= vs
                if vs:
                    column2annotated[c] += 1
            if n.is_leaf():
                num_tips += 1

    if column2annotated:
        c, num_annotated = min(column2annotated.items(), key=lambda _: _[1])
    else:
        c, num_annotated = columns[0], 0

    # Calculate the percentage of unknown tip annotations
    percentage_unknown = (num_tips - num_annotated) / num_tips
    if percentage_unknown >= 0.9:
        raise ValueError(
            '{:.1f}% of tip annotations for character "{}" are unknown, '
            "not enough data to infer ancestral states. "
            "{}".format(
                percentage_unknown * 100,
                c,
                "Check your annotation file and if its ids correspond to the tree tip/node names."
                if data
                else "You tree file should contain character state annotations, "
                "otherwise consider specifying a metadata file.",
            )
        )
    c, states = min(column2annotated_states.items(), key=lambda _: len(_[1]))

    # Check if the number of unique states is too high for the given number of tips
    if len(states) > num_tips * 0.75:
        raise ValueError(
            'Character "{}" has {} unique states annotated in this tree: {}, '
            "which is too much to infer on a {} with only {} tips. "
            "Make sure the character you are analysing is discrete, and if yes use a larger tree.".format(
                c,
                len(states),
                states,
                "tree" if len(roots) == 1 else "forest",
                num_tips,
            )
        )

    # Convert column2states to numpy arrays and sort the states
    column2states = {c: np.array(sorted(states)) for c, states in column2states.items()}

    # Name the trees in the forest
    for i, tree in enumerate(roots):
        name_tree(tree, suffix="" if len(roots) == 1 else "_{}".format(i))

    return roots, columns, column2states

def marginal(
    tree, data, prediction_method="MPPA", model="F81", threads=0, single_tree_file=False
):
    """tree           = Represents the path to the Newick file containing the tree or a string with the tree itself.
    data              = Represents the path to the data file or DataFrame used for annotation with leaf states.
    prediction_method (default: 'MPPA') = Specifies the ancestral character prediction method.
    model             (default: 'F81')  = Specifies the evolutionary model used for reconstruction.
    threads           (default: 0)      = Specifies the number of threads to use for the analysis.
    single_tree_file  (default: False)  = Boolean value that specifies whether the input tree is provided as a single file."""

    # Set the number of threads based on the available CPU cores
    if threads < 1:
        threads = max(os.cpu_count(), 1)

    # Validate the input and get the roots, columns, and column2states
    roots, columns, column2states = _validate_input(
        tree_nwk=tree, data=data, data_sep=",", single_tree_file=single_tree_file
    )
    #print("perform acr")
    # Perform the ancestral character reconstruction (ACR) analysis
    acr_results = acr(
        forest=roots,
        columns=columns,
        column2states=column2states,
        prediction_method=prediction_method,
        model=model,
        threads=threads,
    )

    # Get the leaf names from the tree
    leaf_names = read_tree(tree).get_leaf_names()

    # Get the marginal probabilities and exclude the leaf nodes
    marginal = np.asarray(acr_results[0]["marginal_probabilities"].drop(leaf_names))

    return marginal

# Delta Inputs
lambda0 = 0.1  # rate parameter of the proposal
se = 0.5  # standard deviation of the proposal
sim = 100000  # number of iterations
thin = 10  # Keep only each xth iterate
burn = 100  # Burned-in iterates
ent_type = "LSE"  # Linear Shannon Entropy


p = "/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/comparisons/phylosig/test.d"
df = pd.read_csv(p, sep=",", index_col=0, header=0, dtype=str)

m = marginal("/home-user/thliao/project/AOB/analysis/20210713_repeat/comparisons/phylosig/input.tre",df, single_tree_file=False,)
max_Delta_Final = delta(x=m, lambda0=lambda0, se=se, sim=sim, burn=burn, thin=thin, ent_type="LSE")

from tqdm import tqdm
import random
_c = []
for _ in tqdm(range(100)):
    _df = df.copy()
    idx = list(df.index)
    random.shuffle(idx)
    _df.index = idx
    m = marginal("/home-user/thliao/project/AOB/analysis/20210713_repeat/comparisons/phylosig/input.tre",_df, single_tree_file=False,)
    delta_random = delta(x=m, lambda0=lambda0, se=se, sim=sim, burn=burn, thin=thin, ent_type="LSE")
    _c.append(delta_random)



