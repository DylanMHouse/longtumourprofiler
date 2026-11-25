#!/usr/bin/env python3
# Adapted from workflow/scripts/plot_phyclone_tree.py in Roth-Lab/PhyClone-Workflow
# (https://github.com/Roth-Lab/PhyClone-Workflow, commit 088355457813b6adc3192e6b81c51ad2463a9208),
# which is licensed under the GNU General Public License v3.0. This file stays under GPL-3.0;
# the licence text is in assets/licenses/GPL-3.0.txt. The rest of the pipeline is MIT licensed.
#
# Modified 2026-11-25: the Snakemake entry point is replaced by command-line arguments, and
# the sys import it needed is removed.

from skbio.tree import TreeNode
import networkx as nx
import pandas as pd
from seaborn import color_palette
import colorcet


def draw_tree(t, filename):
    draw_graph = nx.nx_agraph.to_agraph(t)

    draw_graph.graph_attr["bgcolor"] = "transparent"
    draw_graph.graph_attr["fontname"] = "Nimbus Sans"
    draw_graph.node_attr["fontname"] = "Nimbus Sans"

    draw_graph.draw(filename, prog="dot")


def scikit_tree_to_networkx(skbio_tree, clone_sizes, clone_colour_dict, outlier_name="-1"):
    nx_graph = nx.DiGraph()

    grey_to_use = "#696969"
    node_label_template = "{}\nNum SNV: {}"

    for node in skbio_tree.traverse():

        node_name = node.name

        if node.is_tip():
            nx_graph.add_node(node_name)

        for child in node.children:
            child_name = child.name
            nx_graph.add_edge(node_name, child_name)

        if node_name in clone_sizes:
            node_colour = clone_colour_dict[node_name]
            fillcolor = node_colour + "1A"
            node_attr_dict = {
                "label": node_label_template.format(node_name, clone_sizes[node_name]),
                "color": node_colour,
                "fillcolor": fillcolor,
                "style": "filled",
                "penwidth": 3,
            }
            nx_graph.nodes[node_name].update(node_attr_dict)
        else:
            node_colour = grey_to_use
            fillcolor = node_colour + "1A"
            node_attr_dict = {
                "color": node_colour,
                "fillcolor": fillcolor,
                "style": "filled",
                "penwidth": 3,
            }
            nx_graph.nodes[node_name].update(node_attr_dict)

    if outlier_name in clone_colour_dict:
        nx_graph.add_node(outlier_name)
        node_colour = clone_colour_dict[outlier_name]
        fillcolor = node_colour + "1A"
        node_attr_dict = {
            "label": node_label_template.format(outlier_name, clone_sizes[outlier_name]),
            "color": node_colour,
            "fillcolor": fillcolor,
            "style": "filled, dashed",
            "penwidth": 3,
        }
        nx_graph.nodes[outlier_name].update(node_attr_dict)

    return nx_graph


def create_clone_colour_dict(node_names):
    num_colours = len(node_names)
    colour_list = color_palette(colorcet.glasbey_bw_minc_20, num_colours).as_hex()
    gt_colour_dict = dict(zip(node_names, colour_list))
    return gt_colour_dict


def main(newick_file, results_table, outfile):
    clone_df = pd.read_table(results_table)
    clone_df = clone_df[["clone_id", "mutation_id"]].drop_duplicates()
    clone_df["clone_id"] = clone_df["clone_id"].astype(str)
    clone_sizes = clone_df["clone_id"].value_counts().to_dict()
    sk_tree = TreeNode.read(newick_file, convert_underscores=False)
    clone_colour_dict = create_clone_colour_dict(clone_df["clone_id"].unique())
    nx_tree = scikit_tree_to_networkx(sk_tree, clone_sizes, clone_colour_dict)
    draw_tree(nx_tree, outfile)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Draw the PhyClone MAP tree with one node per clone")
    parser.add_argument("--newick", required=True, help="MAP tree from phyclone map, in newick format")
    parser.add_argument("--results-table", required=True, help="MAP results table from phyclone map")
    parser.add_argument("--output", required=True, help="Image to write; the extension sets the format")
    args = parser.parse_args()

    main(args.newick, args.results_table, args.output)
