from pathlib import Path

import networkx as nx


def load_oncogenes(file_path):
    oncogenes = []
    with Path(file_path).open("r", encoding="utf-8-sig") as file:
        for line in file:
            gene = line.strip()
            if gene:
                oncogenes.append(gene.upper())
    return oncogenes


def load_interactions(file_path):
    interactions = []
    with Path(file_path).open("r", encoding="utf-8-sig") as file:
        for line in file:
            protein1, protein2 = line.strip().split()
            interactions.append((protein1.strip().upper(), protein2.strip().upper()))
    return interactions


def build_graph(interactions):
    graph = nx.Graph()
    graph.add_edges_from(interactions)
    return graph


def load_graph_and_oncogenes(onco_path, interactions_path):
    oncogenes = load_oncogenes(onco_path)
    interactions = load_interactions(interactions_path)
    graph = build_graph(interactions)
    return graph, oncogenes
