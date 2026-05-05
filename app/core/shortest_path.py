import random

import networkx as nx


def _build_gene_pairs(genes):
    pairs = []
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            pairs.append((genes[i], genes[j]))
    return pairs


def average_shortest_path_between_genes(graph, genes):
    pairs = _build_gene_pairs(genes)
    distances = []

    for gene1, gene2 in pairs:
        if nx.has_path(graph, gene1, gene2):
            length = nx.shortest_path_length(graph, source=gene1, target=gene2)
            distances.append(length)

    if not distances:
        return None

    return sum(distances) / len(distances)


def random_shortest_path_distribution(graph, oncogenes, num_samples=1000):
    non_oncogenes = sorted(set(graph.nodes()) - set(oncogenes))
    sample_size = len(oncogenes)
    background_scores = []

    for _ in range(num_samples):
        sample_genes = random.sample(non_oncogenes, sample_size)
        avg_distance = average_shortest_path_between_genes(graph, sample_genes)
        if avg_distance is not None:
            background_scores.append(avg_distance)

    return background_scores


def calculate_shortest_path_p_value(actual_score, random_scores):
    num_extreme = sum(score <= actual_score for score in random_scores)
    return num_extreme / len(random_scores)
