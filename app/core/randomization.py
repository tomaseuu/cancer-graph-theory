import random

import networkx as nx

from app.core.rwr import random_walk_with_restart


def randomize_graph_preserve_degree(graph, num_swaps=None):
    randomized_graph = graph.copy()

    if num_swaps is None:
        num_swaps = 10 * randomized_graph.number_of_edges()

    max_tries = 100 * randomized_graph.number_of_edges()
    nx.double_edge_swap(randomized_graph, nswap=num_swaps, max_tries=max_tries)
    return randomized_graph


def run_randomized_graph_rwr_analysis(
    graph,
    seed_genes,
    restart_probability=0.3,
    num_steps=100000,
    num_samples=1000,
):
    randomized_graph = randomize_graph_preserve_degree(graph)

    random_oncogenes = []
    for gene in seed_genes:
        if gene in randomized_graph.nodes:
            random_oncogenes.append(gene)

    random_proximity_scores = []

    for gene in random_oncogenes:
        stationary = random_walk_with_restart(
            randomized_graph,
            gene,
            restart_probability=restart_probability,
            num_steps=num_steps,
        )
        others = []
        for other_gene in random_oncogenes:
            if other_gene != gene:
                others.append(other_gene)

        total = 0
        for other_gene in others:
            if other_gene in stationary:
                total += stationary[other_gene]
            else:
                total += 0

        average_frequency = total / len(others)
        random_proximity_scores.append(average_frequency)

    random_onco_rwr_score = sum(random_proximity_scores) / len(random_proximity_scores)

    sample_size = len(random_oncogenes)
    random_background_rwr_scores = []
    random_non_oncogenes = sorted(set(randomized_graph.nodes()) - set(random_oncogenes))

    for _ in range(num_samples):
        sample_genes = random.sample(random_non_oncogenes, sample_size)
        proximity_scores = []

        for gene in sample_genes:
            stationary = random_walk_with_restart(
                randomized_graph,
                gene,
                restart_probability=restart_probability,
                num_steps=num_steps,
            )
            others = []
            for other_gene in random_oncogenes:
                if other_gene != gene:
                    others.append(other_gene)

            total = 0
            for other_gene in others:
                if other_gene in stationary:
                    total += stationary[other_gene]
                else:
                    total += 0

            average_frequency = total / len(others)
            proximity_scores.append(average_frequency)

        random_background_rwr_scores.append(sum(proximity_scores) / len(proximity_scores))

    num_extreme = sum(score >= random_onco_rwr_score for score in random_background_rwr_scores)
    rand_p_value = num_extreme / len(random_background_rwr_scores)

    return {
        "graph": randomized_graph,
        "seed_genes": random_oncogenes,
        "rwr_score": random_onco_rwr_score,
        "background_scores": random_background_rwr_scores,
        "p_value": rand_p_value,
    }
