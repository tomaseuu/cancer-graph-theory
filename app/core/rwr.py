import random


def random_walk_with_restart(
    graph,
    start_node,
    restart_probability=0.3,
    num_steps=100000,
):
    visit_counts = {node: 0 for node in graph.nodes()}
    current = start_node

    for _ in range(num_steps):
        visit_counts[current] += 1
        if random.random() < restart_probability:
            current = start_node
        else:
            neighbors = list(graph.neighbors(current))
            if neighbors:
                weights = [graph[current][neighbor].get("weight", 1) for neighbor in neighbors]
                current = random.choices(neighbors, weights=weights)[0]

    total_visits = sum(visit_counts.values())
    stationary = {
        node: count / total_visits
        for node, count in visit_counts.items()
    }
    return stationary


def _filter_seed_genes(graph, seed_genes):
    filtered_seed_genes = []
    for gene in seed_genes:
        if gene in graph.nodes():
            filtered_seed_genes.append(gene)
    return filtered_seed_genes


def _average_proximity_for_gene_set(
    graph,
    genes,
    restart_probability=0.3,
    num_steps=100000,
):
    proximity_scores = []

    for gene in genes:
        stationary = random_walk_with_restart(
            graph,
            gene,
            restart_probability=restart_probability,
            num_steps=num_steps,
        )
        other_genes = []
        for other_gene in genes:
            if other_gene != gene:
                other_genes.append(other_gene)

        total = 0
        for other_gene in other_genes:
            total += stationary.get(other_gene, 0)

        avg_frequency = total / len(other_genes)
        proximity_scores.append(avg_frequency)

    return sum(proximity_scores) / len(proximity_scores)


def calculate_rwr_proximity(
    graph,
    seed_genes,
    restart_probability=0.3,
    num_steps=100000,
):
    filtered_seed_genes = _filter_seed_genes(graph, seed_genes)
    return _average_proximity_for_gene_set(
        graph,
        filtered_seed_genes,
        restart_probability=restart_probability,
        num_steps=num_steps,
    )


def random_rwr_distribution(
    graph,
    seed_genes,
    num_samples=1000,
    restart_probability=0.3,
    num_steps=100000,
):
    filtered_seed_genes = _filter_seed_genes(graph, seed_genes)
    sample_size = len(filtered_seed_genes)
    background_scores = []
    non_oncogenes = sorted(set(graph.nodes()) - set(filtered_seed_genes))

    for _ in range(num_samples):
        sample_genes = random.sample(non_oncogenes, sample_size)
        score = _average_proximity_for_gene_set(
            graph,
            sample_genes,
            restart_probability=restart_probability,
            num_steps=num_steps,
        )
        background_scores.append(score)

    return background_scores


def calculate_rwr_p_value(actual_score, random_scores):
    num_extreme = sum(score >= actual_score for score in random_scores)
    return num_extreme / len(random_scores)
