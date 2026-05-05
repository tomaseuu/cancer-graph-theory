from pathlib import Path

from app.core.ml_ranker import rank_genes_with_ml
from app.core.rwr import random_walk_with_restart
from app.utils.helpers import ensure_output_dir


def compute_average_rwr_scores(
    graph,
    seed_genes,
    restart_probability=0.3,
    num_steps=100000,
):
    filtered_seed_genes = []
    for gene in seed_genes:
        if gene in graph.nodes():
            filtered_seed_genes.append(gene)

    combined_frequencies = {node: 0 for node in graph.nodes()}

    for gene in filtered_seed_genes:
        stationary = random_walk_with_restart(
            graph,
            gene,
            restart_probability=restart_probability,
            num_steps=num_steps,
        )
        for node, frequency in stationary.items():
            combined_frequencies[node] += frequency

    average_frequencies = {}
    for node, frequency in combined_frequencies.items():
        average_frequencies[node] = frequency / len(filtered_seed_genes)

    return average_frequencies, filtered_seed_genes


def rank_candidate_genes(
    graph,
    seed_genes,
    restart_probability=0.3,
    num_steps=100000,
    top_n=30,
):
    average_frequencies, filtered_seed_genes = compute_average_rwr_scores(
        graph,
        seed_genes,
        restart_probability=restart_probability,
        num_steps=num_steps,
    )

    candidate_scores = {}
    for gene, score in average_frequencies.items():
        if gene not in filtered_seed_genes:
            candidate_scores[gene] = score

    ranked_genes = sorted(
        candidate_scores.items(),
        key=lambda item: item[1],
        reverse=True,
    )
    return ranked_genes[:top_n]


def rank_candidate_genes_with_ml(
    graph,
    seed_genes,
    restart_probability=0.3,
    num_steps=100000,
    top_n=30,
    random_state=42,
):
    average_frequencies, filtered_seed_genes = compute_average_rwr_scores(
        graph,
        seed_genes,
        restart_probability=restart_probability,
        num_steps=num_steps,
    )
    return rank_genes_with_ml(
        graph,
        filtered_seed_genes,
        rwr_scores=average_frequencies,
        top_n=top_n,
        random_state=random_state,
    )


def save_candidate_genes(candidates, output_path):
    ensure_output_dir(output_path)

    with Path(output_path).open("w", encoding="utf-8") as file:
        for candidate in candidates:
            if isinstance(candidate, tuple):
                gene_name = candidate[0]
            elif isinstance(candidate, dict):
                gene_name = candidate["gene_name"]
            else:
                gene_name = candidate
            file.write(f"{gene_name}\n")
