from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from app.core.gene_ranking import rank_candidate_genes, save_candidate_genes
from app.core.graph_loader import load_graph_and_oncogenes
from app.core.randomization import run_randomized_graph_rwr_analysis
from app.core.rwr import calculate_rwr_p_value, calculate_rwr_proximity, random_rwr_distribution
from app.core.shortest_path import (
    average_shortest_path_between_genes,
    calculate_shortest_path_p_value,
    random_shortest_path_distribution,
)
from app.utils.helpers import ensure_output_dir, set_random_seed


DATA_DIR = Path("data")
OUTPUT_FILE = Path("outputs/Thomas_Le_onco_predictions.txt")


def plot_shortest_path_distribution(actual_score, random_scores):
    plt.figure(figsize=(10, 6))
    plt.hist(random_scores, bins=30, color="lightblue", edgecolor="black")
    plt.axvline(
        actual_score,
        color="red",
        linestyle="dashed",
        linewidth=2,
        label="Oncogene Score",
    )
    plt.title("Distribution of Average Shortest Path (Random Sets)")
    plt.xlabel("Average Shortest Path")
    plt.ylabel("Frequency")
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_rwr_distribution(actual_score, random_scores, title):
    x_max = max(np.percentile(random_scores, 95), actual_score * 1.05)
    plt.figure(figsize=(10, 6))
    plt.hist(random_scores, bins=30, color="lightblue", edgecolor="black")
    plt.axvline(
        actual_score,
        color="red",
        linestyle="dashed",
        linewidth=2,
        label="Oncogene RWR Score",
    )
    plt.xlim(0, x_max)
    plt.title(title)
    plt.xlabel("Average RWR Proximity Score")
    plt.ylabel("Frequency")
    plt.legend()
    plt.grid(True)
    plt.show()


def main():
    set_random_seed(42)
    ensure_output_dir(OUTPUT_FILE)

    print("Loading graph and oncogene data...")
    graph, oncogenes = load_graph_and_oncogenes(
        DATA_DIR / "onco_genes.txt",
        DATA_DIR / "interacting_proteins.txt",
    )
    print(
        f"Loaded {graph.number_of_nodes()} proteins, "
        f"{graph.number_of_edges()} interactions, and {len(oncogenes)} oncogenes."
    )

    print("\nRunning shortest path analysis...")
    avg_onco_distance = average_shortest_path_between_genes(graph, oncogenes)
    if avg_onco_distance is None:
        raise ValueError("No connected oncogene pairs were found.")
    background_shortest_paths = random_shortest_path_distribution(graph, oncogenes)
    shortest_path_p_value = calculate_shortest_path_p_value(
        avg_onco_distance,
        background_shortest_paths,
    )
    print(f"Average shortest path between tumorigenic genes: {avg_onco_distance:.3f}")
    print(f"Shortest path p-value: {shortest_path_p_value:.3f}")
    plot_shortest_path_distribution(avg_onco_distance, background_shortest_paths)

    print("\nRunning Random Walk with Restart analysis...")
    oncogene_rwr_score = calculate_rwr_proximity(graph, oncogenes)
    background_rwr_scores = random_rwr_distribution(graph, oncogenes)
    rwr_p_value = calculate_rwr_p_value(oncogene_rwr_score, background_rwr_scores)
    print(f"RWR proximity score for oncogenes: {oncogene_rwr_score:.6f}")
    print(f"RWR p-value: {rwr_p_value:.4f}")
    plot_rwr_distribution(
        oncogene_rwr_score,
        background_rwr_scores,
        "RWR Proximity Scores (Background vs Oncogenes)",
    )

    print("\nRanking candidate genes...")
    top_candidates = rank_candidate_genes(graph, oncogenes, top_n=30)
    save_candidate_genes(top_candidates, OUTPUT_FILE)
    print(f"Top 30 novel tumorigenic genes written to {OUTPUT_FILE}")

    print("\nRunning randomized graph RWR analysis...")
    randomized_results = run_randomized_graph_rwr_analysis(graph, oncogenes)
    print(
        "Randomized graph RWR proximity score for oncogenes: "
        f"{randomized_results['rwr_score']:.6f}"
    )
    print(f"RWR p-value (randomized graph): {randomized_results['p_value']:.4f}")
    plot_rwr_distribution(
        randomized_results["rwr_score"],
        randomized_results["background_scores"],
        "RWR Proximity Scores (Randomized Graph)",
    )

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
