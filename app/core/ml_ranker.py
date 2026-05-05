import random

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def normalize_scores(values):
    values_array = np.asarray(values, dtype=float)
    if values_array.size == 0:
        return []

    min_value = float(values_array.min())
    max_value = float(values_array.max())
    if max_value == min_value:
        return [0.0 for _ in values_array]

    normalized = (values_array - min_value) / (max_value - min_value)
    return normalized.tolist()


def compute_graph_features(graph, seed_genes, rwr_scores=None, random_state=42):
    filtered_seed_genes = [gene for gene in seed_genes if gene in graph.nodes()]
    seed_gene_set = set(filtered_seed_genes)
    num_nodes = graph.number_of_nodes()

    pagerank_scores = nx.pagerank(graph)
    betweenness_scores = nx.betweenness_centrality(
        graph,
        k=min(500, num_nodes),
        seed=random_state,
    )
    closeness_scores = nx.closeness_centrality(graph)

    if filtered_seed_genes:
        shortest_paths = nx.multi_source_dijkstra_path_length(graph, filtered_seed_genes)
    else:
        shortest_paths = {}

    rows = []
    for gene_name in graph.nodes():
        neighbors = set(graph.neighbors(gene_name))
        rows.append(
            {
                "gene_name": gene_name,
                "rwr_score": float((rwr_scores or {}).get(gene_name, 0.0)),
                "degree": int(graph.degree(gene_name)),
                "pagerank": float(pagerank_scores.get(gene_name, 0.0)),
                "betweenness_centrality": float(betweenness_scores.get(gene_name, 0.0)),
                "closeness_centrality": float(closeness_scores.get(gene_name, 0.0)),
                "shortest_path_to_nearest_seed": float(shortest_paths.get(gene_name, 999)),
                "oncogene_neighbor_count": int(len(neighbors & seed_gene_set)),
            }
        )

    return pd.DataFrame(rows)


def build_training_dataset(features_df, seed_genes, negative_sample_size=None, random_state=42):
    filtered_seed_genes = [gene for gene in seed_genes if gene in set(features_df["gene_name"])]
    if not filtered_seed_genes:
        raise ValueError("No seed genes were found in the graph for ML training.")

    non_seed_genes = [
        gene_name
        for gene_name in features_df["gene_name"].tolist()
        if gene_name not in filtered_seed_genes
    ]
    if not non_seed_genes:
        raise ValueError("No non-seed genes are available for ML training.")

    if negative_sample_size is None:
        negative_sample_size = min(len(filtered_seed_genes) * 5, len(non_seed_genes))
    negative_sample_size = min(negative_sample_size, len(non_seed_genes))

    random_generator = random.Random(random_state)
    negative_genes = random_generator.sample(non_seed_genes, negative_sample_size)
    if not negative_genes:
        raise ValueError("No negative training genes could be sampled.")

    training_df = features_df[features_df["gene_name"].isin(filtered_seed_genes + negative_genes)].copy()
    training_df["label"] = training_df["gene_name"].isin(filtered_seed_genes).astype(int)
    return training_df


def train_gene_ranker(training_df):
    feature_columns = [
        "rwr_score",
        "degree",
        "pagerank",
        "betweenness_centrality",
        "closeness_centrality",
        "shortest_path_to_nearest_seed",
        "oncogene_neighbor_count",
    ]

    model = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            ("classifier", LogisticRegression(max_iter=1000)),
        ]
    )
    model.fit(training_df[feature_columns], training_df["label"])
    return model, feature_columns


def rank_genes_with_ml(graph, seed_genes, rwr_scores=None, top_n=30, random_state=42):
    filtered_seed_genes = [gene for gene in seed_genes if gene in graph.nodes()]
    features_df = compute_graph_features(
        graph,
        filtered_seed_genes,
        rwr_scores=rwr_scores,
        random_state=random_state,
    )
    training_df = build_training_dataset(
        features_df,
        filtered_seed_genes,
        random_state=random_state,
    )
    model, feature_columns = train_gene_ranker(training_df)

    features_df["ml_probability"] = model.predict_proba(features_df[feature_columns])[:, 1]
    features_df["normalized_rwr_score"] = normalize_scores(features_df["rwr_score"].tolist())
    features_df["final_score"] = (
        0.7 * features_df["normalized_rwr_score"] + 0.3 * features_df["ml_probability"]
    )

    candidate_df = features_df[~features_df["gene_name"].isin(filtered_seed_genes)].copy()
    candidate_df = candidate_df.sort_values("final_score", ascending=False).head(top_n)

    candidates = []
    for rank, row in enumerate(candidate_df.itertuples(index=False), start=1):
        candidates.append(
            {
                "gene_name": row.gene_name,
                "score": float(row.final_score),
                "rwr_score": float(row.rwr_score),
                "ml_probability": float(row.ml_probability),
                "final_score": float(row.final_score),
                "rank": rank,
                "degree": int(row.degree),
                "pagerank": float(row.pagerank),
                "betweenness_centrality": float(row.betweenness_centrality),
                "closeness_centrality": float(row.closeness_centrality),
                "shortest_path_to_nearest_seed": float(row.shortest_path_to_nearest_seed),
                "oncogene_neighbor_count": int(row.oncogene_neighbor_count),
            }
        )

    return candidates
