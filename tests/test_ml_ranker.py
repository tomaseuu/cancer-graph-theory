import unittest

import networkx as nx

from app.core.ml_ranker import compute_graph_features, normalize_scores, rank_genes_with_ml


class MlRankerTests(unittest.TestCase):
    def setUp(self):
        self.graph = nx.Graph()
        self.graph.add_edges_from(
            [
                ("A", "B"),
                ("A", "C"),
                ("B", "D"),
                ("C", "D"),
                ("D", "E"),
                ("E", "F"),
            ]
        )
        self.seed_genes = ["A", "B"]

    def test_normalize_scores(self):
        self.assertEqual(normalize_scores([2.0, 4.0, 6.0]), [0.0, 0.5, 1.0])
        self.assertEqual(normalize_scores([3.0, 3.0]), [0.0, 0.0])

    def test_compute_graph_features(self):
        features_df = compute_graph_features(
            self.graph,
            self.seed_genes,
            rwr_scores={"A": 0.4, "C": 0.2},
        )

        self.assertEqual(len(features_df), self.graph.number_of_nodes())
        row_c = features_df[features_df["gene_name"] == "C"].iloc[0]
        self.assertEqual(row_c["degree"], 2)
        self.assertEqual(row_c["shortest_path_to_nearest_seed"], 1.0)
        self.assertEqual(row_c["oncogene_neighbor_count"], 1)
        self.assertEqual(row_c["rwr_score"], 0.2)

    def test_rank_genes_with_ml_excludes_seed_genes(self):
        candidates = rank_genes_with_ml(
            self.graph,
            self.seed_genes,
            rwr_scores={"A": 0.4, "B": 0.35, "C": 0.2, "D": 0.1},
            top_n=3,
            random_state=42,
        )

        candidate_names = [candidate["gene_name"] for candidate in candidates]
        self.assertNotIn("A", candidate_names)
        self.assertNotIn("B", candidate_names)
        self.assertTrue(all("ml_probability" in candidate for candidate in candidates))


if __name__ == "__main__":
    unittest.main()
