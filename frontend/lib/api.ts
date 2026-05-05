const API_BASE_URL = (process.env.NEXT_PUBLIC_API_BASE_URL || "http://localhost:8000").replace(
  /\/+$/,
  "",
);


export type AnalyzePayload = {
  seed_genes?: string[];
  restart_probability: number;
  num_steps: number;
  num_random_sets: number;
  top_n: number;
  use_ml_ranking: boolean;
};

export type AnalyzeResponse = {
  run_id: string;
  status: string;
  message: string;
};

export type StatusResponse = {
  run_id: string;
  status: string;
  message: string;
};

export type CandidateGene = {
  gene_name: string;
  score: number;
  rank: number;
  rwr_score?: number | null;
  ml_probability?: number | null;
  final_score?: number | null;
  degree?: number | null;
  pagerank?: number | null;
  betweenness_centrality?: number | null;
  closeness_centrality?: number | null;
  shortest_path_to_nearest_seed?: number | null;
  oncogene_neighbor_count?: number | null;
};

export type AnalysisResult = {
  run_id: string;
  status: string;
  seed_genes: string[];
  rwr_score?: number | null;
  p_value?: number | null;
  top_genes: CandidateGene[];
  message: string;
  error_message?: string | null;
  created_at?: string | null;
  completed_at?: string | null;
};


async function requestJson<T>(path: string, init?: RequestInit): Promise<T> {
  let response: Response;

  try {
    response = await fetch(`${API_BASE_URL}${path}`, {
      ...init,
      headers: {
        "Content-Type": "application/json",
        ...(init?.headers ?? {}),
      },
      cache: "no-store",
    });
  } catch {
    throw new Error(
      "Could not connect to the OncoGraph API. Make sure the FastAPI server is running on http://localhost:8000.",
    );
  }

  if (!response.ok) {
    let message = "Request failed.";
    try {
      const body = await response.json();
      message = body.detail || body.message || message;
    } catch {
      const text = await response.text();
      if (text) {
        message = text;
      }
    }
    throw new Error(message);
  }

  return response.json() as Promise<T>;
}


export async function startAnalysis(payload: AnalyzePayload) {
  return requestJson<AnalyzeResponse>("/analyze", {
    method: "POST",
    body: JSON.stringify(payload),
  });
}


export async function getStatus(runId: string) {
  return requestJson<StatusResponse>(`/status/${runId}`);
}


export async function getResults(runId: string) {
  return requestJson<AnalysisResult>(`/results/${runId}`);
}
