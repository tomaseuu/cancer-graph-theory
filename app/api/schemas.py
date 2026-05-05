from pydantic import BaseModel, Field


class AnalyzeRequest(BaseModel):
    seed_genes: list[str] | None = None
    restart_probability: float = Field(default=0.3, ge=0.0, le=1.0)
    num_steps: int = Field(default=10000, gt=0)
    num_random_sets: int = Field(default=200, gt=0)
    top_n: int = Field(default=30, gt=0)


class AnalyzeResponse(BaseModel):
    run_id: str
    status: str
    message: str


class GraphSummaryResponse(BaseModel):
    num_nodes: int
    num_edges: int
    num_oncogenes: int
    oncogenes_found_in_graph: list[str]
    oncogenes_missing_from_graph: list[str]


class CandidateGene(BaseModel):
    gene_name: str
    score: float
    rank: int


class AnalysisResultResponse(BaseModel):
    run_id: str
    status: str
    seed_genes: list[str] = Field(default_factory=list)
    rwr_score: float | None = None
    p_value: float | None = None
    top_genes: list[CandidateGene] = Field(default_factory=list)
    message: str
