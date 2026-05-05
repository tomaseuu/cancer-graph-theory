from datetime import datetime

from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Text
from sqlalchemy.orm import relationship

from app.db.database import Base


class AnalysisRun(Base):
    __tablename__ = "analysis_runs"

    id = Column(String, primary_key=True, index=True)
    status = Column(String, nullable=False)
    seed_genes = Column(Text, nullable=False)
    restart_probability = Column(Float, nullable=False)
    num_steps = Column(Integer, nullable=False)
    num_random_sets = Column(Integer, nullable=False)
    rwr_score = Column(Float, nullable=True)
    p_value = Column(Float, nullable=True)
    message = Column(Text, nullable=False)
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    completed_at = Column(DateTime, nullable=True)

    candidate_genes = relationship(
        "CandidateGene",
        back_populates="analysis_run",
        cascade="all, delete-orphan",
        order_by="CandidateGene.rank",
    )


class CandidateGene(Base):
    __tablename__ = "candidate_genes"

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(String, ForeignKey("analysis_runs.id"), nullable=False, index=True)
    gene_name = Column(String, nullable=False)
    score = Column(Float, nullable=False)
    rank = Column(Integer, nullable=False)
    rwr_score = Column(Float, nullable=True)
    ml_probability = Column(Float, nullable=True)
    final_score = Column(Float, nullable=True)
    degree = Column(Integer, nullable=True)
    pagerank = Column(Float, nullable=True)
    betweenness_centrality = Column(Float, nullable=True)
    closeness_centrality = Column(Float, nullable=True)
    shortest_path_to_nearest_seed = Column(Float, nullable=True)
    oncogene_neighbor_count = Column(Integer, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    analysis_run = relationship("AnalysisRun", back_populates="candidate_genes")
