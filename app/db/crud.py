import json
from datetime import datetime

from app.api.schemas import AnalysisResultResponse, CandidateGene as CandidateGeneResponse
from app.db.models import AnalysisRun, CandidateGene


def create_analysis_run(db, run_id, request, seed_genes):
    analysis_run = AnalysisRun(
        id=run_id,
        status="running",
        seed_genes=json.dumps(seed_genes),
        restart_probability=request.restart_probability,
        num_steps=request.num_steps,
        num_random_sets=request.num_random_sets,
        rwr_score=None,
        p_value=None,
        message="Analysis is running",
        error_message=None,
        created_at=datetime.utcnow(),
        completed_at=None,
    )
    db.add(analysis_run)
    db.commit()
    db.refresh(analysis_run)
    return analysis_run


def update_analysis_run_status(db, run_id, status, message=None, error_message=None):
    analysis_run = get_analysis_run(db, run_id)
    if analysis_run is None:
        return None

    analysis_run.status = status
    if message is not None:
        analysis_run.message = message
    analysis_run.error_message = error_message
    if status in {"completed", "failed"}:
        analysis_run.completed_at = datetime.utcnow()

    db.commit()
    db.refresh(analysis_run)
    return analysis_run


def complete_analysis_run(db, run_id, rwr_score, p_value, top_genes):
    analysis_run = get_analysis_run(db, run_id)
    if analysis_run is None:
        return None

    analysis_run.status = "completed"
    analysis_run.rwr_score = rwr_score
    analysis_run.p_value = p_value
    analysis_run.message = "Analysis completed successfully"
    analysis_run.error_message = None
    analysis_run.completed_at = datetime.utcnow()
    analysis_run.candidate_genes.clear()

    for gene in top_genes:
        analysis_run.candidate_genes.append(
            CandidateGene(
                gene_name=gene["gene_name"],
                score=gene["score"],
                rank=gene["rank"],
                created_at=datetime.utcnow(),
            )
        )

    db.commit()
    db.refresh(analysis_run)
    return analysis_run


def get_analysis_run(db, run_id):
    return db.query(AnalysisRun).filter(AnalysisRun.id == run_id).first()


def get_top_candidate_genes(db, run_id):
    return (
        db.query(CandidateGene)
        .filter(CandidateGene.run_id == run_id)
        .order_by(CandidateGene.rank.asc())
        .all()
    )


def analysis_run_to_response(run):
    seed_genes = json.loads(run.seed_genes) if run.seed_genes else []
    top_genes = [candidate_gene_to_response(candidate) for candidate in run.candidate_genes]

    return AnalysisResultResponse(
        run_id=run.id,
        status=run.status,
        seed_genes=seed_genes,
        rwr_score=run.rwr_score,
        p_value=run.p_value,
        top_genes=top_genes,
        message=run.message,
        error_message=run.error_message,
        created_at=run.created_at.isoformat() if run.created_at else None,
        completed_at=run.completed_at.isoformat() if run.completed_at else None,
    )


def candidate_gene_to_response(candidate):
    return CandidateGeneResponse(
        gene_name=candidate.gene_name,
        score=candidate.score,
        rank=candidate.rank,
    )
