from datetime import datetime
from pathlib import Path
from threading import Lock
from uuid import uuid4

from fastapi import APIRouter, BackgroundTasks, HTTPException

from app.api.schemas import (
    AnalysisResultResponse,
    AnalyzeRequest,
    AnalyzeResponse,
    CandidateGene,
    GraphSummaryResponse,
)
from app.core.gene_ranking import rank_candidate_genes, rank_candidate_genes_with_ml
from app.core.graph_loader import load_graph_and_oncogenes
from app.core.rwr import calculate_rwr_p_value, calculate_rwr_proximity, random_rwr_distribution
from app.db.crud import (
    analysis_run_to_response,
    candidate_gene_to_response,
    complete_analysis_run,
    create_analysis_run,
    get_analysis_run,
    get_top_candidate_genes,
    update_analysis_run_status,
)
from app.db.database import DATABASE_ENABLED, SessionLocal
from app.utils.helpers import set_random_seed


router = APIRouter()

BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / "data"
ONCOGENE_FILE = DATA_DIR / "onco_genes.txt"
INTERACTIONS_FILE = DATA_DIR / "interacting_proteins.txt"

ANALYSIS_RESULTS: dict[str, dict] = {}
ANALYSIS_RESULTS_LOCK = Lock()
RUNNING_MESSAGE = "Running graph diffusion and candidate gene ranking..."
COMPLETED_MESSAGE = "Analysis completed successfully."
FAILED_MESSAGE = "Analysis failed. Try a smaller configuration."


def _load_graph_data():
    return load_graph_and_oncogenes(ONCOGENE_FILE, INTERACTIONS_FILE)


def _dump_model(model):
    if hasattr(model, "model_dump"):
        return model.model_dump()
    return model.dict()


def _utc_now_isoformat():
    return datetime.utcnow().isoformat()


def _set_result(run_id: str, result: dict):
    with ANALYSIS_RESULTS_LOCK:
        ANALYSIS_RESULTS[run_id] = result


def _update_result(run_id: str, **updates):
    with ANALYSIS_RESULTS_LOCK:
        if run_id in ANALYSIS_RESULTS:
            ANALYSIS_RESULTS[run_id].update(updates)


def _get_result(run_id: str):
    with ANALYSIS_RESULTS_LOCK:
        result = ANALYSIS_RESULTS.get(run_id)
        return dict(result) if result is not None else None


def _get_db_session():
    if not DATABASE_ENABLED or SessionLocal is None:
        return None
    return SessionLocal()


def _normalize_seed_genes(seed_genes: list[str] | None, default_genes: list[str]) -> list[str]:
    if not seed_genes:
        return default_genes
    return [gene.strip().upper() for gene in seed_genes if gene.strip()]


def _format_rwr_candidates(ranked_genes):
    return [
        _dump_model(
            CandidateGene(
                gene_name=gene_name,
                score=score,
                rank=index,
                rwr_score=score,
            )
        )
        for index, (gene_name, score) in enumerate(ranked_genes, start=1)
    ]


def _mark_analysis_failed(run_id: str, error: Exception):
    error_message = str(error)
    print(f"Analysis failed for {run_id}... {error_message}")

    if DATABASE_ENABLED:
        db = _get_db_session()
        if db is None:
            return
        try:
            update_analysis_run_status(
                db,
                run_id,
                status="failed",
                message=FAILED_MESSAGE,
                error_message=error_message,
            )
        finally:
            db.close()
        return

    _update_result(
        run_id,
        status="failed",
        message=FAILED_MESSAGE,
        error_message=error_message,
        completed_at=_utc_now_isoformat(),
    )


def run_analysis_background(run_id: str, request: AnalyzeRequest, seed_genes: list[str]):
    set_random_seed(42)
    try:
        print(f"Starting analysis for {run_id}...")

        if DATABASE_ENABLED:
            db = _get_db_session()
            if db is None:
                raise RuntimeError("Database session could not be created.")
            try:
                update_analysis_run_status(
                    db,
                    run_id,
                    status="running",
                    message=RUNNING_MESSAGE,
                    error_message=None,
                )

                graph, _ = _load_graph_data()

                rwr_score = calculate_rwr_proximity(
                    graph,
                    seed_genes,
                    restart_probability=request.restart_probability,
                    num_steps=request.num_steps,
                )
                print("Finished RWR proximity")

                random_scores = random_rwr_distribution(
                    graph,
                    seed_genes,
                    num_samples=request.num_random_sets,
                    restart_probability=request.restart_probability,
                    num_steps=request.num_steps,
                )
                print("Finished random distribution")

                p_value = calculate_rwr_p_value(rwr_score, random_scores)

                if request.use_ml_ranking:
                    top_genes = rank_candidate_genes_with_ml(
                        graph,
                        seed_genes,
                        restart_probability=request.restart_probability,
                        num_steps=request.num_steps,
                        top_n=request.top_n,
                    )
                else:
                    ranked_genes = rank_candidate_genes(
                        graph,
                        seed_genes,
                        restart_probability=request.restart_probability,
                        num_steps=request.num_steps,
                        top_n=request.top_n,
                    )
                    top_genes = _format_rwr_candidates(ranked_genes)
                print("Finished candidate ranking")

                complete_analysis_run(db, run_id, rwr_score, p_value, top_genes)
                print(f"Completed analysis for {run_id}...")
            finally:
                db.close()
            return

        _update_result(
            run_id,
            status="running",
            message=RUNNING_MESSAGE,
            error_message=None,
        )

        graph, _ = _load_graph_data()

        rwr_score = calculate_rwr_proximity(
            graph,
            seed_genes,
            restart_probability=request.restart_probability,
            num_steps=request.num_steps,
        )
        print("Finished RWR proximity")

        random_scores = random_rwr_distribution(
            graph,
            seed_genes,
            num_samples=request.num_random_sets,
            restart_probability=request.restart_probability,
            num_steps=request.num_steps,
        )
        print("Finished random distribution")

        p_value = calculate_rwr_p_value(rwr_score, random_scores)

        if request.use_ml_ranking:
            top_genes = rank_candidate_genes_with_ml(
                graph,
                seed_genes,
                restart_probability=request.restart_probability,
                num_steps=request.num_steps,
                top_n=request.top_n,
            )
        else:
            ranked_genes = rank_candidate_genes(
                graph,
                seed_genes,
                restart_probability=request.restart_probability,
                num_steps=request.num_steps,
                top_n=request.top_n,
            )
            top_genes = _format_rwr_candidates(ranked_genes)
        print("Finished candidate ranking")

        _update_result(
            run_id,
            status="completed",
            seed_genes=seed_genes,
            rwr_score=rwr_score,
            p_value=p_value,
            top_genes=top_genes,
            message=COMPLETED_MESSAGE,
            error_message=None,
            completed_at=_utc_now_isoformat(),
        )
        print(f"Completed analysis for {run_id}...")
    except Exception as exc:
        try:
            _mark_analysis_failed(run_id, exc)
        except Exception as update_error:
            print(f"Analysis failed for {run_id}... {exc}")
            print(f"Failed to persist error state for {run_id}... {update_error}")


@router.get("/")
def root():
    return {"message": "OncoGraph API is running"}


@router.get("/health")
def health():
    return {"status": "healthy"}


@router.get("/graph/summary", response_model=GraphSummaryResponse)
def graph_summary():
    graph, oncogenes = _load_graph_data()
    found_genes = [gene for gene in oncogenes if gene in graph.nodes()]
    missing_genes = [gene for gene in oncogenes if gene not in graph.nodes()]

    return GraphSummaryResponse(
        num_nodes=graph.number_of_nodes(),
        num_edges=graph.number_of_edges(),
        num_oncogenes=len(oncogenes),
        oncogenes_found_in_graph=found_genes,
        oncogenes_missing_from_graph=missing_genes,
    )


@router.post("/analyze", response_model=AnalyzeResponse)
def analyze(request: AnalyzeRequest, background_tasks: BackgroundTasks):
    run_id = str(uuid4())

    graph, default_oncogenes = _load_graph_data()
    normalized_seed_genes = _normalize_seed_genes(request.seed_genes, default_oncogenes)
    seed_genes_in_graph = [gene for gene in normalized_seed_genes if gene in graph.nodes()]

    if not seed_genes_in_graph:
        raise HTTPException(status_code=400, detail="None of the provided seed genes were found in the graph.")

    if len(seed_genes_in_graph) < 2:
        raise HTTPException(status_code=400, detail="At least two seed genes found in the graph are required.")

    if DATABASE_ENABLED:
        db = _get_db_session()
        if db is None:
            raise HTTPException(status_code=500, detail="Database session could not be created.")
        try:
            create_analysis_run(db, run_id, request, seed_genes_in_graph)
        finally:
            db.close()
    else:
        running_result = AnalysisResultResponse(
            run_id=run_id,
            status="running",
            seed_genes=seed_genes_in_graph,
            rwr_score=None,
            p_value=None,
            top_genes=[],
            message=RUNNING_MESSAGE,
            error_message=None,
            created_at=_utc_now_isoformat(),
            completed_at=None,
        )
        _set_result(run_id, _dump_model(running_result))

    background_tasks.add_task(run_analysis_background, run_id, request, seed_genes_in_graph)

    return AnalyzeResponse(
        run_id=run_id,
        status="running",
        message="Analysis started",
    )


@router.get("/status/{run_id}", response_model=AnalyzeResponse)
def get_status(run_id: str):
    if DATABASE_ENABLED:
        db = _get_db_session()
        if db is None:
            raise HTTPException(status_code=500, detail="Database session could not be created.")
        try:
            run = get_analysis_run(db, run_id)
            if run is None:
                raise HTTPException(status_code=404, detail="Run ID not found.")
            return AnalyzeResponse(
                run_id=run.id,
                status=run.status,
                message=run.message,
            )
        finally:
            db.close()

    result = _get_result(run_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Run ID not found.")
    return AnalyzeResponse(
        run_id=run_id,
        status=result["status"],
        message=result["message"],
    )


@router.get("/results/{run_id}", response_model=AnalysisResultResponse)
def get_results(run_id: str):
    if DATABASE_ENABLED:
        db = _get_db_session()
        if db is None:
            raise HTTPException(status_code=500, detail="Database session could not be created.")
        try:
            run = get_analysis_run(db, run_id)
            if run is None:
                raise HTTPException(status_code=404, detail="Run ID not found.")
            return analysis_run_to_response(run)
        finally:
            db.close()

    result = _get_result(run_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Run ID not found.")
    return AnalysisResultResponse(**result)


@router.get("/genes/top-candidates/{run_id}", response_model=list[CandidateGene])
def get_top_candidates(run_id: str):
    if DATABASE_ENABLED:
        db = _get_db_session()
        if db is None:
            raise HTTPException(status_code=500, detail="Database session could not be created.")
        try:
            run = get_analysis_run(db, run_id)
            if run is None:
                raise HTTPException(status_code=404, detail="Run ID not found.")
            if run.status != "completed":
                raise HTTPException(status_code=400, detail="Analysis is not completed yet.")
            candidates = get_top_candidate_genes(db, run_id)
            return [candidate_gene_to_response(candidate) for candidate in candidates]
        finally:
            db.close()

    result = _get_result(run_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Run ID not found.")
    if result["status"] != "completed":
        raise HTTPException(status_code=400, detail="Analysis is not completed yet.")
    return [CandidateGene(**gene) for gene in result.get("top_genes", [])]
