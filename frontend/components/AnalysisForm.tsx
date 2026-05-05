"use client";

import { FormEvent, useEffect, useMemo, useRef, useState } from "react";

import { getResults, getStatus, startAnalysis, type AnalysisResult } from "../lib/api";
import { ResultSummary } from "./ResultSummary";
import { ResultsChart } from "./ResultsChart";
import { ResultsTable } from "./ResultsTable";
import { StatusBadge } from "./StatusBadge";


const READY_MESSAGE = "Ready to run analysis.";
const COMPLETED_MESSAGE = "Analysis completed successfully.";
const FAILED_MESSAGE = "Analysis failed. Please try a smaller configuration.";
const TIMEOUT_HELP_MESSAGE =
  "This is taking longer than the estimate, but the backend may still be processing. Try a smaller preset if it does not finish.";
const API_ERROR_MESSAGE =
  "Could not connect to the OncoGraph API. Make sure the FastAPI server is running on http://localhost:8000.";
const POLL_INTERVAL_MS = 1000;
const POLL_TIMEOUT_MS = 120000;


export type AnalysisFormValues = {
  seedGenesInput: string;
  restartProbability: string;
  numSteps: string;
  numRandomSets: string;
  topN: string;
  useMlRanking: boolean;
};

type AnalysisFormProps = {
  values: AnalysisFormValues;
  onChange: (updates: Partial<AnalysisFormValues>) => void;
};


function formatMetric(value: number | null | undefined, digits = 4) {
  if (value === null || value === undefined) {
    return "-";
  }
  return value.toFixed(digits);
}


function parseSeedGenes(value: string) {
  return value
    .split(/[\n,]+/)
    .map((gene) => gene.trim())
    .filter(Boolean);
}


function getEstimatedRuntime(values: AnalysisFormValues) {
  const numSteps = Number(values.numSteps);
  const numRandomSets = Number(values.numRandomSets);

  if (!values.useMlRanking) {
    if (numSteps <= 1000 && numRandomSets <= 20) {
      return "Estimated time: a few seconds.";
    }
    if (numSteps <= 10000 && numRandomSets <= 200) {
      return "Estimated time: 15–60 seconds.";
    }
    return "Estimated time: 1–3+ minutes.";
  }

  if (numSteps <= 1000 && numRandomSets <= 20) {
    return "Estimated time: 30–90 seconds because ML ranking computes extra graph features.";
  }

  return "Estimated time: 2–5+ minutes because ML ranking computes extra graph features.";
}


function getRunningMessage(values: AnalysisFormValues) {
  if (values.useMlRanking) {
    return "Running ML-enhanced graph ranking. This computes extra graph features and may take 30–90 seconds.";
  }

  return "Running graph diffusion and candidate gene ranking. This usually takes a few seconds for quick settings.";
}


export function AnalysisForm({ values, onChange }: AnalysisFormProps) {
  const [runId, setRunId] = useState<string | null>(null);
  const [status, setStatus] = useState<string | null>(null);
  const [message, setMessage] = useState<string>(READY_MESSAGE);
  const [error, setError] = useState<string | null>(null);
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [result, setResult] = useState<AnalysisResult | null>(null);
  const [elapsedSeconds, setElapsedSeconds] = useState(0);
  const [lastSubmittedUseMlRanking, setLastSubmittedUseMlRanking] = useState(false);
  const intervalRef = useRef<number | null>(null);
  const pollingStartedAtRef = useRef<number | null>(null);
  const isPollingRef = useRef(false);

  const seedGenesPreview = useMemo(
    () => parseSeedGenes(values.seedGenesInput),
    [values.seedGenesInput],
  );
  const estimatedRuntime = useMemo(() => getEstimatedRuntime(values), [values]);
  const hasExceededEstimate = status === "running" && elapsedSeconds >= POLL_TIMEOUT_MS / 1000;

  function stopPolling() {
    if (intervalRef.current !== null) {
      window.clearInterval(intervalRef.current);
      intervalRef.current = null;
    }
    isPollingRef.current = false;
    pollingStartedAtRef.current = null;
  }

  function syncElapsedTime() {
    if (pollingStartedAtRef.current !== null) {
      setElapsedSeconds(Math.floor((Date.now() - pollingStartedAtRef.current) / 1000));
    }
  }

  useEffect(() => {
    return () => {
      stopPolling();
    };
  }, []);

  async function fetchFinalResult(currentRunId: string) {
    const fullResult = await getResults(currentRunId);
    setResult(fullResult);

    if (fullResult.status === "completed") {
      setStatus("completed");
      setMessage(COMPLETED_MESSAGE);
      setError(null);
      return;
    }

    if (fullResult.status === "failed") {
      setStatus("failed");
      setMessage(FAILED_MESSAGE);
      setError(fullResult.error_message || FAILED_MESSAGE);
    }
  }

  async function pollRun(currentRunId: string) {
    if (isPollingRef.current) {
      return;
    }

    isPollingRef.current = true;

    try {
      const currentStatus = await getStatus(currentRunId);
      setStatus(currentStatus.status);

      if (currentStatus.status === "completed") {
        syncElapsedTime();
        stopPolling();
        setMessage(COMPLETED_MESSAGE);
        await fetchFinalResult(currentRunId);
        return;
      }

      if (currentStatus.status === "failed") {
        syncElapsedTime();
        stopPolling();
        setStatus("failed");
        setMessage(FAILED_MESSAGE);
        await fetchFinalResult(currentRunId);
        return;
      }

    } catch (pollError) {
      stopPolling();
      setStatus("failed");
      setMessage(FAILED_MESSAGE);
      setError(pollError instanceof Error ? pollError.message : API_ERROR_MESSAGE);
    } finally {
      isPollingRef.current = false;
    }
  }

  function startPolling(currentRunId: string) {
    stopPolling();
    pollingStartedAtRef.current = Date.now();
    setElapsedSeconds(0);
    void pollRun(currentRunId);
    intervalRef.current = window.setInterval(() => {
      if (pollingStartedAtRef.current !== null) {
        setElapsedSeconds(Math.floor((Date.now() - pollingStartedAtRef.current) / 1000));
      }
      void pollRun(currentRunId);
    }, POLL_INTERVAL_MS);
  }

  async function handleSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    stopPolling();
    setIsSubmitting(true);
    setError(null);
    setResult(null);
    setRunId(null);
    setStatus(null);
    setMessage(READY_MESSAGE);
    setElapsedSeconds(0);

    try {
      const seedGenes = parseSeedGenes(values.seedGenesInput);
      setLastSubmittedUseMlRanking(values.useMlRanking);
      const response = await startAnalysis({
        seed_genes: seedGenes.length > 0 ? seedGenes : undefined,
        restart_probability: Number(values.restartProbability),
        num_steps: Number(values.numSteps),
        num_random_sets: Number(values.numRandomSets),
        top_n: Number(values.topN),
        use_ml_ranking: values.useMlRanking,
      });

      setRunId(response.run_id);
      setStatus(response.status);
      setMessage(getRunningMessage(values));
      startPolling(response.run_id);
    } catch (submitError) {
      setError(submitError instanceof Error ? submitError.message : API_ERROR_MESSAGE);
      setStatus("failed");
      setMessage(FAILED_MESSAGE);
    } finally {
      setIsSubmitting(false);
    }
  }

  const isRunning = status === "running" || isSubmitting;
  const elapsedLabel =
    status === "completed"
      ? `Completed in ${elapsedSeconds}s`
      : status === "failed"
        ? `Failed after ${elapsedSeconds}s`
        : `Elapsed time: ${elapsedSeconds}s`;

  return (
    <div className="flex flex-col gap-6">
      <section className="rounded-2xl border border-line bg-panel p-6 shadow-sm md:p-8">
        <div className="mb-6">
          <h2 className="text-2xl font-semibold text-ink">Run Analysis</h2>
          <p className="mt-2 text-sm text-muted">
            Leave seed genes blank to use the backend&apos;s default oncogene seed set.
          </p>
        </div>

        <form className="grid gap-4" onSubmit={handleSubmit}>
          <label className="grid gap-2">
            <span className="text-xs font-semibold uppercase tracking-wide text-muted">Seed genes</span>
            <textarea
              className="min-h-[96px] rounded-xl border border-line bg-white px-3 py-2.5 text-sm text-ink outline-none transition focus:border-accent focus:ring-2 focus:ring-accent/10"
              placeholder="ARID1A, KIF1B, DST, PKHD1, TP53, COL5A3, SMARCA4"
              value={values.seedGenesInput}
              onChange={(event) => onChange({ seedGenesInput: event.target.value })}
            />
          </label>

          <div className="grid grid-cols-1 gap-x-5 gap-y-4 sm:grid-cols-2 lg:[grid-template-columns:repeat(4,minmax(0,220px))] lg:justify-between">
            <label className="grid min-w-0 max-w-[220px] gap-2">
              <span className="text-[11px] font-semibold uppercase tracking-[0.14em] text-muted">Restart probability</span>
              <input
                className="h-10 w-full rounded-xl border border-line bg-white px-3 py-2 text-sm text-ink outline-none transition focus:border-accent focus:ring-2 focus:ring-accent/10"
                type="number"
                step="0.1"
                value={values.restartProbability}
                onChange={(event) => onChange({ restartProbability: event.target.value })}
              />
            </label>

            <label className="grid min-w-0 max-w-[220px] gap-2">
              <span className="text-[11px] font-semibold uppercase tracking-[0.14em] text-muted">Num steps</span>
              <input
                className="h-10 w-full rounded-xl border border-line bg-white px-3 py-2 text-sm text-ink outline-none transition focus:border-accent focus:ring-2 focus:ring-accent/10"
                type="number"
                value={values.numSteps}
                onChange={(event) => onChange({ numSteps: event.target.value })}
              />
            </label>

            <label className="grid min-w-0 max-w-[220px] gap-2">
              <span className="text-[11px] font-semibold uppercase tracking-[0.14em] text-muted">Num random sets</span>
              <input
                className="h-10 w-full rounded-xl border border-line bg-white px-3 py-2 text-sm text-ink outline-none transition focus:border-accent focus:ring-2 focus:ring-accent/10"
                type="number"
                value={values.numRandomSets}
                onChange={(event) => onChange({ numRandomSets: event.target.value })}
              />
            </label>

            <label className="grid min-w-0 max-w-[220px] gap-2">
              <span className="text-[11px] font-semibold uppercase tracking-[0.14em] text-muted">Top N</span>
              <input
                className="h-10 w-full rounded-xl border border-line bg-white px-3 py-2 text-sm text-ink outline-none transition focus:border-accent focus:ring-2 focus:ring-accent/10"
                type="number"
                value={values.topN}
                onChange={(event) => onChange({ topN: event.target.value })}
              />
            </label>
          </div>

          <label className="flex items-center gap-3 rounded-xl border border-line bg-slate-50 px-3 py-3 text-sm text-ink">
            <input
              className="h-4 w-4 rounded border-line text-accent focus:ring-accent"
              type="checkbox"
              checked={values.useMlRanking}
              onChange={(event) => onChange({ useMlRanking: event.target.checked })}
            />
            Use ML ranking
          </label>

          <div className="flex flex-col gap-3 sm:flex-row sm:items-end sm:justify-between">
            <div className="space-y-1">
              <p className="text-sm text-muted">
                {seedGenesPreview.length > 0
                  ? `${seedGenesPreview.length} custom seed genes provided`
                  : "Using default seed genes from the backend"}
              </p>
              <p className="text-sm text-muted">{estimatedRuntime}</p>
            </div>
            <div className="flex items-center justify-start sm:justify-end">
              <button
                className="inline-flex h-11 min-w-36 items-center justify-center rounded-xl bg-ink px-5 py-2.5 text-sm font-medium text-white transition hover:bg-slate-800 focus:outline-none focus:ring-2 focus:ring-accent/20 disabled:cursor-not-allowed disabled:opacity-60"
                type="submit"
                disabled={isRunning}
              >
                {isRunning ? "Running..." : "Run Analysis"}
              </button>
            </div>
          </div>
        </form>
      </section>

      <section className="rounded-2xl border border-line bg-panel p-6 shadow-sm md:p-8">
        <div className="flex flex-col gap-3 md:flex-row md:items-center md:justify-between">
          <div>
            <h3 className="text-xl font-semibold text-ink">Run Status</h3>
            <p className="mt-2 text-sm text-muted">{message}</p>
          </div>
          {status && <StatusBadge status={status} />}
        </div>

        {status === "running" && (
          <div className="mt-4 space-y-2 text-sm text-muted">
            <div className="flex items-center gap-3">
              <span className="h-2.5 w-2.5 rounded-full bg-amber-500" />
              <span>Analysis is still running...</span>
            </div>
            <p>{elapsedLabel}</p>
          </div>
        )}

        {(status === "completed" || status === "failed") && (
          <div className="mt-4 text-sm text-muted">{elapsedLabel}</div>
        )}

        {hasExceededEstimate && (
          <div className="mt-4 rounded-xl border border-orange-200 bg-orange-50 px-4 py-3 text-sm text-orange-700">
            {TIMEOUT_HELP_MESSAGE}
          </div>
        )}

        {error && (
          <div
            className={`mt-4 rounded-xl px-4 py-3 text-sm ${
              hasExceededEstimate
                ? "border border-orange-200 bg-orange-50 text-orange-700"
                : "border border-rose-200 bg-rose-50 text-rose-700"
            }`}
          >
            {error}
          </div>
        )}
      </section>

      {result?.status === "completed" && (
        <>
          <section className="grid gap-4 md:grid-cols-2 xl:grid-cols-4">
            <div className="rounded-2xl border border-line bg-panel p-5 shadow-sm">
              <p className="text-sm text-muted">RWR score</p>
              <p className="mt-2 text-2xl font-semibold text-ink">
                {formatMetric(result.rwr_score, 6)}
              </p>
            </div>
            <div className="rounded-2xl border border-line bg-panel p-5 shadow-sm">
              <p className="text-sm text-muted">p-value</p>
              <p className="mt-2 text-2xl font-semibold text-ink">
                {formatMetric(result.p_value, 6)}
              </p>
            </div>
            <div className="rounded-2xl border border-line bg-panel p-5 shadow-sm">
              <p className="text-sm text-muted">Seed genes</p>
              <p className="mt-2 text-2xl font-semibold text-ink">{result.seed_genes.length}</p>
            </div>
            <div className="rounded-2xl border border-line bg-panel p-5 shadow-sm">
              <p className="text-sm text-muted">Top candidates</p>
              <p className="mt-2 text-2xl font-semibold text-ink">{result.top_genes.length}</p>
            </div>
          </section>

          <ResultSummary result={result} useMlRanking={lastSubmittedUseMlRanking} />

          <ResultsChart topGenes={result.top_genes} />

          <section className="flex flex-col gap-4">
            <div>
              <h3 className="text-xl font-semibold text-ink">Candidate Genes</h3>
              <p className="mt-2 text-sm text-muted">
                Ranked genes returned by the OncoGraph API.
              </p>
            </div>
            <ResultsTable genes={result.top_genes} />
          </section>
        </>
      )}
    </div>
  );
}
