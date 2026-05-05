"use client";

import type { AnalysisResult } from "../lib/api";
import { buildResultSummary } from "../lib/summary";


type ResultSummaryProps = {
  result: AnalysisResult;
  useMlRanking: boolean;
};


function getBadgeStyles(label: string) {
  if (label === "Strong network signal") {
    return "border-emerald-200 bg-emerald-50 text-emerald-800";
  }

  if (label === "Statistically significant network signal") {
    return "border-sky-200 bg-sky-50 text-sky-800";
  }

  return "border-slate-200 bg-slate-100 text-slate-700";
}


export function ResultSummary({ result, useMlRanking }: ResultSummaryProps) {
  const summary = buildResultSummary(result, useMlRanking);

  return (
    <section className="rounded-2xl border border-line bg-panel p-6 shadow-sm">
      <div className="flex flex-col gap-3 md:flex-row md:items-start md:justify-between">
        <div>
          <h3 className="text-xl font-semibold text-ink">What this result means</h3>
          <p className="mt-2 text-sm text-muted">
            A plain-English summary of the network signal and the top-ranked candidate.
          </p>
        </div>
        <span
          className={`inline-flex rounded-lg border px-3 py-1 text-sm font-medium ${getBadgeStyles(summary.significanceLabel)}`}
        >
          {summary.significanceLabel}
        </span>
      </div>

      <div className="mt-5 space-y-5 text-sm leading-7 text-muted">
        <div>
          <p className="font-medium text-ink">Main result</p>
          <p className="mt-1">
            This run compared the selected seed genes against random gene sets in the protein
            interaction network. <strong className="font-semibold text-ink">RWR score:</strong>{" "}
            {summary.formattedRwrScore}. <strong className="font-semibold text-ink">p-value:</strong>{" "}
            {summary.formattedPValue}.{" "}
            <strong className="font-semibold text-ink">{summary.significance}</strong>
          </p>
        </div>

        <div>
          <p className="font-medium text-ink">Top candidate</p>
          <p className="mt-1">
            {summary.topCandidate ? (
              <>
                <strong className="font-semibold text-ink">Top candidate:</strong>{" "}
                {summary.topCandidate}. {summary.topCandidateSummary}
              </>
            ) : (
              summary.topCandidateSummary
            )}
          </p>
        </div>

        <div>
          <p className="font-medium text-ink">How to interpret it</p>
          <p className="mt-1">{summary.interpretation}</p>
        </div>
      </div>

      <p className="mt-6 text-xs leading-6 text-muted">
        These results are exploratory and are not medical or clinical recommendations.
      </p>
    </section>
  );
}
