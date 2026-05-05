import type { AnalysisResult } from "./api";


export function formatScore(value: number | null | undefined, digits = 4) {
  if (value === null || value === undefined || Number.isNaN(value)) {
    return "not available";
  }

  return value.toFixed(digits);
}


export function getSignificanceLabel(pValue: number | null | undefined) {
  if (pValue === null || pValue === undefined || Number.isNaN(pValue)) {
    return "Exploratory / not statistically significant";
  }

  if (pValue < 0.01) {
    return "Strong network signal";
  }

  if (pValue < 0.05) {
    return "Statistically significant network signal";
  }

  return "Exploratory / not statistically significant";
}


export function buildResultSummary(result: AnalysisResult, useMlRanking: boolean) {
  const rwrScore = formatScore(result.rwr_score, 6);
  const pValue = formatScore(result.p_value, 6);
  const topCandidate = result.top_genes[0]?.gene_name ?? null;

  const mainResult = `This run compared the selected seed genes against random gene sets in the protein interaction network. The RWR proximity score was ${rwrScore}, with a p-value of ${pValue}.`;

  const significance =
    result.p_value !== null && result.p_value !== undefined && !Number.isNaN(result.p_value)
      ? result.p_value < 0.05
        ? "The p-value suggests the seed genes are more connected than expected by random chance."
        : "The p-value does not show strong evidence that the seed genes are more connected than random sets."
      : "The p-value was not available for this run, so the network signal should be treated as exploratory.";

  const topCandidateSummary = topCandidate
    ? `The highest-ranked candidate gene was ${topCandidate}, meaning it received the strongest score from the ranking pipeline in this run.`
    : "No candidate genes were returned.";

  const interpretation = useMlRanking
    ? "ML mode combines Random Walk with Restart scores with graph features such as PageRank, centrality, and distance to seed genes."
    : "Scores are based mainly on Random Walk with Restart, which spreads influence from seed genes through the protein interaction network.";

  return {
    formattedRwrScore: rwrScore,
    formattedPValue: pValue,
    topCandidate,
    mainResult,
    significance,
    topCandidateSummary,
    interpretation,
    significanceLabel: getSignificanceLabel(result.p_value),
  };
}
