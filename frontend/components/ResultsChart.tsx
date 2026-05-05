"use client";

import {
  Bar,
  BarChart,
  CartesianGrid,
  Legend,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

import type { CandidateGene } from "../lib/api";


type ResultsChartProps = {
  topGenes: CandidateGene[];
};

type TooltipValue = number | string;

type ChartDatum = {
  gene_name: string;
  score: number;
  rwr_score?: number | null;
  ml_probability?: number | null;
  final_score?: number | null;
};


function formatScore(value: number | null | undefined) {
  if (value === null || value === undefined || Number.isNaN(value)) {
    return "-";
  }
  return value.toFixed(4);
}


export function ResultsChart({ topGenes }: ResultsChartProps) {
  if (topGenes.length === 0) {
    return null;
  }

  const chartData: ChartDatum[] = topGenes.slice(0, 10).map((gene) => ({
    gene_name: gene.gene_name,
    score: gene.score,
    rwr_score: gene.rwr_score ?? null,
    ml_probability: gene.ml_probability ?? null,
    final_score: gene.final_score ?? null,
  }));

  const showMlBreakdown = chartData.some(
    (gene) => gene.ml_probability !== null || gene.final_score !== null,
  );

  return (
    <section className="flex flex-col gap-4">
      <div className="rounded-2xl border border-line bg-panel p-6 shadow-sm">
        <div className="mb-5">
          <h3 className="text-xl font-semibold text-ink">Top Candidate Genes by Score</h3>
          <p className="mt-2 text-sm text-muted">
            Higher scores indicate genes that are more strongly prioritized by the graph
            diffusion and ranking pipeline.
          </p>
        </div>

        <div className="h-80 w-full">
          <ResponsiveContainer width="100%" height="100%">
            <BarChart data={chartData} margin={{ top: 8, right: 12, left: 0, bottom: 12 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#dbe3dd" />
              <XAxis
                dataKey="gene_name"
                tick={{ fill: "#5f6d65", fontSize: 12 }}
                axisLine={{ stroke: "#dbe3dd" }}
                tickLine={false}
              />
              <YAxis
                tick={{ fill: "#5f6d65", fontSize: 12 }}
                axisLine={{ stroke: "#dbe3dd" }}
                tickLine={false}
                tickFormatter={(value: number) => value.toFixed(4)}
              />
              <Tooltip formatter={(value: TooltipValue) => formatScore(Number(value))} />
              <Bar dataKey="score" fill="#1f4f46" radius={[8, 8, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>
      </div>

      {showMlBreakdown && (
        <div className="rounded-2xl border border-line bg-panel p-6 shadow-sm">
          <div className="mb-5">
            <h3 className="text-xl font-semibold text-ink">RWR vs ML Ranking Scores</h3>
            <p className="mt-2 text-sm text-muted">
              ML mode combines Random Walk with Restart scores with graph-based features such
              as degree, PageRank, centrality, and seed-gene distance.
            </p>
          </div>

          <div className="h-80 w-full">
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={chartData} margin={{ top: 8, right: 12, left: 0, bottom: 12 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="#dbe3dd" />
                <XAxis
                  dataKey="gene_name"
                  tick={{ fill: "#5f6d65", fontSize: 12 }}
                  axisLine={{ stroke: "#dbe3dd" }}
                  tickLine={false}
                />
                <YAxis
                  tick={{ fill: "#5f6d65", fontSize: 12 }}
                  axisLine={{ stroke: "#dbe3dd" }}
                  tickLine={false}
                  tickFormatter={(value: number) => value.toFixed(4)}
                />
                <Tooltip formatter={(value: TooltipValue) => formatScore(Number(value))} />
                <Legend />
                <Bar dataKey="rwr_score" name="RWR score" fill="#3f7d73" radius={[6, 6, 0, 0]} />
                <Bar
                  dataKey="ml_probability"
                  name="ML probability"
                  fill="#c98d3f"
                  radius={[6, 6, 0, 0]}
                />
                <Bar
                  dataKey="final_score"
                  name="Final score"
                  fill="#16211d"
                  radius={[6, 6, 0, 0]}
                />
              </BarChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}
    </section>
  );
}
