import type { CandidateGene } from "../lib/api";


type ResultsTableProps = {
  genes: CandidateGene[];
};


function formatNumber(value: number | null | undefined, digits = 4) {
  if (value === null || value === undefined) {
    return "-";
  }
  return value.toFixed(digits);
}


export function ResultsTable({ genes }: ResultsTableProps) {
  if (genes.length === 0) {
    return null;
  }

  return (
    <div className="overflow-hidden rounded-2xl border border-line bg-panel shadow-sm">
      <div className="overflow-x-auto">
        <table className="min-w-full text-left text-sm">
          <thead className="bg-slate-50 text-muted">
            <tr>
              <th className="px-4 py-3 font-medium">Rank</th>
              <th className="px-4 py-3 font-medium">Gene</th>
              <th className="px-4 py-3 font-medium">Score</th>
              <th className="px-4 py-3 font-medium">RWR Score</th>
              <th className="px-4 py-3 font-medium">ML Probability</th>
              <th className="px-4 py-3 font-medium">Final Score</th>
              <th className="px-4 py-3 font-medium">Degree</th>
            </tr>
          </thead>
          <tbody>
            {genes.map((gene) => (
              <tr key={`${gene.gene_name}-${gene.rank}`} className="border-t border-line">
                <td className="px-4 py-3 text-ink">{gene.rank}</td>
                <td className="px-4 py-3 font-medium text-ink">{gene.gene_name}</td>
                <td className="px-4 py-3 text-muted">{formatNumber(gene.score)}</td>
                <td className="px-4 py-3 text-muted">{formatNumber(gene.rwr_score)}</td>
                <td className="px-4 py-3 text-muted">{formatNumber(gene.ml_probability)}</td>
                <td className="px-4 py-3 text-muted">{formatNumber(gene.final_score)}</td>
                <td className="px-4 py-3 text-muted">
                  {gene.degree === null || gene.degree === undefined ? "-" : gene.degree}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
