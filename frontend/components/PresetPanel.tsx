import type { AnalysisFormValues } from "./AnalysisForm";


type PresetPanelProps = {
  onSelectPreset: (values: AnalysisFormValues) => void;
};

type Preset = {
  name: string;
  description: string;
  runtimeLabel: string;
  tags: string[];
  values: AnalysisFormValues;
};


const PRESETS: Preset[] = [
  {
    name: "Quick RWR Test",
    description: "Fast RWR-only run for checking that the pipeline works.",
    runtimeLabel: "Usually a few seconds",
    tags: ["Fast", "RWR"],
    values: {
      seedGenesInput: "",
      restartProbability: "0.3",
      numSteps: "1000",
      numRandomSets: "20",
      topN: "10",
      useMlRanking: false,
    },
  },
  {
    name: "ML Ranking Demo",
    description: "Runs graph feature engineering and ML-enhanced ranking.",
    runtimeLabel: "Usually 30–90 seconds",
    tags: ["ML", "Demo"],
    values: {
      seedGenesInput: "",
      restartProbability: "0.3",
      numSteps: "1000",
      numRandomSets: "20",
      topN: "10",
      useMlRanking: true,
    },
  },
  {
    name: "Higher Confidence RWR",
    description: "More simulations for a slower but more stable RWR result.",
    runtimeLabel: "Usually 15–60 seconds",
    tags: ["Stable", "RWR"],
    values: {
      seedGenesInput: "",
      restartProbability: "0.3",
      numSteps: "10000",
      numRandomSets: "200",
      topN: "15",
      useMlRanking: false,
    },
  },
  {
    name: "TP53 Focus Test",
    description: "Uses TP53-centered seeds to test custom seed behavior.",
    runtimeLabel: "Usually a few seconds",
    tags: ["Custom Seeds", "Focused"],
    values: {
      seedGenesInput: "TP53, SMARCA4, ARID1A",
      restartProbability: "0.3",
      numSteps: "1000",
      numRandomSets: "20",
      topN: "10",
      useMlRanking: false,
    },
  },
  {
    name: "Local Exploration",
    description: "Higher restart probability to keep diffusion closer to seed genes.",
    runtimeLabel: "Usually a few seconds",
    tags: ["Local", "RWR"],
    values: {
      seedGenesInput: "",
      restartProbability: "0.6",
      numSteps: "1000",
      numRandomSets: "20",
      topN: "10",
      useMlRanking: false,
    },
  },
  {
    name: "Global Exploration",
    description: "Lower restart probability to let the random walk explore farther across the network.",
    runtimeLabel: "Usually a few seconds",
    tags: ["Global", "RWR"],
    values: {
      seedGenesInput: "",
      restartProbability: "0.1",
      numSteps: "1000",
      numRandomSets: "20",
      topN: "10",
      useMlRanking: false,
    },
  },
];


export function PresetPanel({ onSelectPreset }: PresetPanelProps) {
  return (
    <aside className="rounded-2xl border border-line bg-panel p-6 shadow-sm md:p-8">
      <div className="mb-6">
        <h2 className="text-2xl font-semibold text-ink">Try a Preset</h2>
        <p className="mt-2 text-sm text-muted">
          Click a preset to fill the form with demo-friendly analysis settings.
        </p>
      </div>

      <div className="grid gap-4">
        {PRESETS.map((preset) => (
          <button
            key={preset.name}
            type="button"
            onClick={() => onSelectPreset(preset.values)}
            className="rounded-xl border border-line bg-white p-4 text-left transition hover:-translate-y-0.5 hover:border-accent hover:shadow-sm focus:outline-none focus:ring-2 focus:ring-accent/10"
          >
            <div className="flex flex-wrap gap-2">
              {preset.tags.map((tag) => (
                <span
                  key={`${preset.name}-${tag}`}
                  className="rounded-lg bg-slate-100 px-2.5 py-1 text-xs font-medium text-muted"
                >
                  {tag}
                </span>
              ))}
            </div>
            <h3 className="mt-3 text-base font-semibold text-ink">{preset.name}</h3>
            <p className="mt-2 text-sm leading-6 text-muted">{preset.description}</p>
            <p className="mt-3 text-xs font-medium uppercase tracking-wide text-accent">
              {preset.runtimeLabel}
            </p>
          </button>
        ))}
      </div>
    </aside>
  );
}
