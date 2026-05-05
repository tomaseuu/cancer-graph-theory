"use client";

import { useState } from "react";

import {
  AnalysisForm,
  type AnalysisFormValues,
} from "../components/AnalysisForm";
import { PresetPanel } from "../components/PresetPanel";


const DEFAULT_FORM_VALUES: AnalysisFormValues = {
  seedGenesInput: "",
  restartProbability: "0.3",
  numSteps: "1000",
  numRandomSets: "20",
  topN: "10",
  useMlRanking: false,
};

export default function HomePage() {
  const [formValues, setFormValues] = useState<AnalysisFormValues>(DEFAULT_FORM_VALUES);

  function handleFormChange(updates: Partial<AnalysisFormValues>) {
    setFormValues((currentValues) => ({
      ...currentValues,
      ...updates,
    }));
  }

  function handleSelectPreset(values: AnalysisFormValues) {
    setFormValues(values);
  }

  return (
    <main className="min-h-screen px-4 py-10 md:px-8">
      <div className="mx-auto flex max-w-6xl flex-col gap-6">
        <section className="rounded-2xl border border-line bg-panel p-8 shadow-sm">
          <p className="text-sm font-medium uppercase tracking-[0.2em] text-accent">
            Bioinformatics - Cancer and Graph Theory
          </p>
          <h1 className="mt-3 text-4xl font-semibold tracking-tight text-ink md:text-5xl">
            OncoGraph
          </h1>
          <p className="mt-4 max-w-3xl text-lg text-muted">
            Graph-based cancer gene discovery using network diffusion and ML
            ranking.
          </p>
          <p className="mt-4 max-w-3xl text-sm text-muted">
            For educational and research exploration only. Not for medical
            decision-making.
          </p>
        </section>

        <section className="grid gap-6 xl:grid-cols-3">
          <div className="xl:col-span-2">
            <AnalysisForm values={formValues} onChange={handleFormChange} />
          </div>
          <div className="xl:col-span-1">
            <PresetPanel onSelectPreset={handleSelectPreset} />
          </div>
        </section>
      </div>
    </main>
  );
}
