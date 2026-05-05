type StatusBadgeProps = {
  status: string;
};


const STATUS_STYLES: Record<string, string> = {
  running: "bg-amber-100 text-amber-800 border-amber-200",
  completed: "bg-emerald-100 text-emerald-800 border-emerald-200",
  failed: "bg-rose-100 text-rose-800 border-rose-200",
  timeout: "bg-orange-100 text-orange-800 border-orange-200",
};


export function StatusBadge({ status }: StatusBadgeProps) {
  const normalizedStatus = status.toLowerCase();
  const style =
    STATUS_STYLES[normalizedStatus] ?? "bg-slate-100 text-slate-700 border-slate-200";

  return (
    <span className={`inline-flex rounded-lg border px-3 py-1 text-sm font-medium ${style}`}>
      {normalizedStatus}
    </span>
  );
}
