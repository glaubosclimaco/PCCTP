"""Descriptive geometry sensitivity for completed studies, never an active gate.

Reconstruct effects from audited seed costs using exact rational arithmetic.
No optimization, active 300-second outcomes, or held-out test inputs are read.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STUDIES = (
    ("TRAIN", "ruin_total60_absolute", 9, (42, 43, 44)),
    ("Reused validation", "ruin_validation60", 7, tuple(range(45, 50))),
)


def summarize(audit, summary, count, seeds):
    assert audit["status"] == "PASS" and audit["batch_complete"]
    assert summary["status"] == "COMPLETE"
    groups = {}
    for row in audit["rows"]:
        assert row["feasible"] and row["mode"] in ("mdm", "hybrid")
        assert row["seed"] in seeds
        records = groups.setdefault(row["instance"], {})
        key = (row["seed"], row["mode"])
        assert key not in records
        records[key] = Fraction(str(row["cost"]))
    assert len(groups) == count
    recorded = {row["instance"]: row for row in summary["rows"]}
    assert set(recorded) == set(groups)
    effects = {}
    rows = []
    for name, records in sorted(groups.items()):
        assert set(records) == {(s, m) for s in seeds for m in ("mdm", "hybrid")}
        baseline = sum(records[(s, "mdm")] for s in seeds)
        hybrid = sum(records[(s, "hybrid")] for s in seeds)
        assert baseline > 0
        gain = 100 * (baseline - hybrid) / baseline
        assert abs(float(gain) - recorded[name]["gain_percent"]) < 1e-10
        effects[name] = gain
        rows.append({"instance": name, "gain_percent": float(gain),
                     "gain_exact": str(gain)})
    gains = list(effects.values())
    mean = sum(gains) / count
    assert abs(float(mean) - summary["mean_gain_percent"]) < 1e-10
    omissions = [{"omitted_instance": name,
                  "remaining_mean_percent": float((sum(gains) - gain) / (count - 1))}
                 for name, gain in effects.items()]
    # Retrospective diagnostic only: never drop these cases in the primary mean.
    top_two = sorted(effects, key=lambda name: (-effects[name], name))[:2]
    remainder = [gain for name, gain in effects.items() if name not in top_two]
    return {"geometries": count, "seeds_per_geometry": len(seeds),
            "mean_percent": float(mean), "median_percent": float(median(gains)),
            "min_percent": float(min(gains)), "max_percent": float(max(gains)),
            "wins": sum(g > 0 for g in gains), "losses": sum(g < 0 for g in gains),
            "ties": sum(g == 0 for g in gains), "geometry_effects": rows,
            "leave_one_geometry_out": omissions,
            "loo_min_percent": min(x["remaining_mean_percent"] for x in omissions),
            "loo_max_percent": max(x["remaining_mean_percent"] for x in omissions),
            "retrospective_two_largest": top_two,
            "remaining_after_two_largest_mean_percent": float(sum(remainder) / len(remainder))}


if __name__ == "__main__":
    studies = {}
    sources = [Path(__file__)]
    for label, prefix, count, seeds in STUDIES:
        audit_path = HERE / (prefix + "_audit.json")
        summary_path = HERE / (prefix + "_summary.json")
        sources.extend((audit_path, summary_path))
        studies[label] = summarize(json.loads(audit_path.read_text()),
                                  json.loads(summary_path.read_text()), count, seeds)
    result = {"status": "PASS", "scope": "Post hoc descriptive sensitivity only; "
              "all primary denominators and failed gates retained; no p-values, "
              "confidence intervals, retuning, or independent generalization claim.",
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources}, "studies": studies}
    (HERE / "completed60_sensitivity.json").write_text(json.dumps(result, indent=2) + "\n")
    table = [r"\begin{table}[htbp]", r"\centering", r"\small",
             r"\begin{tabular}{lrrrr}", r"\toprule",
             r"Study & $n$ & Mean (\%) & Median (\%) & LOO mean range (\%) \\",
             r"\midrule"]
    for name, row in studies.items():
        table.append(f"{name} & {row['geometries']} & {row['mean_percent']:.3f} & "
                     f"{row['median_percent']:.3f} & "
                     f"[{row['loo_min_percent']:.3f}, {row['loo_max_percent']:.3f}] " + r"\\")
    table.extend((r"\bottomrule", r"\end{tabular}",
                  r"\caption{Descriptive sensitivity of completed 60-second comparisons. "
                  r"The leave-one-geometry-out (LOO) range recomputes the equally weighted "
                  r"mean after each possible single omission. It is not a confidence interval "
                  r"or an alternative primary endpoint. All original geometries remain "
                  r"in the reported primary means.}",
                  r"\label{tab:geometry-sensitivity}", r"\end{table}"))
    (ROOT / "submission/tables/geometry_sensitivity.tex").write_text("\n".join(table) + "\n")
    print(json.dumps({name: {key: value for key, value in row.items()
                           if key not in ("geometry_effects", "leave_one_geometry_out")}
                      for name, row in studies.items()}, indent=2))
