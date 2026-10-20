# Reproducing the completed 300-second screen

This guide accompanies a partial manuscript. The completed evidence is a
selected TRAIN screen, not a completed ten-seed iteration-limited evaluation.
Read `../diario.md` and `../research_notes/CONTINUATION.md` for the live state
before launching any experiment. Paths below are relative to the repository
root; use branch `codex/pcctp-learning-submission`.

## Evidence and denominators

The eligible dataset has nine geometries, seeds 70–72 and three methods:
MDM300, MDM10+repair and MDM60+repair, all under 300 seconds total. The endpoint
dataset contains 81 final records and 27 environmentally accepted blocks.
H60 passes the original exploratory screen (+0.979441%, 8 wins/1 loss); H10
fails (−1.279343%, 6 wins/3 losses). The 12 reserved test geometries remain closed.

The first collection affected by host sleep is retained separately. In the
replacement collection, ts225/seed72 also failed clock-monitor thresholds;
all three attempts were archived and repeated together. The archived costs
do not enter the 81 eligible endpoints. The review manifest records why the
whole block was repeated, with hashes of the preserved files.

Important artifacts:

| Artifact | Role |
|---|---|
| `experiments/lkh/ruin_total300_awake_train/protocol.json` | Cases, seeds, arms, frozen sources/binaries and environmental limits |
| `experiments/lkh/ruin_total300_awake_train/g_*.seed*.{mdm,hybrid10,hybrid60}.json` | Eligible final records, including initialization and search traces |
| `experiments/lkh/ruin_total300_awake_train/environment.g_*.json` | Per-block clocks, sleep checks and result hashes |
| `experiments/lkh/ruin_total300_awake_train_audit.json` | Scalar costs/feasibility and endpoint eligibility audit |
| `experiments/lkh/ruin_total300_awake_environment_audit.json` | Independent accepted-block audit |
| `experiments/lkh/awake300_closing_analysis.json` | Exact rational reconstruction and descriptive geometry sensitivity |
| `experiments/lkh/awake300_completed_snapshots/` | Immutable copies of the completed upstream reports |
| `experiments/lkh/ruin_total300_awake_interrupted_blocks/` | Excluded failed-block evidence and review manifest |
| `submission/figures/awake300_effects.*` | Both arms and all paired seeds, with numerical verification/provenance |
| `submission/compilation_audit.json` | Current partial PDF, source hashes and visual review |

## Runtime used

The completed screen used Apple M2 arm64, 16 GiB RAM, the isolated corrected
MDM executable and its FPMax helpers, and LKH 3.0.6. The exact binaries are
identified by SHA256 in the protocol; the baseline includes the declared
decoder initialization and occupied-elite-count corrections. It is not
claimed byte-identical to the originally published program.

The local analysis runtime is `/usr/bin/python3` 3.9.6 with NumPy 1.24.3,
SciPy 1.13.0 and Matplotlib 3.8.4. HiGHS 1.13.1 was available for the separate
integer-repair investigations; MDM/LKH are the optimizers used in this
completed 300-second screen. Local pdfLaTeX/TeX Live 2026 produced the current
PDF. Render and font-library differences may change image/PDF bytes even
when numerical results and typeset content agree. SHA256 identifies the
recorded artifacts rather than guaranteeing bit-identical rerendering on
another machine.

## Read-only and derived-output commands

Run commands one at a time from the repository root. They launch no optimizer.
Observe which commands regenerate derived files before using them alongside
another assistant or the running queue.

```bash
/usr/bin/python3 experiments/lkh/verify_submission_manifest.py
/usr/bin/python3 experiments/lkh/calibration_progress.py
```

The first command checks the manifest scope and reports an external LKH
dependency separately. It reads no unlisted dataset files. The second updates
the passive nine-slot calibration progress report; checkpoint costs and
iterations remain provisional and do not become MDM70 endpoints.

After the complete-batch observers have finished, the following commands
reconstruct the closing analysis and render its figure:

```bash
/usr/bin/python3 experiments/lkh/close_awake300.py
MPLCONFIGDIR=/private/tmp/pcctp-awake300-mpl XDG_CACHE_HOME=/private/tmp/pcctp-awake300-fontcache /usr/bin/python3 experiments/lkh/plot_awake300.py
```

The closing analysis refuses incomplete batches and a different overwrite.
The figure has 18 mean bars and 54 paired-seed markers on one common axis;
markers are not confidence intervals. Its means use the reduction of
seed-mean costs, not the average of seedwise percentages. Omission ranges
are descriptive sensitivity measures and do not replace the primary screen.

`observe_awake300.py` also launches no search, but rewrites several reports.
Run it only when no other observer or queue is writing those outputs. A
successful completion marker is checked together with scalar, deadline,
environment and hash evidence, rather than being accepted alone.

## Commands that launch optimization

`awake300_train.py`, `published70_after_awake.py`,
`published70_calibration.py`, `paired70_train.py` and `cached30_train.py`
can launch optimizers. They are not audit commands. Check actual processes,
kernel locks and the diary before invoking them. In particular, the active
queue owns the seed80 calibration; invoking its worker directly could duplicate
work. Do not edit its frozen dependencies or executable binaries.

The calibration has nine TRAIN references, seed80, 70 iterations, positive time
budget disabled and a 43200-second safety limit interpreted as censoring.
It precedes the separately designed ten-seed comparison. The cache study is
conditional on 27 saved starts and 30 seconds per phase; its generation cost is
excluded and it has no real quality result yet. The shared paired70 kernel
lock prevents the two prospective studies from running together. Each still
requires its declared prerequisites and source-bound certificates.

## Current artifact scope

`reproducibility_manifest.json` binds the completed 300-second data, its inputs,
frozen algorithm sources/runtime binaries and current manuscript artifacts.
It excludes changing calibration checkpoints, user edits to `include/Utils.h`
and reserved test inputs. Separate 60-second and MIP studies retain their own
audits and selection histories. The manifest is an integrity record, not an
independent scientific replication or a claim that the manuscript is ready
for submission. It does not prescribe rerunning closed experiments while an
optimizer is active.

When a later audited manuscript revision changes a listed artifact, preserve
this manifest as historical provenance and issue a new reviewed version.
Recheck numerical and visual consistency before updating expected hashes;
changing a hash merely to suppress a mismatch is not verification.
