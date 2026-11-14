# BioWordle v0.4 research performance comparison

Baseline: checked-out `03febdf790c39b6d14202833af671af469146ba2` (78 Chemle records). Comparison: the isolated, candidate-conditional `tools/v04_compiled_preview.json` (120), whose exact SHA-256 is recorded in `tools/v04_perf_results.json`. It is a research preview built from unapproved candidate records, not a release. Production assets were not changed.

After generating the preview with the solver benchmark, reproduce the measurements from the repository root:

```powershell
$env:PYTHONUTF8='1'
.\.venv\Scripts\python.exe tools/v04_perf.py --preview tools/v04_compiled_preview.json --output tools/v04_perf_results.json --browser --runs 2
```

The saved preview JSON was normalized from CRLF to LF after measurement for cross-platform Git checkouts; the results record both the measured-file hash and the checked-in-file hash, and the parsed JSON is unchanged.

The browser step needs the installed project Playwright Chromium and permission for loopback and browser isolation. Without `--browser`, the same command runs the compiler, asset-size, golden-vector and Ranked checks. The baseline `build_all` timing includes its normal production gates. The 120 compilation is measured separately by the solver tool under its explicitly documented research gate-capture path, so timings are not strictly equivalent. Each is a single wall-clock sample on this machine, not a CI estimate.

| Metric | 78 baseline | 120 research preview |
|---|---:|---:|
| Compiler `build_all` wall time | 3.68 s | 5.22 s, research capture |
| Generated data JavaScript, raw | 3,133,800 B | 3,907,536 B |
| Generated data JavaScript, gzip level 9 | 342,731 B | 435,073 B |
| Self-contained standalone, raw | 3,284,629 B | 4,058,365 B |
| Self-contained standalone, gzip level 9 | 376,466 B | 468,771 B |
| Chemle golden vectors across three profiles | 18,252 | 43,200 |
| Golden JSON, raw | 1,037,622 B | 2,317,686 B |
| Golden generation | 0.130 s | 0.405 s |
| Full Node golden and DOM verifier | 0.631 s, pass | 0.721 s, pass |

Golden growth follows `profiles × N²`; these vectors are test artifacts and are not embedded in the standalone bundle. Sizes are regenerated using the existing compiler serialization and bundling format; gzip is offline level 9 and is not a deployed-server measurement. No full Playwright suite was run against the preview; its browser checks are smoke and latency measurements.

Chromium used a 390 × 844 mobile viewport, device scale 2, 4× CPU slowdown and, for local gzip HTTP, simulated Fast 3G at 1.6 Mbps down and 150 ms latency. `file://` is unthrottled. Each number below is the median of just **two** cold contexts; variations of this size cannot support fine timing conclusions. Ready means the Learn search input was visible; latency measures a search for the existing m-xylene record and submission. Transfer is browser Resource Timing `transferSize`, including its accounting overhead. Heap is CDP JavaScript heap at ready, not total browser/process memory. There were no observed page errors or external requests.

| Measure | 78 baseline | 120 research preview |
|---|---:|---:|
| gzip HTTP ready | 2,511 ms | 2,941 ms |
| gzip HTTP transfer | 378,724 B | 471,066 B |
| gzip HTTP search / submit | 143 / 179 ms | 151 / 140 ms |
| JS heap at ready | 4,471,098 B | 5,089,576 B |
| offline file ready | 1,154 ms | 1,396 ms |

Actual device performance, deployment compression/cache headers, variability and accessibility remain untested. Keep standalone E2E and test a functional 120 build on real mobile devices and the CI runner before release. Existing Challenge UI strings and E2E expectations explicitly mention 78 and require careful updates while preserving fixed baseline coverage.

The 36 Ranked answer and guess IDs match; all 1,296 Ranked guess/target feedback signatures and Ranked clue ladders match; and the Chemle schedule hash remains `9257bd168fa6`. **Nevertheless the all-catalog rules hash changes from `7edd1b8864c9` to `2c81484214e1`.** The current Ranked storage key includes that hash, so this preview would archive and restart existing sessions. Preserving Ranked schedules alone does not meet the historical-session requirement; a separately designed and tested profile-scoped rules identity and migration is required before implementation.
