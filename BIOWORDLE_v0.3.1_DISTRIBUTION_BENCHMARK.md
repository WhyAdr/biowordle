# BioWordle distribution benchmark

Measured 2026-11-14 against commit `83f05f96102abea5d2fcce0ed75dbf3535127f71`. Reproduce from the repository root with `node tools/distribution_benchmark.mjs 3 > distribution_benchmark_results.json`. The script reads the HTTP assets directly from that commit. It tests the working-tree standalone only if it matches that commit byte for byte. On Windows the existing Playwright Chromium installation may require setting `PLAYWRIGHT_BROWSERS_PATH` to its actual location, and Chromium loopback access may require the repository's documented sandbox bypass. Install the project's Playwright browser if absent. Raw observations are in `distribution_benchmark_results.json`.

This is a **local HTTP delivery simulation**, not a measurement of a deployed hosting provider: no deployed URL or production response headers were established. A local server supplied either uncompressed assets with `no-store`, gzip with `no-store`, or gzip with `Cache-Control: public, max-age=3600`. The standalone is opened using `file://`. Chromium used a 390 × 844 mobile viewport, device scale factor 2, 4× CPU slowdown, and simulated Fast 3G (1.6 Mbps down, 150 ms latency) for HTTP. File loading had no network throttle. Each cold observation used a new browser context. The cache result is a second navigation in the same context after a cold navigation and interactions; it is not a cold start. Three samples are too few to estimate tail latency or compare small interaction differences reliably.

| Configuration | Median initial ready (ms) | Median transferred bytes | Switch profile / search / submit guess median (ms) | Median JS heap used at ready (bytes) |
| --- | ---: | ---: | --- | ---: |
| HTTP raw, cold | 16,139 | 3,237,010 | 174 / 101 / 102 | 4,492,820 |
| HTTP gzip, cold | 2,338 | 373,504 | 168 / 120 / 104 | 4,492,820 |
| HTTP gzip + cache, cold | 2,345 | 373,504 | 147 / 91 / 132 | 4,492,820 |
| HTTP gzip + cache, repeat navigation | 146 | 0 | Not measured | Not measured |
| Standalone `file://`, cold | 757 | No network transfer | 109 / 103 / 106 | 4,450,756 |

Ready means navigation and the visible Chemle search input; interactions use the Learn switch, search for Xylene, and selection of the result. Transfer is the sum of browser Navigation/Resource Timing `transferSize`, including its accounting of response overhead; it is not a packet capture. Standalone Resource Timing can report local file bytes, which are **not** network transfer and are excluded above. Memory is the Chromium CDP `Runtime.getHeapUsage` used JavaScript heap at ready, **not total browser or process memory**. No external requests or page errors were observed. The samples are emulation on this machine, not measurements on a physical mobile device, and only Chemle interaction latency was exercised.

## Asset experiments

The script also serializes hypothetical subsets of the exact v0.3 data and compresses them. These are **size projections**, not working builds, timing tests, or additive on-wire totals. Current boot initializes and renders both games, so simply deleting a game's records or SVGs breaks required behavior. `svg_sha256` metadata is retained in the no-SVG projection.

| Data configuration | Raw bytes | gzip bytes | Brotli bytes |
| --- | ---: | ---: | ---: |
| Current `puzzles_data.js` | 3,100,664 | 341,122 | 193,909 |
| Remove SVG values | 514,379 | 69,529 | 50,639 |
| Chemle records and similarity, no Seqle records | 1,451,820 | 174,343 | 107,513 |
| Seqle records, no Chemle records or similarity | 1,659,216 | 169,301 | 90,200 |
| Full standalone HTML | 3,235,746 | 371,293 | 219,337 |

The gzip result demonstrates a substantial cold-network improvement for the unmodified application in this simulation. Brotli numbers are offline compression sizes at quality 11, not a server, browser, latency, or CPU benchmark. Serving gzip and Brotli should negotiate `Accept-Encoding` and send `Vary: Accept-Encoding`; deployment headers and behavior require separate verification. Caching can improve repeat visits, but a one-hour cache on mutable paths is unsafe for scientific cutovers: serve content-addressed, version-consistent assets with long-lived immutable caching and revalidate HTML, retaining older bundles for existing sessions. Verify rollback and mixed-version behavior before deploying.

Splitting games or lazily fetching SVGs may reduce initial bytes, but would require game-specific initialization, loading/error states, on-demand guess/scaffold/reveal assets, integrity/version checks, and tests for offline and mid-session behavior. The current SVG properties include duplicate peptide renderings, so investigate deduplication as well. Benchmark real functional builds and Seqle tab switching on physical mobile devices before choosing. Preserve the existing self-contained HTML generator and its offline E2E checks; it cannot depend on later network fetches. No application or hosting configuration was changed for these experiments.

## Working-tree measurement (v0.3.2 with Challenge Practice)

Re-measured 2026-11-14 against the updated working tree (`node tools/distribution_benchmark.mjs 3 --working`) under identical Chromium mobile emulation (390 × 844, 4× CPU; Fast 3G loopback for HTTP, unthrottled for `file://`):

| Configuration | Median initial ready (ms) | Median transferred bytes | Switch profile / search / submit guess median (ms) | Median JS heap used at ready (bytes) |
| --- | ---: | ---: | --- | ---: |
| HTTP raw, cold | 16,428 | 3,285,935 | 146 / 103 / 122 | 4,539,644 |
| HTTP gzip, cold | 2,328 | 378,743 | 123 / 89 / 99 | 4,539,644 |
| HTTP gzip + cache, cold | 2,386 | 378,743 | 132 / 90 / 132 | 4,539,644 |
| HTTP gzip + cache, repeat navigation | 172 | 0 | Not measured | Not measured |
| Standalone `file://`, cold | 903 | No network transfer | 139 / 106 / 106 | 4,496,776 |

Zero page errors and zero external requests were observed across all 12 samples (`3` runs × `4` modes). Working-tree asset sizes: `puzzles_data.js` = `3,133,804` B raw / `342,735` B gzip / `194,225` B Brotli; `dist/biowordle_standalone.html` = `3,284,671` B raw / `376,481` B gzip / `222,986` B Brotli.
