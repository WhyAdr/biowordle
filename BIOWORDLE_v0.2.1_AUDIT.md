# BioWordle v0.2.1 independent release audit

**Reviewed:** `2a0797e06ac53a75f9ac0e7a1a3c385336d5b50e` (2026-11-14). The checkout and local `origin/main` both resolved to this commit, four commits after `a30fb92`. The release owner confirmed the live remote and Actions results. This audit reran the release checks locally. The companion [v0.2.2 patch plan](BIOWORDLE_v0.2.2_PATCH_PLAN.md) maps each finding to a change and acceptance tests.

## Severity-ranked findings

| ID | Severity | Finding |
| --- | --- | --- |
| F1 | High | Ranked Seqle's short solves are a catalog and information-budget problem. The current solver gate accepts 26 of 28 solves within two guesses and omits legal free-text probes. |
| F2 | High | Two already-open tabs can discard a recorded win. |
| F3 | High | Opening structure wells show target-dependent silhouettes. Seqle puts the full structure in the pre-guess DOM, and Chemle puts the gated scaffold there. |
| F4 | High | Two Seqle structure citations misrepresent their evidence, one exemplar teaches incorrect terminal spacing, and the compiler accepts new Ranked motifs with no provenance entry. |
| F5 | Medium | The browser suite passes on Windows Edge and is absent from the Ubuntu workflow. Its setup, random branches and UTC test make a straight CI transplant fragile. |
| F6 | Medium | The v0.2.1 legacy persistence fallback misses the actual v0.2.0 keys. The transition can orphan a same-day game or an existing classroom seed. |
| F7 | Medium | The reported opening solver metric ignores opening clues, and solver/UI clue observations can diverge after later guesses. |
| F8 | Low | The rebuild job compares only a subset hash after overwriting shipped assets, which leaves stale checked-in artifacts outside its assertion. |

## F1 · Ranked Seqle solves

The local solver and rebuilt metadata agree. The 28 Ranked targets split into charge groups of **8, 8 and 12**. Greedy entropy over the 30 catalog guesses solves them in **2.00** guesses on average, median **2**, p95 **3**, maximum **3**. The histogram is **2 at 1, 24 at 2, 2 at 3**. The gate requires mean ≥1.90 and fewer than 95% within two guesses. The reported maximum of 21 distinct signatures is over the whole pool. A player can already use the opening charge group to select a probe.

### Controlled benchmarks

The table holds the 28 targets, the catalog probe universe and tie breaking constant. Each row is a separate variant using the shipped greedy solver, with the charge and architecture clues available as shipped except where stated. Distribution columns count targets.

| Ranked variant | 1 | 2 | 3 | Mean | ≤2 guesses |
| --- | ---: | ---: | ---: | ---: | ---: |
| Shipped | 2 | 24 | 2 | 2.00 | 92.86% |
| Remove Q and GRAVY comparison strips | 2 | 24 | 2 | 2.00 | 92.86% |
| Delay opening charge to after guess 2 | 1 | 20 | 7 | 2.21 | 75.00% |
| Enable positional biochemical class tiles | 3 | 25 | 0 | 1.89 | 100.00% |
| Exact positions plus total misplaced count, with strips | 2 | 22 | 4 | 2.07 | 85.71% |
| Exact and misplaced totals only, with strips | 3 | 20 | 5 | 2.07 | 82.14% |
| Exact and misplaced totals only, without strips | 3 | 18 | 7 | 2.14 | 75.00% |

Reproduce from [the shipped evaluator and solver](chem_engine.py#L658) with `profiles.ranked` from [the feedback spec](data/feedback_spec.json). Remove strips with `seqle_strip_fields=[]`. Enable class via `show_positional_class_match=True`. For the delay, set Ranked level-0 text to the same constant on every record and append its former text to the level-1 text, preserving the level-1 unlock at two guesses. The shipped fully unconditioned run happened to yield the same distribution. For the count variants, replace only the tile segment of each evaluator signature with either `exact_count:present_count` or `exact_position_bitmask:present_count`. Preserve the evaluator's duplicate handling and comparison strips, except in the stated ablation. Call `run_seqle_greedy_solver(records, aa_properties, spec, condition_on_clues=True)`; its deterministic tie break is rounded entropy, then candidate membership, then lexical maximum.

**Legal probe diagnostic.** The UI accepts arbitrary valid five-residue guesses ([submit handler](app.js#L1783)). The catalog-only solver does not. With Ranked class tiles mapped to absent and the strips removed, best catalog probe `LPETG` gives 18 outcome buckets and an expected 62/28 remaining candidates across the 28 targets. Adding 10,000 unique random legal probes yielded a best result of 21 buckets and 46/28 remaining. In the opening groups of 12, 8 and 8, the best bucket counts were 10/12, 8/8 and 8/8 using catalog probes, versus **12/12, 8/8 and 8/8** with the sample. Those partitions establish a two-guess strategy for the current target pool given the opening charge group, even without strips. This is a lower bound on the quality of legal probes, from a sample rather than an exhaustive vocabulary search.

```python
import json, random
from collections import Counter
from chem_engine import evaluate_seqle_guess

d = json.load(open("data/puzzles.json", encoding="utf-8"))
aa = d["aa_properties"]
ranked = [r for r in d["seqle"] if "ranked" in r["game_eligibility"]]
catalog = {r["seq"] for r in d["seqle"]}
groups = {}
for r in ranked:
    groups.setdefault(r["clue_ladder"]["ranked"][0]["text"], []).append(r["seq"])
rng = random.Random(20261114)
probes = set()
while len(probes) < 10_000:
    probes.add("".join(rng.choices(sorted(aa), k=5)))

def score(g, group):
    buckets = Counter("".join("A" if t == "class" else t[0].upper()
                              for t in evaluate_seqle_guess(g, target, aa))
                      for target in group)
    return sum(n*n for n in buckets.values()), -len(buckets), g

for group in [[r["seq"] for r in ranked], *groups.values()]:
    for vocabulary in [catalog, catalog | probes]:
        print(len(group), min(score(g, group) for g in vocabulary))
```

### Scientifically meaningful design options

1. **Expand the answer pool with documented natural instances.** Keep the familiar positional feedback and introduce more variation within experimentally supported motif families as well as across families. A matched set of naturally occurring variants makes a local residue observation useful without immediately identifying a unique family representative. Screen for duplicate 5-mers, functional evidence, precursor/isoform coordinates, processing or terminal context, and opening clue bucket sizes. Compare candidate catalogs of 60, 90 and 120 targets against the legal-probe benchmark. These quantities are experiment points. Expansion has no measured depth gain until an independently curated corpus is simulated. Avoid filling the catalog with arbitrary protein windows or one representative per maximally different family.
2. **Stage biochemical context as a hint.** Moving charge to guess 2 costs little and gives a measured +0.21 guesses, but every target still solves within three in the catalog benchmark. On a larger corpus, test staged charge, followed by terminal processing or functional-class clues that retain several documented candidates. Late biochemical cues can teach interpretation and rescue hard runs. Earlier or more granular cues further partition the pool, as the class-tile ablation demonstrates. Keep the full comparative cues in Learn. Label modeled Q and pI as estimates for an isolated peptide with free ends, distinct from properties of a motif embedded in a protein or chemically modified.
3. **Offer an aggregate feedback experiment only if it serves a teaching goal.** Exact positions plus a misplaced-residue total retains positional anchors and some compositional reasoning. The measured gain is +0.07 guesses. Count-only feedback trades away positional reasoning for a small gain. Preserve ordinary positional feedback for the first expansion prototype. Broader free-text guess access by itself can shorten solves.

## F2 · Shared stats can lose a win across tabs

On the same HTTP origin in one browser context, open tabs A and B before either wins. Set `window.__BIOWORDLE_TEST_UTC_DATE__ = "2026-11-14"` in an init script before each tab boots. Solve Ranked Chemle with Artemisinin in A, then solve Ranked Seqle with `SYSGY` in B. Inspect `JSON.parse(localStorage.getItem("biowordle_stats_v2"))` after each win. Stored totals changed from `{chemWins:1, seqWins:0}` to `{chemWins:0, seqWins:1}`. Reloading both tabs retained only the Seqle win. The submit handlers mutate cached stats and replace the shared record with a complete snapshot. The game-progress key and stats key are separate writes, so same-puzzle competing edits and interruption between writes also warrant regression tests. [Stats read/write](app.js#L507), [Chemle win](app.js#L1314), [Seqle win](app.js#L1830).

## F3 · Rendered opening leaks structure beyond the calibrated clue

In a real browser at UTC `2026-11-14`, the opening Chemle well contains 44 SVG paths with the same path geometry as Artemisinin's target scaffold. The opening Seqle well contains 104 paths matching the entire `SYSGY` target SVG. Both wells say Level 0 and CSS visually blurs the target-dependent geometry. With the same date fixture as F2, inspect `#chemle-structure-well svg path` and `#seqle-structure-well svg path` before submitting. Compare the ordered `d` attributes with the corresponding asset parsed into a detached DOM element. This reproduces an opening-design defect in the shipped build even with the v0.2.1 hint button lockout working. The standalone catalog and deterministic daily chooser remain public by design. Use a target-independent pre-unlock well to make the rendered clue consistent with the modeled one. [Chemle well](app.js#L1051), [Seqle well](app.js#L1674), [CSS filter](styles.css#L268).

## F4 · Biological provenance and teaching claims

**All 24 `verified_uniprot` sequences passed** the independent FASTA, parent range, motif range and offset checks using UniProtKB release **2026_04**. The release data's registry has 30 entries. I checked the cited RCSB entries at the deposited polymer-entity sequence level, including unmapped short peptides, and checked the four consensus examples separately below. Deposited sequence membership need not establish ordered coordinates or biological function. Direct ELM class pages timed out, so the statements about recognition rules below use the cited structures and primary experiments. Other mechanistic prose outside the identified issues received spot checks.

### Native sequence register

Coordinates are 1-based and inclusive. In the structure column, `entity:offset` is the first 1-based occurrence in the named deposited entity's canonical sequence. Sequence confirmation and same-protein mapping are distinct checks.

| Five-mer | UniProtKB | Parent range | Motif range | RCSB evidence |
| --- | --- | --- | --- | --- |
| GHSMG | [G9BY57](https://www.uniprot.org/uniprotkb/G9BY57/entry) | 160–170 | 163–167 | [4EB0](https://www.rcsb.org/structure/4EB0) 1:128, mapped |
| EKDEL | [P11021](https://www.uniprot.org/uniprotkb/P11021/entry) | 648–654 | 650–654 | [6I6H](https://www.rcsb.org/structure/6I6H) 2:2, exact peptide contained in receptor-bound AEKDEL |
| DEVDG | [P09874](https://www.uniprot.org/uniprotkb/P09874/entry) | 210–219 | 211–215 | No PDB cited |
| GRGDS | [P02751](https://www.uniprot.org/uniprotkb/P02751/entry) | 1612–1622 | 1614–1618 | [1FNF](https://www.rcsb.org/structure/1FNF) 1:351, mapped |
| LPETG | [P02976](https://www.uniprot.org/uniprotkb/P02976/entry) | 480–489 | 482–486 | [1T2W](https://www.rcsb.org/structure/1T2W) 2:1, exact substrate peptide |
| HEIGH | [P08254](https://www.uniprot.org/uniprotkb/P08254/entry) | 215–226 | 218–222 | [1B3D](https://www.rcsb.org/structure/1B3D) 1:119, mapped |
| YIGSR | [P07942](https://www.uniprot.org/uniprotkb/P07942/entry) | 946–956 | 950–954 | No PDB cited |
| KCVIM | [P01116-2](https://www.uniprot.org/uniprotkb/P01116-2/entry) | 175–188 | 184–188 | [1D8D](https://www.rcsb.org/structure/1D8D) 3:7, KRAS4B peptide, RCSB maps to base P01116 |
| DEADE | [P60842](https://www.uniprot.org/uniprotkb/P60842/entry) | 180–190 | 182–186 | [2ZU6](https://www.rcsb.org/structure/2ZU6) 1:164, mapped |
| GVGKS | [P01116](https://www.uniprot.org/uniprotkb/P01116/entry) | 10–17 | 13–17 | [4OBE](https://www.rcsb.org/structure/4OBE) 1:14, mapped |
| VDFGL | [P41279](https://www.uniprot.org/uniprotkb/P41279/entry) | 265–277 | 269–273 | [4Y83](https://www.rcsb.org/structure/4Y83) 1:206, mapped |
| HRDLA | [P10721](https://www.uniprot.org/uniprotkb/P10721/entry) | 789–797 | 790–794 | [1PKG](https://www.rcsb.org/structure/1PKG) 1:184, mapped |
| WSEWS | [P32927](https://www.uniprot.org/uniprotkb/P32927/entry) | 421–433 | 425–429 | [1EGJ](https://www.rcsb.org/structure/1EGJ) 1:88, mapped |
| PPPVP | [Q07889](https://www.uniprot.org/uniprotkb/Q07889/entry) | 1149–1158 | 1150–1154 | [1GBR](https://www.rcsb.org/structure/1GBR) entity 2 is mouse Sos1 Q02384 and lacks this window |
| INFEK | [P01012](https://www.uniprot.org/uniprotkb/P01012/entry) | 258–265 | 260–264 | [1VAC](https://www.rcsb.org/structure/1VAC) 3:3, mapped. The displayed source at `build_puzzles.py:683` incorrectly gives 259–263 |
| YQRLN | [P19814](https://www.uniprot.org/uniprotkb/P19814/entry) | 348–357 | 350–354 | [1BXX](https://www.rcsb.org/structure/1BXX) 2:2, exact TGN38 peptide window |
| CAACH | [P07143](https://www.uniprot.org/uniprotkb/P07143/entry) | 97–109 | 101–105 | [3CX5](https://www.rcsb.org/structure/3CX5) 4:40, mapped |
| DSGLD | [P25963](https://www.uniprot.org/uniprotkb/P25963/entry) | 27–39 | 31–35 | [1P22](https://www.rcsb.org/structure/1P22) entity 3 is β-catenin. [9T9W](https://www.rcsb.org/structure/9T9W) entity 2:4 maps to P25963 and contains the IκBα window |
| KKKRK | [P03070](https://www.uniprot.org/uniprotkb/P03070/entry) | 126–132 | 127–131 | [1BK6](https://www.rcsb.org/structure/1BK6) 2:1, exact peptide window |
| FGGFT | [Q13519](https://www.uniprot.org/uniprotkb/Q13519/entry) | 130–146 | 130–134 | No PDB cited. Coordinates use the precursor |
| YGGFM | [P01210](https://www.uniprot.org/uniprotkb/P01210/entry) | 100–104 | 100–104 | [1PLW](https://www.rcsb.org/structure/1PLW) 1:1, mapped peptide. The 5-mer occurs elsewhere in the precursor too |
| SYSGY | [P35637](https://www.uniprot.org/uniprotkb/P35637/entry) | 34–43 | 37–41 | [5W3N](https://www.rcsb.org/structure/5W3N) 1:64, mapped |
| NPLYK | [P05106](https://www.uniprot.org/uniprotkb/P05106/entry) | 763–788 | 770–774 | [1MK7](https://www.rcsb.org/structure/1MK7) 1:10, mapped |
| GHHPH | [P04196](https://www.uniprot.org/uniprotkb/P04196/entry) | 383–392 | 383–387 | No PDB cited. Multiple exact occurrences in this protein |

### Consensus exemplars and educational abstractions

| Record | Verified observations and correction |
| --- | --- |
| `NPVYK` | Human LDLR P01130 carries `NPVYQ` at 825–829. The RCSB 3SO6 entity 2 contains this LDLR window. The structural paper supports an FxNPxY internalization signal with downstream context in LDLR. The evidence supplied does not support a required Lys/Arg immediately after Tyr or an observed `NPVYK` example. Teach the supported core and show `NPVYK` as a pedagogical sequence, or use the verified LDLR window in a versioned pool. [UniProt P01130](https://www.uniprot.org/uniprotkb/P01130/entry), [RCSB 3SO6](https://www.rcsb.org/structure/3SO6) |
| `KKMPQ` | RCSB 4J77 entity 1 contains `AKEKSD`, a C-terminal KxKxx example. For the canonical COPI motifs, lysines occupy positions −4/−3 (KKxx) or −5/−3 (KxKxx) relative to the C terminus, with the final residue counted as −1. The shipped KKxxx sequence places lysines at −5/−4 and misses −3. This is an incorrect representative of the cited rule, so simply relabeling it as synthetic is insufficient. Use a validated terminal window, such as the five-residue suffix of the co-crystal peptide, after checking source identity and assay context. [RCSB 4J77](https://www.rcsb.org/structure/4J77), [primary coatomer study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3616288/) |
| `NLTGY` | The first three residues fit the usual N–non-Pro–S/T sequon. The appended Y is part of the game five-mer, and the shipped string (N-X-[ST]-X-Y) drops the non-Pro constraint while imposing a Tyr requirement. RCSB 3WAK is apo archaeal AglB O29867, with no acceptor peptide. The UniProt O29867 sequence also lacks this exact five-mer. Use a correctly stated sequon and context-only structural labeling, or curate an experimentally mapped eukaryotic site. [RCSB 3WAK](https://www.rcsb.org/structure/3WAK), [UniProt O29867](https://www.uniprot.org/uniprotkb/O29867/entry), [primary eukaryotic sequon experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC3047074/) |
| `KENIP` | Q08981 Acm1 carries `KENLS` at 98–102, also present in 4BH6 entity 2. P34244 Hsl1 has `NKENE` at 774–778. Both support the KEN core, neither supports a universal downstream Pro, and Acm1 supplies pseudosubstrate binding context. Use the KEN rule with a separately labeled teaching five-mer, or use a directly documented substrate window. [RCSB 4BH6](https://www.rcsb.org/structure/4BH6), [UniProt Q08981](https://www.uniprot.org/uniprotkb/Q08981/entry), [UniProt P34244](https://www.uniprot.org/uniprotkb/P34244/entry), [primary functional study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2493364/) |
| `CCHHC` | Learn-only pedagogical abstraction. PROSITE PS50158 and RCSB 1MFS give domain context. No deposited polymer sequence has the registered contiguous five-mer. Keep the abstraction label in the reveal view. [PROSITE PS50158](https://prosite.expasy.org/PS50158), [RCSB 1MFS](https://www.rcsb.org/structure/1MFS) |
| `LDLAS` | Learn-only pedagogical abstraction. RCSB 3NBY context is nuclear export. None of the three returned polymer sequences has the exact five-mer. Use as context evidence and retain the abstraction label. [RCSB 3NBY](https://www.rcsb.org/structure/3NBY) |

**Compiler and runtime.** In [the validator](chem_engine.py#L1223), missing registry keys bypass the provenance branch even with enforcement enabled. Reproduced by copying the `KENIP` raw record and setting `seq='KENIA', parent_sequence='KENIA', consensus='K-E-N-X-X', exemplar_source='UniProt:NOTREAL'`; the ranked-eligible record passed validation. The existing tests check local assertions, which is why their pass is compatible with the RCSB findings. [The reveal renderer](app.js#L1718) generates badges from `kind` alone, and renders “Parent Sequence Window” for all four consensus records even though their status identifies them as teaching examples. Use explicit evidence roles for native windows, exact experimental peptides, motif-family structures and abstractions. Remove the misleading bound-native claim in the DSGLD clue when citing 1P22, or cite the verified IκBα structure 9T9W.

**Reproduce database checks.** Extract `SEQLE_REFERENCE_REGISTRY` from [chem_engine.py](chem_engine.py#L827). For each native entry retrieve `https://rest.uniprot.org/uniprotkb/{accession}.fasta` using the suffixed accession for isoform 2B. Join sequence lines and assert both `seq[parent_start-1:parent_end] == verified_parent_sequence` and `seq[motif_start-1:motif_end] == registry_key`, then check `motif_start-parent_start == window_offset`. Record the UniProt release response header. For each PDB, fetch `https://data.rcsb.org/rest/v1/core/entry/{pdb}` to enumerate `polymer_entity_ids`. Fetch `polymer_entity/{pdb}/{id}` on the same base URL and inspect the whitespace-stripped `entity_poly.pdbx_seq_one_letter_code_can`, `uniprot_ids` and author chain IDs. Use a chain-specific deposited residue map or mmCIF for any eventual claim about observed atoms or alternative numbering.

## F5 · Portable browser coverage is missing

On Windows, `node tests/test_playwright_e2e.mjs` passed all five scenarios plus standalone smoke under installed Edge. The script requires Antigravity's private Playwright path, `LOCALAPPDATA`, the Edge channel and `file://`. Ubuntu CI has a Python and Node verify job, and no browser step. Current Learn filtering checks allow loose inequalities and random retry behavior. The UTC case changes the date override before several automated interactions and can let the interval fire first. It asserts consecutive targets differ, which the schedule does not guarantee. Use a lockfile, managed Chromium, localhost HTTP and controlled randomness/time. Add one standalone `file://` smoke to retain the advertised offline path. [Current E2E](tests/test_playwright_e2e.mjs), [current CI](.github/workflows/ci.yml), [Playwright CI](https://playwright.dev/docs/ci), [Playwright Clock](https://playwright.dev/docs/clock).

## F6 · Legacy key fallback uses the wrong version

At `a30fb92`, v0.2.0 stored daily and seed progress under dataset hash `919ad8afb190`. v0.2.1 first checks the schedule hash `9257bd168fa6`, then its **own** dataset hash `8bb96ac156cc`. The actual v0.2.0 suffix is absent from the lookup. The later target-ID check also prevents a blind key migration when the answer changed. To reproduce, load `app.js` and `data/puzzles.json` independently using `git show <commit>:<path>` and use each release's exported `getAnswerPool` and `selectDailyIndex` on its own data. For `selectDailyIndex("chemle", "ranked", "2026-11-14", versionHash, pool.length, seed)`, pass the old dataset hash or current schedule hash and use `null` or `"CHEM101-WEEK3"`. Daily Chemle changes from `nicotine` to `artemisinin`. The classroom seed changes from `penicillin_g` to `pyocyanin`. The schedule hash includes shared schema version, all Seqle records and their kind/eligibility, so a Learn-only change can invalidate both games' keys. The compiler sorts records today and the browser indexes in payload order. [Storage fallback](app.js#L539), [target validation](app.js#L842), [schedule hash](build_puzzles.py#L994).

## F7 · Solver and actual observations can diverge

Seqle conditions on level 0 and at most level 1. Ranked level 2 and 3 become visible at guesses 4 and 5. Chemle simulates opening text, superclass and scaffold, without the other visible clue-text filters. The gates' “min opening expected remaining” comes from partitions over the full answer pool. For Seqle Ranked, it stays 1.64 even with conditioning disabled. Finally, both solvers stop at an internal safety cap of 12; the common summary knows ≤6 but does not distinguish failure from a solve step at the cap or state the Learn success rate at its actual limit of 8. These mismatches matter when calibrating a deeper pool. They do not invalidate the reproduced short v0.2.1 Seqle histogram. [Seqle simulator](chem_engine.py#L719), [Chemle simulator](chem_engine.py#L479), [summary](chem_engine.py#L511), [runtime clue gating](app.js#L1018).

## F8 · Checked-in artifacts and user-facing release copy

I rebuilt in an isolated scratch directory using Python 3.12.10 and pinned RDKit 2026.09.2. JSON, browser data JS, golden vectors and standalone HTML matched the four tracked files byte for byte. CI currently checks a dataset hash after the build, which does not cover e.g. a stale inline runtime in the previously committed standalone HTML. README also still presents the old v0.2 thresholds as if active in the historical section, claims Learn packs generally despite the Seqle answer pool using only profile eligibility, and its v0.2.1 unlock summary obscures Seqle's level 1 after guess 2. Align release docs with per-game behavior when adding portable tests. [Build outputs](build_puzzles.py#L1379), [CI rebuild](.github/workflows/ci.yml#L30), [README](README.md), [Seqle answer pool](app.js#L292).

## Reproduction summary

- `git rev-parse HEAD origin/main` returned the pinned commit twice. Source review and local tests started from a clean tracked checkout.
- `& .\.venv\Scripts\python.exe tools/solver_bot.py` passed four gates. `& .\.venv\Scripts\python.exe -B -m unittest discover -v -s . -p 'test_*.py'` passed 11 tests. Use the pinned [requirements](requirements.txt) in a clean environment.
- `node tests/verify_golden_and_dom.mjs` independently passed 2,592 Chemle and 1,800 Seqle vectors plus runtime invariants. `node tests/test_playwright_e2e.mjs` passed five real-browser flows and standalone smoke under Edge.
- Linux Chromium CI, an expanded provenance-curated catalog, and human solve-depth measurements are acceptance work for the patch cycle.
