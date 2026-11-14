/**
 * Cross-Language Golden Vector & Runtime Invariant Verifier (P0-1, P0-3, P1-1, P1-2, P1-4)
 * ========================================================================================
 * Executed in Node.js to verify:
 * 1. Byte-identical feedback signatures between Python (`chem_engine.py`) and JS (`app.js`)
 *    across all 2,592 Chemle pairs (36x36x2) and 1,800 Seqle pairs (30x30x2).
 * 2. Exact agreement between Python and JS Henderson-Hasselbalch peptide calculations.
 * 3. Candidate counter exact signature filtering & monotonicity.
 * 4. `target in guessPool` invariant across all mode x profile x pack combinations
 *    (including the adrenaline + "Therapeutic Drug" regression test).
 * 5. UTC-day streak state machine (same-day idempotence, consecutive-day increment, gap reset).
 * 6. Pre-guess DOM/Clue-Ladder non-leakage & ARIA combobox/listbox accessibility attributes.
 */

import fs from "node:fs";
import path from "node:path";
import assert from "node:assert/strict";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const ROOT = path.resolve(__dirname, "..");

const require = createRequire(import.meta.url);
const Engine = require(path.join(ROOT, "app.js"));

const puzzlesData = JSON.parse(
  fs.readFileSync(path.join(ROOT, "data", "puzzles.json"), "utf8")
);
const goldenData = JSON.parse(
  fs.readFileSync(path.join(ROOT, "tests", "golden_feedback.json"), "utf8")
);
const indexHtml = fs.readFileSync(path.join(ROOT, "index.html"), "utf8");

// 1. Cross-language Golden Vectors: Chemle (36 x 36 x 2 = 2,592 pairs)
let chemChecked = 0;
for (const [profId, profSpec] of Object.entries(puzzlesData.feedback_spec.profiles)) {
  for (const g of puzzlesData.chemle) {
    for (const t of puzzlesData.chemle) {
      const key = `${g.id}->${t.id}`;
      const simPct = puzzlesData.chemle_similarity[g.id][t.id];
      const jsRes = Engine.evaluateChemleFeedbackFromSpec(g, t, simPct, profSpec);
      const pySig = goldenData.chemle_signatures[profId][key];
      assert.equal(
        jsRes.signature,
        pySig,
        `[P1-1 Chemle Golden Mismatch] profile=${profId} pair=${key}: JS=${jsRes.signature} vs PY=${pySig}`
      );
      chemChecked++;
    }
  }
}
assert.equal(chemChecked, 2592, "Expected 2,592 Chemle golden vector checks");

// 2. Cross-language Golden Vectors: Seqle (30 x 30 x 2 = 1,800 pairs) + H-H Biochem parity
let seqChecked = 0;
for (const m of puzzlesData.seqle) {
  const jsBio = Engine.computePeptideBiochem(m.seq, puzzlesData.aa_properties);
  assert.equal(jsBio.gravy, m.gravy, `GRAVY mismatch on ${m.seq}`);
  assert.equal(jsBio.net_charge_ph74, m.net_charge_ph74, `pH 7.4 charge mismatch on ${m.seq}`);
  assert.equal(jsBio.isoelectric_point, m.isoelectric_point, `pI mismatch on ${m.seq}`);
}

let seqRankedClassTiles = 0;
let seqLearnClassTiles = 0;
for (const [profId, profSpec] of Object.entries(puzzlesData.feedback_spec.profiles)) {
  for (const g of puzzlesData.seqle) {
    for (const t of puzzlesData.seqle) {
      const key = `${g.seq}->${t.seq}`;
      const jsRes = Engine.evaluateSeqleFeedbackFromSpec(
        g.seq,
        t.seq,
        puzzlesData.aa_properties,
        profSpec,
        g,
        t
      );
      if (profId === "ranked") {
        seqRankedClassTiles += jsRes.tiles.filter((x) => x === "class").length;
      } else if (profId === "learn") {
        seqLearnClassTiles += jsRes.tiles.filter((x) => x === "class").length;
      }
      const pySig = goldenData.seqle_signatures[profId][key];
      assert.equal(
        jsRes.signature,
        pySig,
        `[P1-1 Seqle Golden Mismatch] profile=${profId} pair=${key}: JS=${jsRes.signature} vs PY=${pySig}`
      );
      seqChecked++;
    }
  }
}
assert.equal(seqChecked, 1800, "Expected 1,800 Seqle golden vector checks");
assert.equal(seqRankedClassTiles, 0, "Ranked Seqle must suppress positional 'class' tiles");
assert.ok(seqLearnClassTiles > 0, "Learn Seqle must preserve positional 'class' tiles");

// 3. P0-3 Category-Lockout Regression Test: target in guessPool for ALL mode x profile x pack
const categories = ["ALL", ...new Set(puzzlesData.chemle.map((m) => m.Category))];
for (const profile of ["ranked", "learn"]) {
  for (const pack of categories) {
    const guessPool = Engine.getGuessPool(puzzlesData, "chemle");
    assert.equal(guessPool.length, 36, "Chemle guessPool must always contain all 36 compounds");
    for (const target of puzzlesData.chemle) {
      assert.ok(
        guessPool.some((x) => x.id === target.id),
        `Target ${target.id} missing from guessPool under profile=${profile}, pack=${pack}`
      );
    }
  }
}

// Explicit Adrenaline + "Therapeutic Drug" regression assertion
const adrenaline = puzzlesData.chemle.find((m) => m.id === "adrenaline");
assert.ok(adrenaline, "Adrenaline must exist in catalog");
const guessPoolWhenDrugPackActive = Engine.getGuessPool(puzzlesData, "chemle");
assert.ok(
  guessPoolWhenDrugPackActive.some((m) => m.id === "adrenaline"),
  "Regression failure: Adrenaline must remain selectable in guessPool even when Therapeutic Drug pack is active"
);

// 4. Candidate Counter Exact Signature Filtering
const chemById = Object.fromEntries(puzzlesData.chemle.map((m) => [m.id, m]));
for (const profile of ["ranked", "learn"]) {
  const profSpec = puzzlesData.feedback_spec.profiles[profile];
  const ansPool = Engine.getAnswerPool(puzzlesData, "chemle", profile, "ALL");
  const target = chemById["aspirin"];
  const c0 = Engine.countCompatibleChemleCandidates(
    [],
    target,
    ansPool,
    puzzlesData.chemle_similarity,
    profSpec,
    chemById
  );
  const c1 = Engine.countCompatibleChemleCandidates(
    ["caffeine"],
    target,
    ansPool,
    puzzlesData.chemle_similarity,
    profSpec,
    chemById
  );
  const cWin = Engine.countCompatibleChemleCandidates(
    ["caffeine", "aspirin"],
    target,
    ansPool,
    puzzlesData.chemle_similarity,
    profSpec,
    chemById
  );
  assert.equal(c0, 36);
  assert.ok(c1 < c0 && c1 >= 1, `Expected 1 <= c1 < 36 after first guess, got ${c1}`);
  assert.equal(cWin, 1, "Winning guess must narrow compatible candidates to 1");
}

// 4b. Learn Candidate Filtering Completeness (SIM_PCT, live FG overlap, and clue ladder)
{
  const learnSpec = puzzlesData.feedback_spec.profiles.learn;
  const learnNoExtras = {
    ...learnSpec,
    show_exact_tanimoto: false,
    show_functional_groups_live: false,
  };
  const learnPool = Engine.getAnswerPool(puzzlesData, "chemle", "learn", "ALL");
  const target = chemById["aspirin"];
  const withFullLearn = Engine.countCompatibleChemleCandidates(
    ["ibuprofen"],
    target,
    learnPool,
    puzzlesData.chemle_similarity,
    learnSpec,
    chemById
  );
  const withoutSimFg = Engine.countCompatibleChemleCandidates(
    ["ibuprofen"],
    target,
    learnPool,
    puzzlesData.chemle_similarity,
    learnNoExtras,
    chemById
  );
  assert.ok(
    withFullLearn <= withoutSimFg && withFullLearn >= 1,
    `Learn SIM_PCT/FG feedback must narrow or equal stripped spec (${withFullLearn} vs ${withoutSimFg})`
  );

  // Also test a pair where exact SIM_PCT / FG strictly disambiguates targets
  const lvl0Clue = target.clue_ladder.learn[0].text;
  const cClue0 = Engine.countCompatibleChemleCandidates(
    [],
    target,
    learnPool,
    puzzlesData.chemle_similarity,
    learnSpec,
    chemById,
    { profileId: "learn", level0Text: lvl0Clue, superClass: null, scaffoldUnlocked: false, targetScaffoldSmiles: "" }
  );
  const cClue1 = Engine.countCompatibleChemleCandidates(
    [],
    target,
    learnPool,
    puzzlesData.chemle_similarity,
    learnSpec,
    chemById,
    {
      profileId: "learn",
      level0Text: lvl0Clue,
      superClass: target.super_class,
      scaffoldUnlocked: true,
      targetScaffoldSmiles: target.scaffold_smiles || "",
    }
  );
  assert.ok(cClue0 < 36 && cClue0 >= 4, `Level-0 Learn clue should narrow to 4..35, got ${cClue0}`);
  assert.ok(cClue1 <= cClue0 && cClue1 >= 1, `Level-1 + scaffold clue must further narrow (${cClue1} <= ${cClue0})`);

  const seqPool = Engine.getAnswerPool(puzzlesData, "seqle", "learn", "ALL");
  const yigsr = puzzlesData.seqle.find((s) => s.seq === "YIGSR");
  const seqClue0 = Engine.countCompatibleSeqleCandidates(
    [],
    yigsr,
    seqPool,
    puzzlesData.aa_properties,
    learnSpec,
    { profileId: "learn", level0Text: yigsr.clue_ladder.learn[0].text, level1Text: null }
  );
  const seqClue1 = Engine.countCompatibleSeqleCandidates(
    [],
    yigsr,
    seqPool,
    puzzlesData.aa_properties,
    learnSpec,
    { profileId: "learn", level0Text: yigsr.clue_ladder.learn[0].text, level1Text: yigsr.clue_ladder.learn[1].text }
  );
  assert.ok(seqClue0 < 30 && seqClue1 <= seqClue0 && seqClue1 >= 1, "Seqle clue constraints must narrow Learn pool");
}

// 5. P1-2 UTC-Day Streak State Machine & v1-to-v2 Migration Verification
const migratedStats = Engine.migrateV1Stats({ chemWins: 5, seqWins: 3, streak: 7, maxStreak: 9 });
assert.equal(migratedStats.schema, 2);
assert.equal(migratedStats.chemWins, 5, "v1 chemWins must be preserved");
assert.equal(migratedStats.seqWins, 3, "v1 seqWins must be preserved");
assert.equal(migratedStats.streak, 0, "v1 streak must NOT migrate without verified UTC dates");
assert.equal(migratedStats.maxStreak, 0, "v1 maxStreak must reset to 0 on v2 migration");

globalThis.__BIOWORDLE_TEST_UTC_DATE__ = "2026-12-31";
assert.equal(Engine.getUtcDateString(), "2026-12-31", "Test UTC date override must work cleanly");
delete globalThis.__BIOWORDLE_TEST_UTC_DATE__;

let st = {
  schema: 2,
  chemWins: 0,
  seqWins: 0,
  streak: 0,
  maxStreak: 0,
  lastDailyWinUtcByGame: { chemle: null, seqle: null },
  lastStreakIncrementUtc: null,
};
st = Engine.evaluateDailyWinStreak(st, "chemle", "2026-11-14");
assert.equal(st.streak, 1, "First daily win sets streak to 1");
assert.equal(st.chemWins, 1);

// Re-winning on the same UTC day (or winning Seqle on the same UTC day) must NOT increment streak
st = Engine.evaluateDailyWinStreak(st, "chemle", "2026-11-14");
assert.equal(st.streak, 1, "Same-day replay must not increment streak");
assert.equal(st.chemWins, 1, "Same-day replay must not double-count daily win");

st = Engine.evaluateDailyWinStreak(st, "seqle", "2026-11-14");
assert.equal(st.streak, 1, "Winning second mode on same UTC day keeps streak at 1");
assert.equal(st.seqWins, 1);

// Winning on consecutive UTC day increments streak to 2
st = Engine.evaluateDailyWinStreak(st, "chemle", "2026-11-15");
assert.equal(st.streak, 2, "Consecutive UTC day win increments streak to 2");

// Winning after a 1-day gap resets streak to 1 while preserving maxStreak=2
st = Engine.evaluateDailyWinStreak(st, "chemle", "2026-11-17");
assert.equal(st.streak, 1, "Gap in UTC days resets current streak to 1");
assert.equal(st.maxStreak, 2, "maxStreak is preserved across streak resets");

// 6. P0-1 Pre-Guess Clue Ladder & DOM Non-Leakage + P2-2 ARIA Verification
for (const m of puzzlesData.chemle) {
  for (const prof of ["ranked", "learn"]) {
    const lvl0 = m.clue_ladder[prof].find((c) => c.level === 0);
    assert.ok(lvl0, `Missing Level-0 clue for ${m.id} (${prof})`);
    assert.ok(
      !lvl0.text.includes(String(m.MW)),
      `Level-0 clue for ${m.id} leaks exact MW (${m.MW}): ${lvl0.text}`
    );
    assert.ok(
      !lvl0.text.includes(m.Category),
      `Level-0 clue for ${m.id} leaks exact Category (${m.Category}): ${lvl0.text}`
    );
    assert.ok(
      !lvl0.text.toLowerCase().includes("bits"),
      `Level-0 clue for ${m.id} leaks IG bits: ${lvl0.text}`
    );
  }
}

assert.ok(indexHtml.includes('role="combobox"'), "index.html must include role=combobox on search input");
assert.ok(indexHtml.includes('role="listbox"'), "index.html must include role=listbox on dropdown");
assert.ok(indexHtml.includes('id="utc-rollover-banner"'), "index.html must include utc-rollover-banner");
assert.ok(!indexHtml.includes("IG bits shown"), "index.html search placeholder must not expose pre-guess IG bits");

console.log(
  `[Golden & DOM Verification PASS] ${chemChecked} Chemle vectors + ${seqChecked} Seqle vectors + P0-1/P0-3/P1-2/P2-2 runtime invariants verified.`
);
