const fs = require("node:fs");
const path = require("node:path");
const { pathToFileURL } = require("node:url");
const { test, expect } = require("@playwright/test");

const ROOT = path.resolve(__dirname, "..", "..");
const BioWordleEngine = require("../../app.js");
const PUZZLES = JSON.parse(fs.readFileSync(path.join(ROOT, "data", "puzzles.json"), "utf8"));
const CHEM_SCHEDULE_HASH = PUZZLES.metadata.chemle_schedule_hash;
const SEQ_SCHEDULE_HASH = PUZZLES.metadata.seqle_schedule_hash;

function dailyTargetForDate(game, utcDate, customSeed = null) {
  const pool = BioWordleEngine.getAnswerPool(PUZZLES, game, "ranked", "ALL");
  const h = game === "chemle" ? CHEM_SCHEDULE_HASH : SEQ_SCHEDULE_HASH;
  const idx = BioWordleEngine.selectDailyIndex(game, "ranked", utcDate, h, pool.length, customSeed);
  return pool[idx];
}

function findRepeatedDailyTargetDates(game) {
  const start = new Date("2026-01-01T00:00:00Z");
  let prevDate = null;
  let prevKey = null;
  for (let i = 0; i < 400; i++) {
    const d = new Date(start.getTime() + i * 86400_000).toISOString().slice(0, 10);
    const t = dailyTargetForDate(game, d);
    const key = game === "chemle" ? t.id : t.seq;
    if (prevKey && key === prevKey) {
      return { dateA: prevDate, dateB: d, target: t };
    }
    prevDate = d;
    prevKey = key;
  }
  throw new Error(`No adjacent repeated ${game} daily target found in 400 days`);
}

function findOpeningGroupPair(game) {
  const pool = BioWordleEngine.getAnswerPool(PUZZLES, game, "ranked", "ALL");
  const byClue = new Map();
  for (const r of pool) {
    const lvl0 = r.clue_ladder.ranked.find((c) => c.level === 0)?.text || "";
    if (!byClue.has(lvl0)) byClue.set(lvl0, []);
    byClue.get(lvl0).push(r);
  }
  for (const [clue, items] of byClue.entries()) {
    if (items.length >= 2) {
      return { clue, first: items[0], second: items[1] };
    }
  }
  throw new Error(`No opening group with >= 2 targets found for ${game}`);
}

function findRankedSeqleTargetBySeq(seq) {
  return PUZZLES.seqle.find((r) => r.seq === seq);
}

function attachPageErrorTrap(page) {
  const errors = [];
  page.on("pageerror", (err) => errors.push(err.message));
  return errors;
}

test.describe("BioWordle v0.2.2 E2E Suite", () => {
  test("1. Chemle pre-guess non-leakage, combobox keyboard nav, scaffold threshold (1 before & at unlock), reveal, reload", async ({ page }) => {
    const pageErrors = attachPageErrorTrap(page);
    const utcDate = "2026-11-14";
    const target = dailyTargetForDate("chemle", utcDate);
    // Find a target with cyclic scaffold, or pick a date whose target has_scaffold=true for threshold verification
    let scaffoldDate = utcDate;
    let scaffoldTarget = target;
    for (let d = 1; d <= 30; d++) {
      const ds = `2026-11-${String(d).padStart(2, "0")}`;
      const cand = dailyTargetForDate("chemle", ds);
      if (cand.has_scaffold) {
        scaffoldDate = ds;
        scaffoldTarget = cand;
        break;
      }
    }
    expect(scaffoldTarget.has_scaffold).toBe(true);

    await page.addInitScript((d) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
    }, scaffoldDate);
    await page.goto("/");

    const well = page.locator("#chemle-structure-well");
    await expect(well.locator('[data-placeholder="chemle-preunlock"]')).toBeVisible();
    const level0Text = await page.locator("#chemle-hints-list").innerText();
    expect(level0Text).not.toContain(String(scaffoldTarget.MW));
    expect(level0Text).not.toContain(scaffoldTarget.Category);
    expect(level0Text.toLowerCase()).not.toContain("bits");

    // Pick 6 distinct wrong guesses for full loss/reveal test (first via combobox ArrowDown + Enter)
    const wrongs = PUZZLES.chemle.filter((m) => m.id !== scaffoldTarget.id).slice(0, 6);
    const searchInput = page.locator("#chemle-search-input");
    await searchInput.focus();
    await searchInput.fill(wrongs[0].name.slice(0, 4));
    await expect(page.locator("#chemle-dropdown")).toHaveClass(/open/);
    // Keyboard ArrowDown / Enter navigation with aria-activedescendant
    await searchInput.press("ArrowDown");
    const activeDesc = await searchInput.getAttribute("aria-activedescendant");
    expect(activeDesc).toBeTruthy();
    await searchInput.press("ArrowUp");
    await expect(searchInput).toHaveAttribute("aria-expanded", "true");
    // Choose exact wrongs[0] reliably using filter + Enter
    await searchInput.fill(wrongs[0].name);
    await searchInput.press("Enter");

    // Guesses 1..3 (before Level 1/2 super-class unlock at after_guesses=4): placeholder & Level 0 stay locked
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);
    await expect(well.locator('[data-placeholder="chemle-preunlock"]')).toBeVisible();
    await expect(page.locator("#chemle-candidates-pill")).not.toContainText("36 / 36");

    for (let i = 1; i <= 2; i++) {
      await searchInput.fill(wrongs[i].name);
      await page.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${wrongs[i].id}"]`).click();
    }
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(3);
    await expect(page.locator("#chemle-clue-stage")).toContainText("Level 0");
    await expect(well.locator('[data-placeholder="chemle-preunlock"]')).toBeVisible();

    // Guess #4: Level 1/2 clues unlock (after_guesses=4), 1 guess BEFORE Ranked scaffold unlock (after_guesses=5)
    await searchInput.fill(wrongs[3].name);
    await page.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${wrongs[3].id}"]`).click();
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(4);
    await expect(page.locator("#chemle-clue-stage")).toContainText("Level 2");
    await expect(well.locator('[data-placeholder="chemle-preunlock"]')).toBeVisible();

    // Guess #5 AT scaffold threshold (scaffold_unlock_after_guesses=5): scaffold unlocks!
    await searchInput.fill(wrongs[4].name);
    await page.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${wrongs[4].id}"]`).click();
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(5);
    await expect(well.locator('[data-placeholder="chemle-preunlock"]')).toHaveCount(0);
    await expect(well.locator(".well-overlay-tag")).toContainText(/Scaffold|Skeleton/);

    // Guess #6 -> Loss & Reveal at cap=6
    await searchInput.fill(wrongs[5].name);
    await page.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${wrongs[5].id}"]`).click();
    await expect(page.locator("#chemle-full-props")).toBeVisible();
    await expect(page.locator("#chemle-clue-stage")).toHaveText("Revealed");

    // Reload preserves completed session
    await page.reload();
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(6);
    await expect(page.locator("#chemle-full-props")).toBeVisible();
    expect(pageErrors).toEqual([]);
  });

  test("2. Ranked Seqle typed + legal non-catalog probe, Level-1 threshold (before & at 2), win & reload", async ({ page }) => {
    const pageErrors = attachPageErrorTrap(page);
    const utcDate = "2026-11-14";
    const target = dailyTargetForDate("seqle", utcDate);

    await page.addInitScript((d) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
    }, utcDate);
    await page.goto("/");
    await page.locator("#tab-seqle").click();

    await expect(page.locator("#seqle-structure-well [data-placeholder='seqle-preunlock']")).toBeVisible();

    // Non-catalog 5-AA legal probe 'AALVS' (skip if by chance equals target)
    const probe1 = target.seq === "AALVS" ? "AACSF" : "AALVS";
    const probe2 = target.seq === "FIALS" ? "VVAAS" : "FIALS";
    // Blur focus from button so document keydown receives typing
    await page.evaluate(() => document.activeElement?.blur());
    await page.keyboard.type(probe1);
    await page.keyboard.press("Enter");

    // Guess 1 (1 before Level 1 threshold at after_guesses=2): Level 1 still locked
    await expect(page.locator("#seqle-board .seqle-row").nth(0).locator(".seqle-tile")).toHaveCount(5);
    await expect(page.locator("#seqle-clue-stage")).toContainText("Level 0");
    await expect(page.locator("#seqle-hints-list")).toContainText("Auto-unlocks after Guess 2");

    // Guess 2 (AT threshold): Level 1 unlocks
    await page.keyboard.type(probe2);
    await page.locator("#seqle-btn-submit").click();
    await expect(page.locator("#seqle-clue-stage")).toContainText("Level 1");
    await expect(page.locator("#seqle-hints-list")).toContainText("Clue 2 · Motif Architecture");

    // Solve on Guess 3 via exact target typing
    await page.keyboard.type(target.seq);
    await page.keyboard.press("Enter");
    await expect(page.locator("#seqle-full-props")).toBeVisible();
    await expect(page.locator("#stat-seq")).toHaveText("1");
    await expect(page.locator("#stat-streak")).toHaveText("1");

    await page.reload();
    await page.locator("#tab-seqle").click();
    await expect(page.locator("#seqle-full-props")).toBeVisible();
    await expect(page.locator("#stat-seq")).toHaveText("1");
    await expect(page.locator("#stat-streak")).toHaveText("1");
    expect(pageErrors).toEqual([]);
  });

  test("3. Learn mode deterministic targets, 8-guess cap, category pack filtering, live sim/FGs & hint button", async ({ page }) => {
    const pageErrors = attachPageErrorTrap(page);
    const aspirin = PUZZLES.chemle.find((m) => m.name.toLowerCase().includes("aspirin"));
    const caffeine = PUZZLES.chemle.find((m) => m.name.toLowerCase().includes("caffeine"));
    const npvyk = PUZZLES.seqle.find((m) => m.seq === "NPVYK");
    expect(aspirin && caffeine && npvyk).toBeTruthy();

    await page.addInitScript(
      ({ chemId, seqCode }) => {
        window.__BIOWORDLE_TEST_LEARN_TARGETS__ = { chemle: chemId, seqle: seqCode };
      },
      { chemId: aspirin.id, seqCode: npvyk.seq }
    );
    await page.goto("/?profile=learn");

    // Category pack dropdown enabled in Learn; switching category changes denominator in candidates pill
    const catSelect = page.locator("#chemle-cat-filter");
    await expect(catSelect).toBeEnabled();
    await expect(page.locator("#chemle-candidates-pill")).toContainText("/ 36");

    // Unlock Level 1 clue manually BEFORE guessing via Learn Hint button
    await expect(page.locator("#chemle-clue-stage")).toContainText("Level 0");
    await page.locator("#btn-hint").click();
    await expect(page.locator("#chemle-clue-stage")).toContainText("Level 1");

    // Probe a different molecule and inspect Learn live rounded Tanimoto + FG overlap
    const searchInput = page.locator("#chemle-search-input");
    await searchInput.fill(caffeine.name);
    await page.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${caffeine.id}"]`).click();
    const cardText = await page.locator("#chemle-guesses .chemle-guess-card").first().innerText();
    expect(cardText).toMatch(/\d+%/);
    await expect(page.locator("#chemle-guesses .fg-chips")).toBeVisible();

    // Switch category pack -> denominator shrinks while target moves to active category & Level-0 compatibility applies
    const targetCat = aspirin.Category;
    const packMols = PUZZLES.chemle.filter((m) => m.Category === targetCat);
    const catCount = packMols.length;
    const aspirinLvl0 = aspirin.clue_ladder.learn.find((c) => c.level === 0).text;
    const compatLvl0 = packMols.filter((m) => m.clue_ladder.learn.find((c) => c.level === 0).text === aspirinLvl0).length;
    await catSelect.selectOption(targetCat);
    await expect(page.locator("#chemle-candidates-pill")).toContainText(`Candidates: ${compatLvl0} / ${catCount}`);

    // Switch to Seqle in Learn -> 8 rows & deterministic NPVYK target reveals Consensus exemplar
    await page.locator("#tab-seqle").click();
    await expect(page.locator("#seqle-board .seqle-row")).toHaveCount(8);
    expect(pageErrors).toEqual([]);
  });

  test("4. Phase 2 F3 target-independent pre-unlock wells & F4 provenance DOM assertions", async ({ page }) => {
    const pageErrors = attachPageErrorTrap(page);
    const seqPair = findOpeningGroupPair("seqle");
    const chemPair = findOpeningGroupPair("chemle");

    // Pre-unlock Chemle well byte equality across two distinct targets in the same Level-0 group
    await page.addInitScript(({ id1, sq1 }) => {
      window.__BIOWORDLE_TEST_LEARN_TARGETS__ = { chemle: id1, seqle: sq1 };
    }, { id1: chemPair.first.id, sq1: seqPair.first.seq });
    await page.goto("/?profile=learn");
    const chemWell1 = await page.locator("#chemle-structure-well").innerHTML();
    await page.locator("#tab-seqle").click();
    const seqWell1 = await page.locator("#seqle-structure-well").innerHTML();

    await page.addInitScript(({ id2, sq2 }) => {
      window.__BIOWORDLE_TEST_LEARN_TARGETS__ = { chemle: id2, seqle: sq2 };
    }, { id2: chemPair.second.id, sq2: seqPair.second.seq });
    await page.goto("/?profile=learn");
    const chemWell2 = await page.locator("#chemle-structure-well").innerHTML();
    await page.locator("#tab-seqle").click();
    const seqWell2 = await page.locator("#seqle-structure-well").innerHTML();

    expect(chemWell1).toBe(chemWell2);
    expect(seqWell1).toBe(seqWell2);

    // Picker omits retired KKMPQ and includes KEKSD
    const pickerOptions = await page.locator("#seqle-motif-picker").innerText();
    expect(pickerOptions).not.toContain("KKMPQ");
    expect(pickerOptions).toContain("KEKSD");

    // Post-reveal provenance badges & wording checks across NPVYK, INFEK, PPPVP, DSGLD
    const checkPostReveal = async (seqCode, expectedRegexList, forbiddenList = []) => {
      await page.addInitScript((s) => {
        window.__BIOWORDLE_TEST_LEARN_TARGETS__ = { seqle: s };
      }, seqCode);
      await page.goto("/?profile=learn");
      await page.locator("#tab-seqle").click();
      await page.locator("#btn-giveup").click();
      const propsText = await page.locator("#seqle-full-props").innerText();
      for (const rx of expectedRegexList) {
        expect(propsText).toMatch(rx);
      }
      for (const bad of forbiddenList) {
        expect(propsText).not.toContain(bad);
      }
    };

    await checkPostReveal("NPVYK", [/Consensus exemplar/, /Teaching Sequence/, /Functional Rule/], ["Parent Sequence Window"]);
    await checkPostReveal("INFEK", [/Verified UniProt Instance/, /Parent Sequence Window/, /260/]);
    await checkPostReveal("PPPVP", [/Context-only evidence/, /1GBR/]);
    await checkPostReveal("DSGLD", [/9T9W/], ["1P22"]);
    expect(pageErrors).toEqual([]);
  });

  test("5. Phase 1 America/New_York page.clock UTC rollover (timer-driven, stale-click, and repeated-target date)", async ({ browser }) => {
    const context = await browser.newContext({ timezoneId: "America/New_York" });
    const page = await context.newPage();
    const pageErrors = attachPageErrorTrap(page);

    // (a) Timer-driven rollover: 23:59:50Z in New York (18:59:50 local) -> ticks past 00:00:00Z
    await page.clock.install({ time: new Date("2026-11-14T23:59:50Z") });
    await page.goto("/");
    const targetDay1 = dailyTargetForDate("chemle", "2026-11-14");
    const targetDay2 = dailyTargetForDate("chemle", "2026-11-15");
    await page.clock.runFor(15_000);
    await expect(page.locator("#utc-rollover-banner")).toHaveClass(/show/);
    await expect(page.locator("#utc-rollover-banner")).toContainText("2026-11-15");
    await context.close();

    // (b) Stale-click rollover with timer paused (user leaves tab idle overnight and clicks without tick)
    const ctx2 = await browser.newContext({ timezoneId: "America/New_York" });
    const p2 = await ctx2.newPage();
    const errs2 = attachPageErrorTrap(p2);
    await p2.clock.install({ time: new Date("2026-11-14T23:55:00Z") });
    await p2.goto("/");
    // Advance system Date without firing setInterval countdown
    await p2.clock.setFixedTime(new Date("2026-11-15T00:05:00Z"));
    const staleGuess = PUZZLES.chemle.find((m) => m.id !== targetDay1.id && m.id !== targetDay2.id);
    await p2.locator("#chemle-search-input").fill(staleGuess.name);
    await p2.locator("#chemle-search-input").press("Enter");
    // Stale click must trigger rollover notice and keep guess count at 0 for 2026-11-15
    await expect(p2.locator("#utc-rollover-banner")).toHaveClass(/show/);
    await expect(p2.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(0);
    expect(errs2).toEqual([]);
    await ctx2.close();

    // (c) Adjacent UTC date where Ranked target repeats must still reset guesses
    const repeatPair = findRepeatedDailyTargetDates("chemle");
    const ctx3 = await browser.newContext({ timezoneId: "America/New_York" });
    const p3 = await ctx3.newPage();
    await p3.clock.install({ time: new Date(`${repeatPair.dateA}T23:59:50Z`) });
    await p3.goto("/");
    const probe = PUZZLES.chemle.find((m) => m.id !== repeatPair.target.id);
    await p3.locator("#chemle-search-input").fill(probe.name);
    await p3.locator("#chemle-search-input").press("Enter");
    await expect(p3.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);
    await p3.clock.runFor(15_000);
    await expect(p3.locator("#utc-rollover-banner")).toContainText(repeatPair.dateB);
    await expect(p3.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(0);
    expect(pageErrors).toEqual([]);
    await ctx3.close();
  });

  test("6. Phase 4 F2 cross-tab sequential & concurrent wins, session lock, pending-win recovery, streak & Learn isolation", async ({ browser }) => {
    const context = await browser.newContext();
    const utcDate = "2026-11-14";
    const nextDay = "2026-11-15";
    const skipDay = "2026-11-16";
    const chemTarget = dailyTargetForDate("chemle", utcDate);
    const seqTarget = dailyTargetForDate("seqle", utcDate);

    await context.addInitScript((d) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
    }, utcDate);

    // Part A: Simultaneous Chemle win in Tab 1 & Seqle win in Tab 2 via Promise.all
    const tab1 = await context.newPage();
    const tab2 = await context.newPage();
    await Promise.all([tab1.goto("/"), tab2.goto("/")]);
    await tab2.locator("#tab-seqle").click();

    await tab1.locator("#chemle-search-input").fill(chemTarget.name);
    await tab2.evaluate(() => document.activeElement?.blur());
    await tab2.keyboard.type(seqTarget.seq);

    await Promise.all([
      tab1.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${chemTarget.id}"]`).click(),
      tab2.keyboard.press("Enter"),
    ]);

    await expect(tab1.locator("#chemle-full-props")).toBeVisible();
    await expect(tab2.locator("#seqle-full-props")).toBeVisible();
    await expect(tab1.locator("#stat-chem")).toHaveText("1");
    await expect(tab1.locator("#stat-seq")).toHaveText("1");
    await expect(tab1.locator("#stat-streak")).toHaveText("1");

    // Verify Learn & seeded mode wins increment win counts without touching streak
    const seedChemTarget = dailyTargetForDate("chemle", utcDate, "CLASS-TEST");
    await tab1.goto("/?seed=CLASS-TEST");
    await tab1.locator("#chemle-search-input").fill(seedChemTarget.name);
    await tab1.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${seedChemTarget.id}"]`).click();
    await expect(tab1.locator("#stat-chem")).toHaveText("2");
    await expect(tab1.locator("#stat-streak")).toHaveText("1");

    await context.close();

    // Part B: Competing-move session lock + interrupted pending-win recovery + missed-day streak reset
    const ctxB = await browser.newContext();
    await ctxB.addInitScript((d) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
    }, utcDate);
    const pA = await ctxB.newPage();
    // Suppress storage listener on pB to simulate stale offline tab that didn't see pA's move yet
    const pB = await ctxB.newPage();
    await pB.addInitScript(() => {
      const origAdd = window.addEventListener.bind(window);
      window.addEventListener = (type, fn, opts) => {
        if (type === "storage") return;
        return origAdd(type, fn, opts);
      };
    });
    await Promise.all([pA.goto("/"), pB.goto("/")]);
    const wA = PUZZLES.chemle.find((m) => m.id !== chemTarget.id);
    const wB = PUZZLES.chemle.find((m) => m.id !== chemTarget.id && m.id !== wA.id);
    await pA.locator("#chemle-search-input").fill(wA.name);
    await pA.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${wA.id}"]`).click();
    await expect(pA.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);

    // Stale pB tries to submit wB -> rejected cleanly, adopts pA's wA guess instead of overwriting
    await pB.locator("#chemle-search-input").fill(wB.name);
    await pB.locator(`#chemle-dropdown .chemle-opt[data-mol-id="${wB.id}"]`).click();
    await expect(pB.locator("#utc-rollover-banner")).toContainText(/Concurrent Tab Sync/);
    await expect(pB.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);
    await expect(pB.locator("#chemle-guesses .chemle-guess-card").first()).toContainText(wA.name);
    await ctxB.close();

    // Part C: Pending-win boot recovery + yesterday streak=1 on nextDay + missed-day reset to 0 on skipDay
    const ctxC = await browser.newContext();
    const pC = await ctxC.newPage();
    await pC.addInitScript((d) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
    }, utcDate);
    await pC.goto("/");
    // Simulate crashed tab: wrote pending win intent before finishing stats flush
    await pC.evaluate((d) => {
      localStorage.setItem("biowordle_stats_v2", JSON.stringify({ schema: 2, chemWins: 0, seqWins: 0, streak: 0, maxStreak: 0 }));
      localStorage.setItem(
        "biowordle_pending_wins_v2",
        JSON.stringify([{ creditId: `daily:chemle:${d}`, game: "chemle", isDaily: true, utcDate: d }])
      );
    }, utcDate);
    await pC.reload();
    await expect(pC.locator("#stat-chem")).toHaveText("1");
    await expect(pC.locator("#stat-streak")).toHaveText("1");

    // Next day before playing: streak preserved at 1
    const pC2 = await ctxC.newPage();
    await pC2.addInitScript((d) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
    }, nextDay);
    await pC2.goto("/");
    await expect(pC2.locator("#stat-streak")).toHaveText("1");

    // Skip nextDay & open on skipDay: streak displays 0
    const pC3 = await ctxC.newPage();
    await pC3.addInitScript((d) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
    }, skipDay);
    await pC3.goto("/");
    await expect(pC3.locator("#stat-streak")).toHaveText("0");
    await expect(pC3.locator("#stat-chem")).toHaveText("1");
    await ctxC.close();
  });

  test("7. Phase 4 F6 legacy key compatible import (919ad8afb190) and incompatible target archiving (KKMPQ / mismatch)", async ({ page }) => {
    const pageErrors = attachPageErrorTrap(page);
    const utcDate = "2026-11-14";
    const chemTarget = dailyTargetForDate("chemle", utcDate);
    const wrongChem = PUZZLES.chemle.find((m) => m.id !== chemTarget.id);
    const seqTarget = dailyTargetForDate("seqle", utcDate);

    // 7a. Compatible legacy Chemle record under 919ad8afb190
    await page.addInitScript(({ d, chemId, guessSeqId, seqCode }) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
      localStorage.setItem(
        `biowordle_daily_v2:chemle:${d}:919ad8afb190`,
        JSON.stringify({
          targetId: chemId,
          guesses: [guessSeqId, chemId],
          won: true,
          gameOver: true,
          hintLevelUnlocked: 0,
          utcDate: d,
        })
      );
      // Incompatible legacy Seqle record citing retired KKMPQ
      localStorage.setItem(
        `biowordle_daily_v2:seqle:${d}:919ad8afb190`,
        JSON.stringify({
          targetId: "KKMPQ",
          guesses: ["KKMPQ"],
          won: true,
          gameOver: true,
          hintLevelUnlocked: 0,
          utcDate: d,
        })
      );
    }, { d: utcDate, chemId: chemTarget.id, guessSeqId: wrongChem.id, seqCode: seqTarget.seq });

    await page.goto("/");
    // Compatible Chemle imported cleanly with 2 guesses & won state (no archive banner)
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(2);
    await expect(page.locator('[data-testid="chemle-legacy-archive"]')).toHaveCount(0);

    // Incompatible KKMPQ Seqle archived cleanly & board fresh at 0 guesses for real seqTarget
    await page.locator("#tab-seqle").click();
    await expect(page.locator('[data-testid="seqle-legacy-archive"]')).toBeVisible();
    await expect(page.locator('[data-testid="seqle-legacy-archive"]')).toContainText("KKMPQ");
    await expect(page.locator("#seqle-board .seqle-row.filled, #seqle-board .seqle-tile.state-exact")).toHaveCount(0);

    // 7b. Mismatched Chemle legacy target under 919ad8afb190 archives cleanly instead of re-scoring
    const ctx2 = await page.context().browser().newContext();
    const page2 = await ctx2.newPage();
    await page2.addInitScript(({ d, badId }) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = d;
      localStorage.setItem(
        `biowordle_daily_v2:chemle:${d}:919ad8afb190`,
        JSON.stringify({
          targetId: badId,
          guesses: [badId],
          won: true,
          gameOver: true,
          hintLevelUnlocked: 0,
          utcDate: d,
        })
      );
    }, { d: utcDate, badId: wrongChem.id });
    await page2.goto("/");
    await expect(page2.locator('[data-testid="chemle-legacy-archive"]')).toBeVisible();
    await expect(page2.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(0);
    await ctx2.close();
    expect(pageErrors).toEqual([]);
  });

  test("8. Phase 1 offline file:// standalone HTML smoke test (zero external requests, interactive Chemle + Seqle)", async ({ page }) => {
    const pageErrors = attachPageErrorTrap(page);
    const externalUrls = [];
    await page.route("**/*", async (route) => {
      const url = route.request().url();
      if (!url.startsWith("file://") && !url.startsWith("data:") && !url.startsWith("blob:")) {
        externalUrls.push(url);
        await route.abort();
        return;
      }
      await route.continue();
    });

    const standaloneUrl = pathToFileURL(path.join(ROOT, "dist", "biowordle_standalone.html")).href;
    await page.goto(standaloneUrl);

    await expect(page.locator("#chemle-structure-well")).toBeVisible();
    await expect(page.locator("#app-meta-subtitle")).toContainText("v0.2.2");

    // Submit 1 Chemle guess in offline file:// mode (verifies lock fallback works without navigator.locks/HTTP origin errors)
    const mol = PUZZLES.chemle[0];
    await page.locator("#chemle-search-input").fill(mol.name);
    await page.locator("#chemle-search-input").press("Enter");
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);

    // Switch to Seqle & submit 1 AA guess
    await page.locator("#tab-seqle").click();
    await page.evaluate(() => document.activeElement?.blur());
    await page.keyboard.type("AALVS");
    await page.keyboard.press("Enter");
    await expect(page.locator("#seqle-board .seqle-row").nth(0).locator(".seqle-tile")).toHaveCount(5);

    expect(externalUrls).toEqual([]);
    expect(pageErrors).toEqual([]);
  });

  test("9. Chemle pool separation, Ranked rejection of Learn-only probes, Learn pack denominators & invalid pack fallback", async ({ page }) => {
    const pageErrors = attachPageErrorTrap(page);
    const clonedData = JSON.parse(JSON.stringify(PUZZLES));
    const baseMol = clonedData.chemle[0];
    const learnMolId = "synthetic_learn_probe";
    const learnMol = {
      ...JSON.parse(JSON.stringify(baseMol)),
      id: learnMolId,
      name: "Synthetic Learn Probe",
      Category: baseMol.Category,
      inchikey: "SYNTHLEARNPROB-UHFFFAOYSA-N",
      game_eligibility: ["learn"],
      guess_eligibility: ["learn"],
    };
    clonedData.chemle.push(learnMol);
    clonedData.chemle_similarity[learnMolId] = {};
    for (const m of clonedData.chemle) {
      const s = m.id === learnMolId ? 100 : (clonedData.chemle_similarity[baseMol.id][m.id] ?? 0);
      clonedData.chemle_similarity[learnMolId][m.id] = s;
      clonedData.chemle_similarity[m.id][learnMolId] = s;
    }
    clonedData.metadata.chemle_count = clonedData.chemle.length;
    clonedData.metadata.chemle_ranked_count = 36;
    clonedData.metadata.chemle_ranked_guess_count = 36;
    clonedData.metadata.chemle_learn_count = 37;
    clonedData.metadata.chemle_learn_guess_count = 37;

    await page.route("**/data/puzzles_data.js", async (route) => {
      await route.fulfill({
        status: 200,
        contentType: "application/javascript; charset=utf-8",
        body: `window.BIOWORDLE_DATA = ${JSON.stringify(clonedData)};`,
      });
    });

    await page.addInitScript(() => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = "2026-11-14";
    });
    await page.goto("/");

    // Ranked mode sees 36 Ranked / 36 denominator and excludes the Learn-only probe
    await expect(page.locator("#badge-chemle-count")).toHaveText("36 Ranked");
    await expect(page.locator("#chemle-candidates-pill")).toHaveText(/^Candidates: \d+ \/ 36$/);
    await page.locator("#chemle-search-input").fill("Synthetic Learn Probe");
    await expect(page.locator("#chemle-dropdown")).toContainText("No remaining molecules match");
    await page.locator("#chemle-search-input").press("Enter");
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(0);

    // Ranked keyboard entry works cleanly
    await page.locator("#chemle-search-input").fill("Aspirin");
    await page.locator("#chemle-search-input").press("ArrowDown");
    await page.locator("#chemle-search-input").press("Enter");
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);

    // Switch to Learn mode -> sees 37 Mols / 37 denominator and accepts Synthetic Learn Probe via mouse click
    await page.locator("#btn-playstyle").click();
    await expect(page.locator("#badge-chemle-count")).toHaveText("37 Mols");
    await expect(page.locator("#chemle-candidates-pill")).toHaveText(/^Candidates: \d+ \/ 37$/);
    await expect(page.locator("#chemle-cat-filter")).toContainText("All Answer Packs (37 Mols)");

    await page.locator("#chemle-search-input").fill("Synthetic Learn Probe");
    const opt = page.locator(`#chemle-opt-${learnMolId}`);
    await expect(opt).toBeVisible();
    await opt.click();
    await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);
    await expect(page.locator("#chemle-guesses")).toContainText("Synthetic Learn Probe");

    // Pack filter shows exact Learn answer denominator and invalid filter selection falls back safely
    const targetCat = baseMol.Category;
    const expectedPackAnsCount = clonedData.chemle.filter(
      (m) => (m.game_eligibility || []).includes("learn") && m.Category === targetCat
    ).length;
    await page.locator("#chemle-cat-filter").selectOption(targetCat);
    await expect(page.locator("#chemle-cat-filter option:checked")).toContainText(`(${expectedPackAnsCount})`);
    await expect(page.locator("#chemle-candidates-pill")).toHaveText(
      new RegExp(`^Candidates: \\d+ / ${expectedPackAnsCount}$`)
    );

    // Selecting an invalid pack via script + change handler falls back cleanly to ALL
    await page.evaluate(() => {
      const sel = document.getElementById("chemle-cat-filter");
      const badOpt = document.createElement("option");
      badOpt.value = "NONEXISTENT_PACK";
      badOpt.textContent = "Broken Pack";
      sel.appendChild(badOpt);
      sel.value = "NONEXISTENT_PACK";
      sel.dispatchEvent(new Event("change", { bubbles: true }));
    });
    await expect(page.locator("#chemle-cat-filter")).toHaveValue("ALL");

    // Target selection near start and end of sorted pool both resolve & accept guesses
    const pool = BioWordleEngine.getAnswerPool(PUZZLES, "chemle", "ranked", "ALL");
    let startSeed = null;
    let endSeed = null;
    for (let i = 0; i < 400 && (!startSeed || !endSeed); i++) {
      const candidateSeed = `boundary-check-${i}`;
      const idx = BioWordleEngine.selectDailyIndex(
        "chemle",
        "ranked",
        "2026-11-14",
        CHEM_SCHEDULE_HASH,
        pool.length,
        candidateSeed
      );
      if (!startSeed && idx <= 2) startSeed = candidateSeed;
      if (!endSeed && idx >= pool.length - 3) endSeed = candidateSeed;
    }
    expect(startSeed && endSeed).toBeTruthy();

    for (const seed of [startSeed, endSeed]) {
      await page.goto(`/?seed=${seed}`);
      await expect(page.locator("#chemle-structure-well")).toBeVisible();
      await page.locator("#chemle-search-input").fill("Caffeine");
      await page.locator("#chemle-search-input").press("Enter");
      await expect(page.locator("#chemle-guesses .chemle-guess-card")).toHaveCount(1);
    }

    expect(pageErrors).toEqual([]);
  });
});
