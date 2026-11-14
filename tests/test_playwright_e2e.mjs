import fs from "node:fs";
import path from "node:path";
import assert from "node:assert/strict";
import { createRequire } from "node:module";
import { fileURLToPath, pathToFileURL } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..");
const screenshotsDir = path.join(__dirname, "screenshots");

if (!process.env.LOCALAPPDATA) {
  throw new Error("LOCALAPPDATA environment variable is required to locate Antigravity playwright-core.");
}

const require = createRequire(import.meta.url);
const playwrightCorePath = path.join(
  process.env.LOCALAPPDATA,
  "Programs",
  "Antigravity IDE",
  "resources",
  "app",
  "node_modules",
  "playwright-core"
);
const { chromium } = require(playwrightCorePath);

const indexUrl = pathToFileURL(path.join(repoRoot, "index.html")).href;
const standaloneUrl = pathToFileURL(path.join(repoRoot, "dist", "biowordle_standalone.html")).href;

async function submitChemleMoleculeByUi(page, molId, searchTerm) {
  const input = page.locator("#chemle-search-input");
  await input.click();
  await input.fill(searchTerm);
  const opt = page.locator(`#chemle-opt-${molId}`);
  await opt.waitFor({ state: "visible", timeout: 5000 });
  await opt.click();
}

async function getRankedDailyTarget(page, game, utcDateOverride) {
  return page.evaluate(
    ({ gameArg, dateArg }) => {
      const data = window.BIOWORDLE_DATA;
      const engine = window.BioWordleEngine;
      const utc = dateArg || engine.getUtcDateString();
      const pool = engine.getAnswerPool(data, gameArg, "ranked", "ALL");
      const hash = (data.metadata && data.metadata.schedule_version_hash) || data.metadata.dataset_version_hash;
      const idx = engine.selectDailyIndex(gameArg, "ranked", utc, hash, pool.length, null);
      return { target: pool[idx], pool };
    },
    { gameArg: game, dateArg: utcDateOverride || null }
  );
}

async function createCleanPage(browser, { utcDate = "2026-11-14", v1Stats = null } = {}) {
  const context = await browser.newContext({ viewport: { width: 1400, height: 960 } });
  await context.addInitScript(
    ({ initialUtc, initialV1 }) => {
      window.__BIOWORDLE_TEST_UTC_DATE__ = initialUtc;
      if (!sessionStorage.getItem("__bw_e2e_seeded")) {
        sessionStorage.setItem("__bw_e2e_seeded", "1");
        localStorage.clear();
        if (initialV1) {
          localStorage.setItem("biowordle_stats_v1", JSON.stringify(initialV1));
        }
      }
    },
    { initialUtc: utcDate, initialV1: v1Stats }
  );
  const page = await context.newPage();
  const pageErrors = [];
  page.on("pageerror", (err) => pageErrors.push(err.message));
  await page.goto(indexUrl);
  await page.waitForLoadState("domcontentloaded");
  return { context, page, pageErrors };
}

async function run() {
  fs.mkdirSync(screenshotsDir, { recursive: true });
  const browser = await chromium.launch({
    channel: "msedge",
    headless: true,
    args: ["--no-sandbox", "--disable-gpu"],
  });

  try {
    // -----------------------------------------------------------------------
    // E2E Test 1: Ranked Manual Hint Lockout (Repeated Clicks Before Guess 1)
    // -----------------------------------------------------------------------
    {
      const { context, page, pageErrors } = await createCleanPage(browser);
      try {
        const hintBtn = page.locator("#btn-hint");
        assert.equal(await hintBtn.isDisabled(), true, "Ranked #btn-hint must be disabled");
        assert.ok(
          (await hintBtn.textContent()).includes("Auto-Clues (Ranked)"),
          "Ranked #btn-hint label must clearly state Auto-Clues"
        );
        assert.equal(
          (await page.locator("#chemle-clue-stage").textContent()).trim(),
          "Level 0 · k-Anonymous"
        );

        for (let i = 0; i < 5; i++) {
          await hintBtn.click({ force: true });
          await page.evaluate(() => document.getElementById("btn-hint").click());
        }
        await page.waitForTimeout(150);

        assert.equal(
          (await page.locator("#chemle-clue-stage").textContent()).trim(),
          "Level 0 · k-Anonymous",
          "5 Hint clicks in Ranked must NOT unlock clues before Guess 1"
        );
        const chemHintsText = await page.locator("#chemle-hints-list").textContent();
        assert.ok(!chemHintsText.includes("Broad Super-Class:"), "Level 1 clue must remain locked");
        assert.ok(!chemHintsText.includes("or click 💡 Hint"), "Ranked locked copy must not invite Hint clicks");

        // Same Ranked Hint lockout check on Seqle tab
        await page.locator("#tab-seqle").click();
        assert.equal(await hintBtn.isDisabled(), true, "Ranked Seqle #btn-hint must stay disabled");
        for (let i = 0; i < 5; i++) {
          await hintBtn.click({ force: true });
        }
        assert.equal(
          (await page.locator("#seqle-clue-stage").textContent()).trim(),
          "Level 0 · k-Anonymous",
          "5 Hint clicks in Ranked Seqle must NOT unlock clues"
        );

        // Switch to Learn mode -> Hint button enables and unlocks progressive clues
        await page.locator("#btn-playstyle").click();
        assert.equal(await hintBtn.isDisabled(), false, "Learn #btn-hint must be enabled");
        await hintBtn.click();
        assert.ok(
          (await page.locator("#seqle-clue-stage").textContent()).startsWith("Level 1"),
          "Clicking Hint in Learn Seqle must unlock Level 1"
        );
        await hintBtn.click();
        assert.ok((await page.locator("#seqle-clue-stage").textContent()).startsWith("Level 2"));
        await hintBtn.click();
        assert.ok((await page.locator("#seqle-clue-stage").textContent()).startsWith("Level 3"));
        assert.equal(await hintBtn.isDisabled(), true, "Learn #btn-hint disables after Level 3");

        await page.screenshot({ path: path.join(screenshotsDir, "e2e_1_hint_lockout.png"), fullPage: true });
        assert.deepEqual(pageErrors, []);
        console.log("[E2E 1 PASS] Ranked manual Hint lockout & Learn Hint progression verified.");
      } finally {
        await context.close();
      }
    }

    // -----------------------------------------------------------------------
    // E2E Test 2: Learn Mode Candidate Filtering Completeness
    // -----------------------------------------------------------------------
    {
      const { context, page, pageErrors } = await createCleanPage(browser);
      try {
        await page.locator("#btn-playstyle").click(); // Switch to Learn Mode (Chemle)
        const pillText0 = (await page.locator("#chemle-candidates-pill").textContent()).trim();
        const m0 = pillText0.match(/Candidates:\s*(\d+)\s*\/\s*(\d+)/);
        assert.ok(m0, `Unexpected candidates pill text: ${pillText0}`);
        const c0 = Number(m0[1]);
        assert.ok(c0 < 36 && c0 >= 4, `Learn Level-0 band must narrow initial candidates (${c0}/36)`);

        // Clicking Hint in Learn unlocks Level 1 and updates Candidates pill immediately
        await page.locator("#btn-hint").click();
        const pillText1 = (await page.locator("#chemle-candidates-pill").textContent()).trim();
        const c1 = Number(pillText1.match(/Candidates:\s*(\d+)/)[1]);
        assert.ok(c1 <= c0 && c1 >= 1, `Learn Level 1 hint must preserve/narrow candidates (${c1} <= ${c0})`);

        // Submit a guaranteed non-target guess in Learn mode and check SIM_PCT % + narrowed candidates
        await submitChemleMoleculeByUi(page, "caffeine", "Caffeine");
        if ((await page.locator("#chemle-banner").textContent()).includes("Solved")) {
          await page.locator("#btn-new-puzzle").click();
          await page.locator("#btn-hint").click();
          await submitChemleMoleculeByUi(page, "caffeine", "Caffeine");
        }
        const guessesHtml = await page.locator("#chemle-guesses").textContent();
        assert.ok(guessesHtml.includes("%"), "Learn mode guess card must render exact Tanimoto percentage");

        const pillText2 = (await page.locator("#chemle-candidates-pill").textContent()).trim();
        const c2 = Number(pillText2.match(/Candidates:\s*(\d+)/)[1]);
        assert.ok(c2 >= 1 && c2 <= c1, `Post-guess Learn candidates must be 1..${c1}, got ${c2}`);

        // Also check Learn Seqle candidate filtering
        await page.locator("#tab-seqle").click();
        await page.locator("#seqle-motif-picker").selectOption("YIGSR");
        await page.locator("#seqle-btn-submit").click();
        if ((await page.locator("#seqle-banner").textContent()).includes("Motif Decoded")) {
          await page.locator("#btn-new-puzzle").click();
          await page.locator("#seqle-motif-picker").selectOption("YIGSR");
          await page.locator("#seqle-btn-submit").click();
        }
        const seqPill = (await page.locator("#seqle-candidates-pill").textContent()).trim();
        const sc = Number(seqPill.match(/Candidates:\s*(\d+)/)[1]);
        assert.ok(sc >= 1 && sc < 30, `Learn Seqle guess must narrow candidates <30, got ${sc}`);
        await page.locator("#tab-chemle").click();

        await page.screenshot({ path: path.join(screenshotsDir, "e2e_2_learn_filtering.png"), fullPage: true });
        assert.deepEqual(pageErrors, []);
        console.log(`[E2E 2 PASS] Learn candidate filtering narrowed ${c0} -> ${c1} -> ${c2} (Seqle -> ${sc}).`);
      } finally {
        await context.close();
      }
    }

    // -----------------------------------------------------------------------
    // E2E Test 3: Puzzle Completion, Persistence Across Reload & v1->v2 Migration
    // -----------------------------------------------------------------------
    {
      const { context, page, pageErrors } = await createCleanPage(browser, {
        utcDate: "2026-11-14",
        v1Stats: { chemWins: 5, seqWins: 2, streak: 9, maxStreak: 12 },
      });
      try {
        // Verify v1->v2 migration preserved win totals and reset unverified streak
        assert.equal((await page.locator("#stat-chem").textContent()).trim(), "5");
        assert.equal((await page.locator("#stat-seq").textContent()).trim(), "2");
        assert.equal((await page.locator("#stat-streak").textContent()).trim(), "0");

        const statsV2BeforeWin = await page.evaluate(() => JSON.parse(localStorage.getItem("biowordle_stats_v2")));
        assert.equal(statsV2BeforeWin.schema, 2);
        assert.equal(statsV2BeforeWin.chemWins, 5);
        assert.equal(statsV2BeforeWin.seqWins, 2);
        assert.equal(statsV2BeforeWin.streak, 0);
        assert.equal(statsV2BeforeWin.maxStreak, 0);

        const { target, pool } = await getRankedDailyTarget(page, "chemle", "2026-11-14");
        const decoy = pool.find((m) => m.id !== target.id);
        await submitChemleMoleculeByUi(page, decoy.id, decoy.name);
        await submitChemleMoleculeByUi(page, target.id, target.name);

        const winBanner = (await page.locator("#chemle-banner").textContent()).trim();
        assert.ok(winBanner.includes("Solved in 2/6!"), `Expected win banner, got: ${winBanner}`);
        assert.equal((await page.locator("#stat-chem").textContent()).trim(), "6");
        assert.equal((await page.locator("#stat-streak").textContent()).trim(), "1");

        // Reload page and confirm persisted daily solve & streak
        await page.reload({ waitUntil: "domcontentloaded" });
        assert.equal((await page.locator("#stat-chem").textContent()).trim(), "6");
        assert.equal((await page.locator("#stat-streak").textContent()).trim(), "1");
        assert.ok((await page.locator("#chemle-banner").textContent()).includes("Solved in 2/6!"));
        assert.equal(await page.locator("#chemle-guesses .chemle-guess-card").count(), 2);

        await page.screenshot({ path: path.join(screenshotsDir, "e2e_3_completion_persistence.png"), fullPage: true });
        assert.deepEqual(pageErrors, []);
        console.log("[E2E 3 PASS] v1->v2 stats migration, Ranked solve, and reload persistence verified.");
      } finally {
        await context.close();
      }
    }

    // -----------------------------------------------------------------------
    // E2E Test 4: Answer-Pool Integrity Under Learn Category Pack Filtering
    // -----------------------------------------------------------------------
    {
      const { context, page, pageErrors } = await createCleanPage(browser);
      try {
        await page.locator("#btn-playstyle").click(); // Learn Mode
        await page.locator("#chemle-cat-filter").selectOption("Therapeutic Drug");

        // Even with Therapeutic Drug pack active, non-pack molecule Adrenaline must remain searchable in guessPool
        const input = page.locator("#chemle-search-input");
        await input.click();
        await input.fill("Adrenaline");
        await page.locator("#chemle-opt-adrenaline").waitFor({ state: "visible", timeout: 5000 });
        await page.locator("#chemle-opt-adrenaline").click();
        assert.equal(await page.locator("#chemle-guesses .chemle-guess-card").count(), 1);

        // Reveal & verify Learn pack target belongs to Therapeutic Drug AND is in guessPool
        await page.locator("#btn-giveup").click();
        const bannerText = await page.locator("#chemle-banner").textContent();
        assert.ok(bannerText.includes("Therapeutic Drug"), `Expected Therapeutic Drug target, got: ${bannerText}`);

        await page.screenshot({ path: path.join(screenshotsDir, "e2e_4_pack_guesspool.png"), fullPage: true });
        assert.deepEqual(pageErrors, []);
        console.log("[E2E 4 PASS] Answer-pool & guessPool integrity under Learn pack filter verified.");
      } finally {
        await context.close();
      }
    }

    // -----------------------------------------------------------------------
    // E2E Test 5: UTC Date Rollover Detection & Stale-Session Prevention
    // -----------------------------------------------------------------------
    {
      const { context, page, pageErrors } = await createCleanPage(browser, { utcDate: "2026-11-14" });
      try {
        const day1 = await getRankedDailyTarget(page, "chemle", "2026-11-14");
        const day2 = await getRankedDailyTarget(page, "chemle", "2026-11-15");
        assert.notEqual(day1.target.id, day2.target.id, "2026-11-14 and 2026-11-15 should have distinct test targets");

        const day1Decoy = day1.pool.find((m) => m.id !== day1.target.id);
        await submitChemleMoleculeByUi(page, day1Decoy.id, day1Decoy.name);
        assert.equal((await page.locator("#chemle-try-counter").textContent()).trim(), "Guesses: 1 / 6");

        // Simulate UTC midnight rollover BEFORE submitting day1.target
        await page.evaluate(() => {
          window.__BIOWORDLE_TEST_UTC_DATE__ = "2026-11-15";
        });
        await submitChemleMoleculeByUi(page, day1.target.id, day1.target.name);

        // Stale win MUST be blocked; rollover notice shown; puzzle reset to 2026-11-15
        const rolloverVisible = await page.locator("#utc-rollover-banner.show").isVisible();
        assert.equal(rolloverVisible, true, "UTC rollover banner must become visible");
        const rolloverText = await page.locator("#utc-rollover-banner").textContent();
        assert.ok(rolloverText.includes("2026-11-14") && rolloverText.includes("2026-11-15"));
        assert.equal((await page.locator("#chemle-try-counter").textContent()).trim(), "Guesses: 0 / 6");
        assert.equal((await page.locator("#stat-chem").textContent()).trim(), "0", "Stale day win must NOT increment stats");
        assert.equal((await page.locator("#stat-streak").textContent()).trim(), "0", "Stale day win must NOT increment streak");

        // Now solve actual 2026-11-15 daily target and verify win is credited strictly to 2026-11-15
        await submitChemleMoleculeByUi(page, day2.target.id, day2.target.name);
        assert.ok((await page.locator("#chemle-banner").textContent()).includes("Solved in 1/6!"));
        const statsAfterDay2 = await page.evaluate(() => JSON.parse(localStorage.getItem("biowordle_stats_v2")));
        assert.equal(statsAfterDay2.lastStreakIncrementUtc, "2026-11-15");
        assert.equal(statsAfterDay2.streak, 1);

        // Also verify automatic 1-second timer detects rollover to 2026-11-16 cleanly without a click
        await page.evaluate(() => {
          window.__BIOWORDLE_TEST_UTC_DATE__ = "2026-11-16";
        });
        await page.waitForTimeout(1300);
        assert.equal((await page.locator("#chemle-try-counter").textContent()).trim(), "Guesses: 0 / 6");
        assert.ok((await page.locator("#utc-rollover-banner").textContent()).includes("2026-11-16"));

        await page.screenshot({ path: path.join(screenshotsDir, "e2e_5_utc_rollover.png"), fullPage: true });
        assert.deepEqual(pageErrors, []);
        console.log("[E2E 5 PASS] UTC rollover guard & timer reset verified.");
      } finally {
        await context.close();
      }
    }

    // -----------------------------------------------------------------------
    // Standalone Bundle Smoke Check (dist/biowordle_standalone.html)
    // -----------------------------------------------------------------------
    {
      const context = await browser.newContext({ viewport: { width: 1400, height: 900 } });
      const page = await context.newPage();
      const pageErrors = [];
      page.on("pageerror", (err) => pageErrors.push(err.message));
      try {
        await page.goto(standaloneUrl);
        await page.waitForLoadState("domcontentloaded");
        assert.equal(await page.locator("#btn-hint").isDisabled(), true);
        assert.equal((await page.locator("#chemle-clue-stage").textContent()).trim(), "Level 0 · k-Anonymous");
        assert.deepEqual(pageErrors, []);
        console.log("[E2E Standalone Smoke PASS] dist/biowordle_standalone.html loaded cleanly.");
      } finally {
        await context.close();
      }
    }

    console.log("[Playwright E2E PASS] All 5 E2E scenarios + standalone bundle check passed.");
  } finally {
    await browser.close();
  }
}

run().catch((err) => {
  console.error("[Playwright E2E FAIL]", err);
  process.exit(1);
});
