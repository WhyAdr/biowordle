/**
 * BioWordle (Seqle + Chemle) — Client Application Engine
 * Powered by the build-time RDKit puzzle database (`window.BIOWORDLE_DATA`).
 */
(function () {
  "use strict";

  const DATA = window.BIOWORDLE_DATA;
  if (!DATA) {
    console.error("BIOWORDLE_DATA not loaded. Run build_puzzles.py first.");
    return;
  }

  const MAX_GUESSES = 6;
  const AA_GROUPS = [
    { name: "Aliphatic", keys: ["A", "V", "I", "L", "M"] },
    { name: "Aromatic", keys: ["F", "Y", "W"] },
    { name: "Polar", keys: ["S", "T", "N", "Q"] },
    { name: "Positive (+)", keys: ["K", "R", "H"] },
    { name: "Negative (-)", keys: ["D", "E"] },
    { name: "Special", keys: ["C", "G", "P"] },
  ];

  // Application State
  const state = {
    mode: "chemle", // "chemle" | "seqle"
    playStyle: "daily", // "daily" | "infinite"
    stats: loadStats(),

    // Chemle state
    chemTarget: null,
    chemGuesses: [], // array of compound objects
    chemStatus: "playing", // "playing" | "won" | "lost"
    chemHintsUnlocked: 0,
    chemCategoryFilter: "All",

    // Seqle state
    seqTarget: null,
    seqGuesses: [], // array of { seq, evals, biochem }
    seqCurrentInput: "",
    seqStatus: "playing", // "playing" | "won" | "lost"
    seqHintsUnlocked: 0,
  };

  function loadStats() {
    try {
      const raw = localStorage.getItem("biowordle_stats_v1");
      if (raw) return JSON.parse(raw);
    } catch (_) {}
    return {
      chemSolved: 0,
      seqSolved: 0,
      streak: 0,
    };
  }

  function saveStats() {
    try {
      localStorage.setItem("biowordle_stats_v1", JSON.stringify(state.stats));
    } catch (_) {}
    renderHeaderStats();
  }

  function getDaySeed() {
    const now = new Date();
    return (
      now.getFullYear() * 10000 +
      (now.getMonth() + 1) * 100 +
      now.getDate()
    );
  }

  function seededIndex(seed, mod) {
    const x = Math.sin(seed * 9301 + 49297) * 233280;
    const frac = x - Math.floor(x);
    return Math.floor(frac * mod);
  }

  // -------------------------------------------------------------------------
  // Initialization & Mode Switching
  // -------------------------------------------------------------------------
  function init() {
    bindTopEvents();
    startNewChemlePuzzle(true);
    startNewSeqlePuzzle(true);
    switchMode("chemle");
    renderHeaderStats();
  }

  function bindTopEvents() {
    document.getElementById("tab-chemle").addEventListener("click", () => switchMode("chemle"));
    document.getElementById("tab-seqle").addEventListener("click", () => switchMode("seqle"));

    document.getElementById("btn-playstyle").addEventListener("click", () => {
      state.playStyle = state.playStyle === "daily" ? "infinite" : "daily";
      const btn = document.getElementById("btn-playstyle");
      btn.textContent = state.playStyle === "daily" ? "Mode: Daily Seed" : "Mode: Infinite Play";
      if (state.mode === "chemle") {
        startNewChemlePuzzle(state.playStyle === "daily");
      } else {
        startNewSeqlePuzzle(state.playStyle === "daily");
      }
    });

    document.getElementById("btn-new-puzzle").addEventListener("click", () => {
      state.playStyle = "infinite";
      document.getElementById("btn-playstyle").textContent = "Mode: Infinite Play";
      if (state.mode === "chemle") {
        startNewChemlePuzzle(false);
      } else {
        startNewSeqlePuzzle(false);
      }
    });

    document.getElementById("btn-hint").addEventListener("click", () => {
      if (state.mode === "chemle" && state.chemStatus === "playing") {
        state.chemHintsUnlocked = Math.min(3, state.chemHintsUnlocked + 1);
        renderChemleInspector();
      } else if (state.mode === "seqle" && state.seqStatus === "playing") {
        state.seqHintsUnlocked = Math.min(3, state.seqHintsUnlocked + 1);
        renderSeqleInspector();
      }
    });

    document.getElementById("btn-giveup").addEventListener("click", () => {
      if (state.mode === "chemle" && state.chemStatus === "playing") {
        state.chemStatus = "lost";
        renderChemleView();
      } else if (state.mode === "seqle" && state.seqStatus === "playing") {
        state.seqStatus = "lost";
        renderSeqleView();
      }
    });

    document.getElementById("btn-theme").addEventListener("click", () => {
      document.documentElement.classList.toggle("dark");
    });

    // Global physical keyboard listener for Seqle
    window.addEventListener("keydown", (e) => {
      if (state.mode !== "seqle" || state.seqStatus !== "playing") return;
      if (document.activeElement && ["INPUT", "SELECT", "TEXTAREA"].includes(document.activeElement.tagName)) {
        return;
      }
      if (e.key === "Backspace") {
        e.preventDefault();
        handleSeqleKey("BACKSPACE");
      } else if (e.key === "Enter") {
        e.preventDefault();
        handleSeqleKey("ENTER");
      } else if (/^[a-zA-Z]$/.test(e.key)) {
        const up = e.key.toUpperCase();
        if (DATA.aa_properties[up]) {
          e.preventDefault();
          handleSeqleKey(up);
        }
      }
    });

    setupChemleSearch();
    setupSeqleQuickPick();
  }

  function switchMode(mode) {
    state.mode = mode;
    document.getElementById("tab-chemle").classList.toggle("active", mode === "chemle");
    document.getElementById("tab-seqle").classList.toggle("active", mode === "seqle");
    document.getElementById("view-chemle").style.display = mode === "chemle" ? "grid" : "none";
    document.getElementById("view-seqle").style.display = mode === "seqle" ? "grid" : "none";
    if (mode === "chemle") {
      renderChemleView();
    } else {
      renderSeqleView();
    }
  }

  function renderHeaderStats() {
    document.getElementById("stat-chem").textContent = state.stats.chemSolved;
    document.getElementById("stat-seq").textContent = state.stats.seqSolved;
    document.getElementById("stat-streak").textContent = state.stats.streak;
  }

  // =========================================================================
  // CHEMLE ENGINE (Small Molecules, Metabolites & Hydrocarbon Intermediates)
  // =========================================================================
  function startNewChemlePuzzle(useDaily) {
    const pool = DATA.chemle;
    let idx;
    if (useDaily) {
      idx = seededIndex(getDaySeed() + 17, pool.length);
    } else {
      const currentId = state.chemTarget ? state.chemTarget.id : null;
      const candidates = pool.filter((c) => c.id !== currentId);
      idx = Math.floor(Math.random() * candidates.length);
      state.chemTarget = candidates[idx];
    }
    if (useDaily || !state.chemTarget) {
      state.chemTarget = pool[idx];
    }
    state.chemGuesses = [];
    state.chemStatus = "playing";
    state.chemHintsUnlocked = 0;
    const input = document.getElementById("chemle-search-input");
    if (input) input.value = "";
    renderChemleView();
  }

  function setupChemleSearch() {
    const input = document.getElementById("chemle-search-input");
    const dropdown = document.getElementById("chemle-dropdown");
    const catSelect = document.getElementById("chemle-cat-filter");

    // Populate categories
    const cats = ["All", ...Array.from(new Set(DATA.chemle.map((c) => c.Category)))];
    catSelect.innerHTML = cats
      .map((c) => `<option value="${escapeHtml(c)}">${escapeHtml(c)}</option>`)
      .join("");

    catSelect.addEventListener("change", () => {
      state.chemCategoryFilter = catSelect.value;
      updateDropdown();
    });

    function updateDropdown() {
      if (state.chemStatus !== "playing") {
        dropdown.classList.remove("open");
        return;
      }
      const q = input.value.trim().toLowerCase();
      const guessedIds = new Set(state.chemGuesses.map((g) => g.id));
      const matches = DATA.chemle.filter((c) => {
        if (guessedIds.has(c.id)) return false;
        if (state.chemCategoryFilter !== "All" && c.Category !== state.chemCategoryFilter) {
          return false;
        }
        if (!q) return true;
        return (
          c.name.toLowerCase().includes(q) ||
          c.formula.toLowerCase().includes(q) ||
          c.Category.toLowerCase().includes(q) ||
          c.functional_groups.some((fg) => fg.toLowerCase().includes(q))
        );
      });

      if (matches.length === 0) {
        dropdown.innerHTML = `<div class="chemle-option" style="color:var(--muted-foreground)">No matching unguessed compounds</div>`;
      } else {
        dropdown.innerHTML = matches
          .map(
            (c) => `
          <div class="chemle-option" data-id="${c.id}">
            <div>
              <strong>${escapeHtml(c.name)}</strong>
              <span style="color:var(--muted-foreground); font-size:11.5px; margin-left:6px;">${escapeHtml(c.formula)} · ${c.MW.toFixed(1)} Da</span>
            </div>
            <span class="badge badge-indigo">${escapeHtml(c.Category)}</span>
          </div>`
          )
          .join("");
      }
      dropdown.classList.add("open");

      dropdown.querySelectorAll(".chemle-option[data-id]").forEach((el) => {
        el.addEventListener("click", () => {
          const id = el.getAttribute("data-id");
          submitChemleGuess(id);
          input.value = "";
          dropdown.classList.remove("open");
        });
      });
    }

    input.addEventListener("focus", updateDropdown);
    input.addEventListener("input", updateDropdown);

    input.addEventListener("keydown", (e) => {
      if (e.key === "Enter") {
        const first = dropdown.querySelector(".chemle-option[data-id]");
        if (first) {
          e.preventDefault();
          submitChemleGuess(first.getAttribute("data-id"));
          input.value = "";
          dropdown.classList.remove("open");
        }
      } else if (e.key === "Escape") {
        dropdown.classList.remove("open");
      }
    });

    document.addEventListener("click", (e) => {
      if (!e.target.closest(".chemle-search-bar")) {
        dropdown.classList.remove("open");
      }
    });
  }

  function submitChemleGuess(compoundId) {
    if (state.chemStatus !== "playing") return;
    const mol = DATA.chemle.find((c) => c.id === compoundId);
    if (!mol) return;
    if (state.chemGuesses.some((g) => g.id === mol.id)) return;

    state.chemGuesses.push(mol);

    if (mol.id === state.chemTarget.id) {
      state.chemStatus = "won";
      state.stats.chemSolved += 1;
      state.stats.streak += 1;
      saveStats();
    } else if (state.chemGuesses.length >= MAX_GUESSES) {
      state.chemStatus = "lost";
      state.stats.streak = 0;
      saveStats();
    }
    renderChemleView();
  }

  function evaluateDelta(guessVal, targetVal, exactTol, closeTol) {
    const diff = guessVal - targetVal;
    const absDiff = Math.abs(diff);
    let cls = "match-far";
    if (absDiff <= exactTol) {
      cls = "match-exact";
    } else if (absDiff <= closeTol) {
      cls = "match-close";
    }
    let arrow = "✓";
    if (absDiff > exactTol) {
      arrow = targetVal > guessVal ? "↑" : "↓";
    }
    return { cls, arrow };
  }

  function renderChemleView() {
    const target = state.chemTarget;
    const banner = document.getElementById("chemle-banner");
    const guessesContainer = document.getElementById("chemle-guesses");
    const counterBadge = document.getElementById("chemle-try-counter");

    counterBadge.textContent = `Guesses: ${state.chemGuesses.length} / ${MAX_GUESSES}`;

    if (state.chemStatus === "won") {
      banner.className = "toast-banner show win";
      banner.innerHTML = `
        <span>Solved in ${state.chemGuesses.length}/${MAX_GUESSES}! The target molecule is <strong>${escapeHtml(target.name)}</strong> (${escapeHtml(target.formula)}).</span>
        <button class="btn btn-primary" id="chemle-play-again">Next Molecule</button>
      `;
      document.getElementById("chemle-play-again").onclick = () => startNewChemlePuzzle(false);
    } else if (state.chemStatus === "lost") {
      banner.className = "toast-banner show lose";
      banner.innerHTML = `
        <span>Out of guesses! The hidden molecule was <strong>${escapeHtml(target.name)}</strong> (${escapeHtml(target.formula)}).</span>
        <button class="btn btn-primary" id="chemle-play-again">Try Another</button>
      `;
      document.getElementById("chemle-play-again").onclick = () => startNewChemlePuzzle(false);
    } else {
      banner.className = "toast-banner";
      banner.innerHTML = "";
    }

    if (state.chemGuesses.length === 0) {
      guessesContainer.innerHTML = `
        <div class="hint-item" style="text-align:center; padding:22px 14px; color:var(--muted-foreground);">
          Search or pick any of the <strong>${DATA.chemle.length} RDKit-validated molecules</strong> above to make your first guess.<br>
          You'll get live feedback on <strong>MW, cLogP, TPSA, HBD, HBA, Rings, Aromaticity, Morgan Tanimoto Similarity (%)</strong>, and shared functional groups!
        </div>
      `;
    } else {
      const cardsHtml = state.chemGuesses
        .slice()
        .reverse()
        .map((g, revIdx) => {
          const guessNum = state.chemGuesses.length - revIdx;
          const sim =
            (DATA.chemle_similarity[g.id] && DATA.chemle_similarity[g.id][target.id]) || 0;

          const mwDelta = evaluateDelta(g.MW, target.MW, 5.0, 35.0);
          const logpDelta = evaluateDelta(g.cLogP, target.cLogP, 0.25, 1.0);
          const tpsaDelta = evaluateDelta(g.TPSA, target.TPSA, 5.0, 25.0);
          const hbdDelta = evaluateDelta(g.HBD, target.HBD, 0, 1);
          const hbaDelta = evaluateDelta(g.HBA, target.HBA, 0, 1);
          const ringsDelta = evaluateDelta(g.Rings, target.Rings, 0, 1);
          const aromDelta = evaluateDelta(g.AromaticRings, target.AromaticRings, 0, 1);

          const targetFGs = new Set(target.functional_groups);
          const fgHtml = g.functional_groups
            .map((fg) => {
              const isShared = targetFGs.has(fg);
              return `<span class="fg-chip ${isShared ? "shared" : ""}" title="${isShared ? "Present in Target!" : "Not in Target"}">${isShared ? "✓ " : ""}${escapeHtml(fg)}</span>`;
            })
            .join("");

          const simBadgeCls =
            sim === 100
              ? "badge-emerald"
              : sim >= 35
              ? "badge-amber"
              : "badge-indigo";

          return `
          <div class="chemle-guess-card">
            <div class="chemle-guess-top">
              <div class="chemle-guess-name">
                <span class="badge badge-indigo">#${guessNum}</span>
                <span>${escapeHtml(g.name)}</span>
                <span style="font-size:11.5px; color:var(--muted-foreground); font-weight:500;">${escapeHtml(g.formula)}</span>
              </div>
              <div style="display:flex; gap:6px; align-items:center;">
                <span class="badge ${ g.Category === target.Category ? "badge-emerald" : "badge-cyan" }">
                  ${g.Category === target.Category ? "✓ " : ""}${escapeHtml(g.Category)}
                </span>
                <span class="badge ${simBadgeCls}">Tanimoto Sim: ${sim}%</span>
              </div>
            </div>

            <div class="chemle-Compare-grid">
              <div class="mini-svg-well" title="${escapeHtml(g.canonical_smiles)}">
                ${g.svg_full}
              </div>
              <div class="delta-cell ${mwDelta.cls}">
                <span class="prop-label">MW (Da)</span>
                <span class="prop-val">${g.MW.toFixed(1)} <span class="delta-arrow">${mwDelta.arrow}</span></span>
              </div>
              <div class="delta-cell ${logpDelta.cls}">
                <span class="prop-label">cLogP</span>
                <span class="prop-val">${g.cLogP.toFixed(2)} <span class="delta-arrow">${logpDelta.arrow}</span></span>
              </div>
              <div class="delta-cell ${tpsaDelta.cls}">
                <span class="prop-label">TPSA (Å²)</span>
                <span class="prop-val">${g.TPSA.toFixed(1)} <span class="delta-arrow">${tpsaDelta.arrow}</span></span>
              </div>
              <div class="delta-cell ${hbdDelta.cls}">
                <span class="prop-label">HBD</span>
                <span class="prop-val">${g.HBD} <span class="delta-arrow">${hbdDelta.arrow}</span></span>
              </div>
              <div class="delta-cell ${hbaDelta.cls}">
                <span class="prop-label">HBA</span>
                <span class="prop-val">${g.HBA} <span class="delta-arrow">${hbaDelta.arrow}</span></span>
              </div>
              <div class="delta-cell ${ringsDelta.cls}">
                <span class="prop-label">Rings</span>
                <span class="prop-val">${g.Rings} <span class="delta-arrow">${ringsDelta.arrow}</span></span>
              </div>
              <div class="delta-cell ${aromDelta.cls}">
                <span class="prop-label">Arom Rings</span>
                <span class="prop-val">${g.AromaticRings} <span class="delta-arrow">${aromDelta.arrow}</span></span>
              </div>
            </div>

            ${
              fgHtml
                ? `<div class="fg-chips" style="margin-top:8px;"><span style="font-size:10px; color:var(--muted-foreground); align-self:center; margin-right:4px;">SMARTS Groups:</span>${fgHtml}</div>`
                : ""
            }
          </div>`;
        })
        .join("");
      guessesContainer.innerHTML = cardsHtml;
    }

    renderChemleInspector();
  }

  function renderChemleInspector() {
    const target = state.chemTarget;
    const well = document.getElementById("chemle-structure-well");
    const statusBadge = document.getElementById("chemle-clue-stage");
    const hintsBox = document.getElementById("chemle-hints-list");
    const fullPropsBox = document.getElementById("chemle-full-props");

    const tries = state.chemGuesses.length;
    const effectiveClueLevel =
      state.chemStatus !== "playing"
        ? 4
        : Math.max(state.chemHintsUnlocked, Math.floor(tries / 2));

    // Render 2D SVG well with appropriate progressive blur or scaffold
    if (state.chemStatus !== "playing") {
      well.className = "chem-svg-well";
      well.innerHTML = `${target.svg_full}<span class="well-overlay-tag">RDKit ACS-1996 2D</span>`;
      statusBadge.textContent = "Revealed";
      statusBadge.className = "badge badge-emerald";
    } else if (effectiveClueLevel >= 2 && target.has_scaffold) {
      well.className = "chem-svg-well";
      well.innerHTML = `${target.svg_scaffold}<span class="well-overlay-tag">Bemis-Murcko Scaffold Core</span>`;
      statusBadge.textContent = "Scaffold Core Unlocked";
      statusBadge.className = "badge badge-amber";
    } else if (effectiveClueLevel === 1) {
      well.className = "chem-svg-well partial-blur";
      well.innerHTML = `${target.svg_scaffold}<span class="well-overlay-tag">Blurred Core Preview</span>`;
      statusBadge.textContent = "Clue Level 1";
      statusBadge.className = "badge badge-cyan";
    } else {
      well.className = "chem-svg-well silhouette";
      well.innerHTML = `${target.svg_full}<span class="well-overlay-tag">Mystery Structure (Unlock via Guesses or Hint)</span>`;
      statusBadge.textContent = "Hidden";
      statusBadge.className = "badge badge-indigo";
    }

    // Progressive textual clues
    const clues = [];
    clues.push(`
      <div class="hint-item">
        <strong>Clue 1 · Chemical Space & Topology</strong>
        Category: <span class="badge badge-indigo">${escapeHtml(target.Category)}</span> ·
        Heavy Atoms: <b>${target.HeavyAtoms}</b> · Fsp³: <b>${target.Fsp3.toFixed(2)}</b> ·
        CNS MPO: <b>${target.CNS_MPO.toFixed(2)}</b>
      </div>
    `);

    if (effectiveClueLevel >= 1 || state.chemStatus !== "playing") {
      clues.push(`
        <div class="hint-item">
          <strong>Clue 2 · Biochemical Role & Mechanism</strong>
          ${escapeHtml(target.hint_role)}
        </div>
      `);
    }
    if (effectiveClueLevel >= 2 || state.chemStatus !== "playing") {
      clues.push(`
        <div class="hint-item">
          <strong>Clue 3 · Biosynthesis / Pathway Origin</strong>
          ${escapeHtml(target.hint_origin)}
        </div>
      `);
    }
    if (effectiveClueLevel >= 3 || state.chemStatus !== "playing") {
      clues.push(`
        <div class="hint-item">
          <strong>Clue 4 · Exact Molecular Formula</strong>
          Formula: <code>${escapeHtml(target.formula)}</code> · Exact Mass: <b>${target.exact_mass.toFixed(2)} Da</b>
        </div>
      `);
    }
    hintsBox.innerHTML = clues.join("");

    // Full canonical descriptor card when game finishes
    if (state.chemStatus !== "playing") {
      fullPropsBox.style.display = "block";
      fullPropsBox.innerHTML = `
        <div style="margin-top:12px; border-top:1px solid var(--border); padding-top:10px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:6px;">
            <strong style="font-size:12.5px;">${escapeHtml(target.name)} (${escapeHtml(target.formula)})</strong>
            <span class="badge ${target.PAINS_Alerts.length === 0 ? "badge-emerald" : "badge-rose"}">
              PAINS: ${target.PAINS_Alerts.length === 0 ? "Clean" : target.PAINS_Alerts.join(", ")}
            </span>
          </div>
          <div class="prop-grid">
            <div class="prop-cell"><span class="prop-label">MW</span><span class="prop-val">${target.MW.toFixed(1)}</span></div>
            <div class="prop-cell"><span class="prop-label">cLogP</span><span class="prop-val">${target.cLogP.toFixed(2)}</span></div>
            <div class="prop-cell"><span class="prop-label">TPSA</span><span class="prop-val">${target.TPSA.toFixed(1)}</span></div>
            <div class="prop-cell"><span class="prop-label">HBD / HBA</span><span class="prop-val">${target.HBD} / ${target.HBA}</span></div>
            <div class="prop-cell"><span class="prop-label">RotB</span><span class="prop-val">${target.RotB}</span></div>
            <div class="prop-cell"><span class="prop-label">Fsp3</span><span class="prop-val">${target.Fsp3.toFixed(2)}</span></div>
            <div class="prop-cell"><span class="prop-label">QED</span><span class="prop-val">${target.QED.toFixed(2)}</span></div>
            <div class="prop-cell"><span class="prop-label">CNS MPO</span><span class="prop-val">${target.CNS_MPO.toFixed(2)}</span></div>
            <div class="prop-cell"><span class="prop-label">SA Score</span><span class="prop-val">${target.SA_Score.toFixed(2)}</span></div>
          </div>
          <div style="margin-top:8px; font-size:11px; color:var(--muted-foreground); word-break:break-all;">
            <strong>Canonical SMILES:</strong> <code>${escapeHtml(target.canonical_smiles)}</code>
          </div>
        </div>
      `;
    } else {
      fullPropsBox.style.display = "none";
      fullPropsBox.innerHTML = "";
    }
  }

  // =========================================================================
  // SEQLE ENGINE (5-Residue Functional Motif Wordle)
  // =========================================================================
  function startNewSeqlePuzzle(useDaily) {
    const pool = DATA.seqle;
    let idx;
    if (useDaily) {
      idx = seededIndex(getDaySeed() + 99, pool.length);
    } else {
      const currentSeq = state.seqTarget ? state.seqTarget.seq : null;
      const candidates = pool.filter((m) => m.seq !== currentSeq);
      idx = Math.floor(Math.random() * candidates.length);
      state.seqTarget = candidates[idx];
    }
    if (useDaily || !state.seqTarget) {
      state.seqTarget = pool[idx];
    }
    state.seqGuesses = [];
    state.seqCurrentInput = "";
    state.seqStatus = "playing";
    state.seqHintsUnlocked = 0;
    renderSeqleView();
  }

  function setupSeqleQuickPick() {
    const sel = document.getElementById("seqle-motif-picker");
    sel.innerHTML =
      `<option value="">— Quick-fill a known 5-AA motif (or type any 5 amino acids) —</option>` +
      DATA.seqle
        .map(
          (m) =>
            `<option value="${m.seq}">${m.seq} · ${escapeHtml(m.title)} (${escapeHtml(m.family)})</option>`
        )
        .join("");

    sel.addEventListener("change", () => {
      if (sel.value && state.seqStatus === "playing") {
        state.seqCurrentInput = sel.value;
        sel.value = "";
        renderSeqleBoard();
      }
    });

    document.getElementById("seqle-btn-backspace").addEventListener("click", () => {
      handleSeqleKey("BACKSPACE");
    });
    document.getElementById("seqle-btn-submit").addEventListener("click", () => {
      handleSeqleKey("ENTER");
    });
  }

  function handleSeqleKey(key) {
    if (state.seqStatus !== "playing") return;
    if (key === "BACKSPACE") {
      state.seqCurrentInput = state.seqCurrentInput.slice(0, -1);
      renderSeqleBoard();
    } else if (key === "ENTER") {
      if (state.seqCurrentInput.length === 5) {
        submitSeqleGuess(state.seqCurrentInput);
      }
    } else if (state.seqCurrentInput.length < 5 && DATA.aa_properties[key]) {
      state.seqCurrentInput += key;
      renderSeqleBoard();
    }
  }

  function computeLivePeptideBiochem(seq) {
    const props = DATA.aa_properties;
    let gravy = 0;
    let charge = 0;
    for (let i = 0; i < seq.length; i++) {
      const aa = seq[i];
      gravy += props[aa].hydropathy;
      charge += props[aa].charge;
    }
    gravy = gravy / seq.length;

    // Henderson-Hasselbalch pI bisection matching build_puzzles.py
    const pka_nterm = 8.0;
    const pka_cterm = 3.1;
    const pka_side = {
      D: [3.9, -1],
      E: [4.1, -1],
      C: [8.3, -1],
      Y: [10.1, -1],
      H: [6.0, 1],
      K: [10.5, 1],
      R: [12.5, 1],
    };
    function netQ(ph) {
      let q =
        1.0 / (1.0 + Math.pow(10, ph - pka_nterm)) -
        1.0 / (1.0 + Math.pow(10, pka_cterm - ph));
      for (let i = 0; i < seq.length; i++) {
        const aa = seq[i];
        if (pka_side[aa]) {
          const [pk, sign] = pka_side[aa];
          if (sign > 0) q += 1.0 / (1.0 + Math.pow(10, ph - pk));
          else q -= 1.0 / (1.0 + Math.pow(10, pk - ph));
        }
      }
      return q;
    }
    let lo = 1.0,
      hi = 13.5;
    for (let iter = 0; iter < 25; iter++) {
      const mid = 0.5 * (lo + hi);
      if (netQ(mid) > 0) lo = mid;
      else hi = mid;
    }
    return {
      gravy: Number(gravy.toFixed(2)),
      net_charge_ph74: Number(charge.toFixed(1)),
      isoelectric_point: Number((0.5 * (lo + hi)).toFixed(2)),
    };
  }

  function evaluateSeqleGuess(guessSeq, targetSeq) {
    // 1. Standard Wordle exact + misplaced pass
    const evals = new Array(5).fill("absent");
    const targetRemaining = targetSeq.split("");

    // Exact matches first
    for (let i = 0; i < 5; i++) {
      if (guessSeq[i] === targetSeq[i]) {
        evals[i] = "exact";
        targetRemaining[i] = null;
      }
    }

    // Misplaced matches second
    for (let i = 0; i < 5; i++) {
      if (evals[i] !== "absent") continue;
      const foundIdx = targetRemaining.indexOf(guessSeq[i]);
      if (foundIdx !== -1) {
        evals[i] = "present";
        targetRemaining[foundIdx] = null;
      }
    }

    // Biochemical class match third (same amino acid group at position i!)
    for (let i = 0; i < 5; i++) {
      if (evals[i] === "absent") {
        const gGroup = DATA.aa_properties[guessSeq[i]].group;
        const tGroup = DATA.aa_properties[targetSeq[i]].group;
        if (gGroup === tGroup) {
          evals[i] = "class";
        }
      }
    }

    return evals;
  }

  function submitSeqleGuess(seq) {
    const target = state.seqTarget;
    const evals = evaluateSeqleGuess(seq, target.seq);
    const biochem = computeLivePeptideBiochem(seq);

    state.seqGuesses.push({ seq, evals, biochem });
    state.seqCurrentInput = "";

    if (seq === target.seq) {
      state.seqStatus = "won";
      state.stats.seqSolved += 1;
      state.stats.streak += 1;
      saveStats();
    } else if (state.seqGuesses.length >= MAX_GUESSES) {
      state.seqStatus = "lost";
      state.stats.streak = 0;
      saveStats();
    }
    renderSeqleView();
  }

  function renderSeqleView() {
    const target = state.seqTarget;
    const banner = document.getElementById("seqle-banner");
    const counterBadge = document.getElementById("seqle-try-counter");
    counterBadge.textContent = `Guesses: ${state.seqGuesses.length} / ${MAX_GUESSES}`;

    if (state.seqStatus === "won") {
      banner.className = "toast-banner show win";
      banner.innerHTML = `
        <span>Motif solved in ${state.seqGuesses.length}/${MAX_GUESSES}! <strong>${target.seq}</strong> — ${escapeHtml(target.title)}.</span>
        <button class="btn btn-primary" id="seqle-play-again">Next Motif</button>
      `;
      document.getElementById("seqle-play-again").onclick = () => startNewSeqlePuzzle(false);
    } else if (state.seqStatus === "lost") {
      banner.className = "toast-banner show lose";
      banner.innerHTML = `
        <span>Out of tries! The target motif was <strong>${target.seq}</strong> (${escapeHtml(target.title)}).</span>
        <button class="btn btn-primary" id="seqle-play-again">Try Another Motif</button>
      `;
      document.getElementById("seqle-play-again").onclick = () => startNewSeqlePuzzle(false);
    } else {
      banner.className = "toast-banner";
      banner.innerHTML = "";
    }

    renderSeqleBoard();
    renderSeqleKeyboard();
    renderSeqleInspector();
  }

  function renderSeqleBoard() {
    const board = document.getElementById("seqle-board");
    const target = state.seqTarget;
    const rowsHtml = [];

    for (let r = 0; r < MAX_GUESSES; r++) {
      const guess = state.seqGuesses[r];
      const isCurrentRow = r === state.seqGuesses.length && state.seqStatus === "playing";
      const tilesHtml = [];

      for (let c = 0; c < 5; c++) {
        let letter = "";
        let sub = "";
        let cls = "seqle-tile";

        if (guess) {
          letter = guess.seq[c];
          sub = DATA.aa_properties[letter].three;
          cls += ` state-${guess.evals[c]}`;
        } else if (isCurrentRow && state.seqCurrentInput[c]) {
          letter = state.seqCurrentInput[c];
          sub = DATA.aa_properties[letter].three;
          cls += " filled";
        }

        tilesHtml.push(`
          <div class="${cls}">
            <span>${letter}</span>
            ${sub ? `<span class="tile-sub">${sub}</span>` : ""}
          </div>
        `);
      }

      let stripHtml = "";
      if (guess) {
        const b = guess.biochem;
        const qDelta = evaluateDelta(b.net_charge_ph74, target.net_charge_ph74, 0.05, 1.0);
        const gDelta = evaluateDelta(b.gravy, target.gravy, 0.15, 0.8);
        const piDelta = evaluateDelta(b.isoelectric_point, target.isoelectric_point, 0.25, 1.2);
        const qStr = (b.net_charge_ph74 > 0 ? "+" : "") + b.net_charge_ph74.toFixed(1);

        stripHtml = `
          <div class="seqle-biochem-strip">
            <div class="strip-metric">
              <span class="prop-label">Charge</span>
              <span class="prop-val" style="color:${qDelta.cls === "match-exact" ? "var(--pal-emerald)" : "var(--foreground)"}">${qStr} ${qDelta.arrow}</span>
            </div>
            <div class="strip-metric">
              <span class="prop-label">GRAVY</span>
              <span class="prop-val" style="color:${gDelta.cls === "match-exact" ? "var(--pal-emerald)" : "var(--foreground)"}">${b.gravy.toFixed(2)} ${gDelta.arrow}</span>
            </div>
            <div class="strip-metric">
              <span class="prop-label">pI</span>
              <span class="prop-val" style="color:${piDelta.cls === "match-exact" ? "var(--pal-emerald)" : "var(--foreground)"}">${b.isoelectric_point.toFixed(1)} ${piDelta.arrow}</span>
            </div>
          </div>
        `;
      } else {
        stripHtml = `
          <div class="seqle-biochem-strip" style="opacity:0.45;">
            <div class="strip-metric"><span class="prop-label">Charge</span><span class="prop-val">—</span></div>
            <div class="strip-metric"><span class="prop-label">GRAVY</span><span class="prop-val">—</span></div>
            <div class="strip-metric"><span class="prop-label">pI</span><span class="prop-val">—</span></div>
          </div>
        `;
      }

      rowsHtml.push(`
        <div class="seqle-row">
          ${tilesHtml.join("")}
          ${stripHtml}
        </div>
      `);
    }

    board.innerHTML = rowsHtml.join("");
  }

  function renderSeqleKeyboard() {
    const kb = document.getElementById("seqle-keyboard");
    const rank = { absent: 1, class: 2, present: 3, exact: 4 };
    const keyStates = {};

    state.seqGuesses.forEach((g) => {
      for (let i = 0; i < 5; i++) {
        const aa = g.seq[i];
        const ev = g.evals[i];
        if (!keyStates[aa] || rank[ev] > rank[keyStates[aa]]) {
          keyStates[aa] = ev;
        }
      }
    });

    kb.innerHTML = AA_GROUPS.map((grp) => {
      const keysHtml = grp.keys
        .map((aa) => {
          const info = DATA.aa_properties[aa];
          const st = keyStates[aa] ? `key-${keyStates[aa]}` : "";
          return `
            <button class="aa-key ${st}" data-aa="${aa}" title="${info.name} (${info.three}) · Group: ${info.group} · Hydropathy: ${info.hydropathy}">
              <span>${aa}</span>
              <small>${info.three}</small>
            </button>
          `;
        })
        .join("");

      return `
        <div class="aa-group-box">
          <div class="aa-group-title"><span>${escapeHtml(grp.name)}</span></div>
          <div class="aa-keys-row">${keysHtml}</div>
        </div>
      `;
    }).join("");

    kb.querySelectorAll(".aa-key[data-aa]").forEach((btn) => {
      btn.addEventListener("click", () => {
        handleSeqleKey(btn.getAttribute("data-aa"));
      });
    });
  }

  function renderSeqleInspector() {
    const target = state.seqTarget;
    const well = document.getElementById("seqle-structure-well");
    const statusBadge = document.getElementById("seqle-clue-stage");
    const hintsBox = document.getElementById("seqle-hints-list");
    const fullPropsBox = document.getElementById("seqle-full-props");

    const tries = state.seqGuesses.length;
    const effectiveClueLevel =
      state.seqStatus !== "playing"
        ? 4
        : Math.max(state.seqHintsUnlocked, Math.floor(tries / 2));

    if (state.seqStatus !== "playing") {
      well.className = "chem-svg-well";
      well.innerHTML = `${target.svg_peptide}<span class="well-overlay-tag">RDKit Pentapeptide (${target.seq})</span>`;
      statusBadge.textContent = "Revealed";
      statusBadge.className = "badge badge-emerald";
    } else if (effectiveClueLevel >= 2) {
      well.className = "chem-svg-well partial-blur";
      well.innerHTML = `${target.svg_peptide}<span class="well-overlay-tag">2D Peptide Backbone Preview</span>`;
      statusBadge.textContent = "Clue Level 2";
      statusBadge.className = "badge badge-amber";
    } else {
      well.className = "chem-svg-well silhouette";
      well.innerHTML = `${target.svg_peptide}<span class="well-overlay-tag">Mystery 5-AA Peptide Structure</span>`;
      statusBadge.textContent = "Hidden";
      statusBadge.className = "badge badge-indigo";
    }

    const clues = [];
    clues.push(`
      <div class="hint-item">
        <strong>Clue 1 · Motif Family & Peptide Mass</strong>
        Family: <span class="badge badge-indigo">${escapeHtml(target.family)}</span> ·
        Pentapeptide MW: <b>${target.MW.toFixed(1)} Da</b> · Aromatic Residues: <b>${target.aromatic_count}</b>
      </div>
    `);

    if (effectiveClueLevel >= 1 || state.seqStatus !== "playing") {
      clues.push(`
        <div class="hint-item">
          <strong>Clue 2 · Biological Function</strong>
          ${escapeHtml(target.hint_role)}
        </div>
      `);
    }

    if (effectiveClueLevel >= 2 || state.seqStatus !== "playing") {
      clues.push(`
        <div class="hint-item">
          <strong>Clue 3 · Structural & Biochemical Mechanism</strong>
          ${escapeHtml(target.hint_detail)}
        </div>
      `);
    }

    if (effectiveClueLevel >= 3 || state.seqStatus !== "playing") {
      clues.push(`
        <div class="hint-item">
          <strong>Clue 4 · Residue Class Pattern</strong>
          Classes: <code>${target.residue_groups.join(" → ")}</code>
        </div>
      `);
    }

    hintsBox.innerHTML = clues.join("");

    if (state.seqStatus !== "playing") {
      fullPropsBox.style.display = "block";
      fullPropsBox.innerHTML = `
        <div style="margin-top:12px; border-top:1px solid var(--border); padding-top:10px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:6px;">
            <strong style="font-size:12.5px;">${target.seq} — ${escapeHtml(target.title)}</strong>
            <span class="badge badge-emerald">${escapeHtml(target.pattern)}</span>
          </div>
          <div class="prop-grid">
            <div class="prop-cell"><span class="prop-label">Formula</span><span class="prop-val">${escapeHtml(target.formula)}</span></div>
            <div class="prop-cell"><span class="prop-label">MW (Da)</span><span class="prop-val">${target.MW.toFixed(1)}</span></div>
            <div class="prop-cell"><span class="prop-label">Net Charge</span><span class="prop-val">${(target.net_charge_ph74 > 0 ? "+" : "") + target.net_charge_ph74.toFixed(1)}</span></div>
            <div class="prop-cell"><span class="prop-label">GRAVY</span><span class="prop-val">${target.gravy.toFixed(2)}</span></div>
            <div class="prop-cell"><span class="prop-label">pI</span><span class="prop-val">${target.isoelectric_point.toFixed(2)}</span></div>
            <div class="prop-cell"><span class="prop-label">TPSA (Å²)</span><span class="prop-val">${target.TPSA.toFixed(1)}</span></div>
          </div>
        </div>
      `;
    } else {
      fullPropsBox.style.display = "none";
      fullPropsBox.innerHTML = "";
    }
  }

  function escapeHtml(str) {
    return String(str ?? "")
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

