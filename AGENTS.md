# BioWordle

Interactive chemical (`Chemle`) and peptide motif (`Seqle`) deduction workspace with reproducible RDKit/UniProt compiler, cross-language parity tests, and Playwright Chromium E2E checks.

## Verification Commands
- **Python compiler & deterministic data generation**: `.\biowordle\.venv\Scripts\python.exe .\biowordle\build_puzzles.py` (normalizes tracked outputs to `LF`)
- **Python audit & unit tests**: `.\biowordle\.venv\Scripts\python.exe -B -m unittest discover -v -s biowordle -p "test_*.py"`
- **Solver difficulty calibration**: `.\biowordle\.venv\Scripts\python.exe .\biowordle\tools\solver_bot.py`
- **Online UniProt / PDB provenance check**: `.\biowordle\.venv\Scripts\python.exe .\biowordle\tools\verify_provenance_online.py` (`BypassSandbox: true`)
- **Phase 3 Seqle prototype benchmark**: `.\biowordle\.venv\Scripts\python.exe .\biowordle\tools\seqle_prototype.py`
- **Node golden vectors & stats ledger parity**: `node biowordle/tests/verify_golden_and_dom.mjs`
- **Portable Playwright Chromium E2E (`127.0.0.1` HTTP + `file://` standalone)**: `Push-Location biowordle; npm.cmd run test:e2e; Pop-Location` (`BypassSandbox: true` on Windows due to loopback/Chromium isolation)
