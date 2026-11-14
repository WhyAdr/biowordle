# BioWordle UI & Browser Testing

## Subagent Delegation
Delegate browser UI, click-flow, persistence, and E2E checks to `invoke_subagent` (`TypeName: "self"`, `Role: "UI Test Runner"`), which drives headless Microsoft Edge via Antigravity's bundled `playwright-core` without requiring the `chrome-devtools` plugin.
- **Done criteria**: run the script, verify assertions pass, save any screenshots under `tests/screenshots/` (gitignored), and return a concise pass/fail summary.

## Headless Edge + `playwright-core`
- **Module**: `createRequire(import.meta.url)(path.join(process.env.LOCALAPPDATA, 'Programs/Antigravity IDE/resources/app/node_modules/playwright-core'))` (no `npm i playwright` needed).
- **Launch**: `chromium.launch({ channel: 'msedge', headless: true, args: ['--no-sandbox', '--disable-gpu'] })` inside `try / finally { await browser.close(); }`.
- **Sandbox**: `run_command` launching `msedge.exe` requires `BypassSandbox: true` (Windows `--remote-debugging-pipe` cannot attach across the terminal sandbox).
- **Template**: run `node tests/test_playwright_e2e.mjs` (`BypassSandbox: true`).
