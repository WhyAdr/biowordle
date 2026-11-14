// Reproducible local delivery benchmark. No changes to application artifacts.
// Usage: node tools/distribution_benchmark.mjs [runs=3]
import { chromium } from '@playwright/test';
import { execFileSync } from 'node:child_process';
import { createServer } from 'node:http';
import { gzipSync, brotliCompressSync, constants } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2);
const useWorking = args.includes('--working');
const positional = args.filter(a => !a.startsWith('--'));
const revision = useWorking ? 'working-tree' : '83f05f96102abea5d2fcce0ed75dbf3535127f71';
const runs = Number(positional[0] || 3);
if (!Number.isInteger(runs) || runs < 1 || runs > 20) throw Error('runs must be 1..20');
const files = ['index.html', 'styles.css', 'app.js', 'data/puzzles_data.js', 'dist/biowordle_standalone.html'];
const { readFileSync } = await import('node:fs');
const base = Object.fromEntries(files.map(f => [
  f,
  useWorking
    ? readFileSync(path.join(root, f))
    : execFileSync('git', ['show', `${revision}:${f}`], { cwd: root, maxBuffer: 16 * 1024 * 1024 }),
]));
const sizes = b => ({ raw: b.length, gzip: gzipSync(b, { level: 9 }).length, brotli: brotliCompressSync(b, { params: { [constants.BROTLI_PARAM_QUALITY]: 11 } }).length });
const payload = JSON.parse(base['data/puzzles_data.js'].toString().replace(/^[\s\S]*?window\.BIOWORDLE_DATA\s*=\s*/, '').replace(/;\s*$/, ''));
const compact = x => Buffer.from('window.BIOWORDLE_DATA = ' + JSON.stringify(x) + ';');
const stripped = JSON.parse(JSON.stringify(payload));
for (const row of [...stripped.chemle, ...stripped.seqle]) for (const k of Object.keys(row)) if (k.startsWith('svg_') && k !== 'svg_sha256') delete row[k];
const splitChem = { ...payload, seqle: [] };
const splitSeq = { ...payload, chemle: [], chemle_similarity: {} };
const projections = { original_data: sizes(base['data/puzzles_data.js']), without_svg_values_nonfunctional: sizes(compact(stripped)), chemle_only_nonfunctional: sizes(compact(splitChem)), seqle_only_nonfunctional: sizes(compact(splitSeq)), standalone: sizes(base['dist/biowordle_standalone.html']) };

const server = createServer((req, res) => {
  const [mode, ...parts] = new URL(req.url, 'http://localhost').pathname.slice(1).split('/');
  const key = parts.join('/') || 'index.html';
  if (!base[key] || !['raw', 'gzip', 'cache'].includes(mode)) { res.writeHead(404); res.end(); return; }
  const body = mode === 'raw' ? base[key] : gzipSync(base[key], { level: 9 });
  const type = key.endsWith('.css') ? 'text/css' : key.endsWith('.js') ? 'application/javascript' : 'text/html';
  res.writeHead(200, { 'Content-Type': `${type}; charset=utf-8`, 'Content-Length': body.length, 'Cache-Control': mode === 'cache' ? 'public, max-age=3600' : 'no-store', ...(mode !== 'raw' ? { 'Content-Encoding': 'gzip', 'Vary': 'Accept-Encoding' } : {}) });
  res.end(body);
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const port = server.address().port;
const browser = await chromium.launch({ headless: true });
const output = { revision, runs, conditions: 'Chromium 390x844 mobile, 4x CPU; HTTP Fast 3G 1.6 Mbps/150 ms; loopback server; file unthrottled network', projections_bytes: projections, samples: {} };
async function sample(mode) {
  const context = await browser.newContext({ viewport: { width: 390, height: 844 }, deviceScaleFactor: 2, isMobile: true, hasTouch: true });
  const page = await context.newPage();
  page.setDefaultTimeout(60000);
  const cdp = await context.newCDPSession(page);
  await cdp.send('Emulation.setCPUThrottlingRate', { rate: 4 });
  if (mode !== 'file') {
    await cdp.send('Network.enable');
    await cdp.send('Network.emulateNetworkConditions', { offline: false, latency: 150, downloadThroughput: Math.floor(1.6 * 1024 * 1024 / 8), uploadThroughput: Math.floor(750 * 1024 / 8) });
  }
  let failures = [], outside = [];
  page.on('pageerror', e => failures.push(e.message));
  page.on('request', r => { if (!r.url().startsWith(`http://127.0.0.1:${port}/`) && !r.url().startsWith('file:') && !r.url().startsWith('data:')) outside.push(r.url()); });
  const url = mode === 'file' ? pathToFileURL(path.join(root, 'dist/biowordle_standalone.html')).href : `http://127.0.0.1:${port}/${mode}/index.html`;
  const navigate = async () => {
    const start = performance.now();
    await page.goto(url, { waitUntil: 'load' });
    await page.locator('#chemle-search-input').waitFor({ state: 'visible' });
    const boot_ms = performance.now() - start;
    const perf = await page.evaluate(() => [...performance.getEntriesByType('navigation'), ...performance.getEntriesByType('resource')].map(x => ({ name: x.name, transfer: x.transferSize, encoded: x.encodedBodySize })));
    return { boot_ms, transfer_bytes: perf.reduce((a, x) => a + x.transfer, 0), encoded_body_bytes: perf.reduce((a, x) => a + x.encoded, 0) };
  };
  const cold = await navigate();
  const heap = await cdp.send('Runtime.getHeapUsage');
  const t0 = performance.now();
  await page.locator('#btn-playstyle').click();
  await page.locator('#badge-chemle-count').filter({ hasText: '78 Mols' }).waitFor();
  const switch_ms = performance.now() - t0;
  const t1 = performance.now();
  await page.locator('#chemle-search-input').fill('Xylene');
  await page.locator('#chemle-dropdown .chemle-opt[data-mol-id="m_xylene"]').waitFor({ state: 'visible' });
  const search_ms = performance.now() - t1;
  const t2 = performance.now();
  await page.locator('#chemle-dropdown .chemle-opt[data-mol-id="m_xylene"]').click();
  await page.locator('#chemle-guesses .chemle-guess-card').first().waitFor({ state: 'visible' });
  const guess_ms = performance.now() - t2;
  const warm = mode === 'cache' ? await navigate() : null;
  await context.close();
  return { cold, warm, switch_ms, search_ms, guess_ms, js_heap_used_bytes: heap.usedSize, page_errors: failures, external_requests: outside };
}
try {
  // file:// points at the working tree: require byte-identical bundle for honest comparison.
  const fileMatches = useWorking || readFileSync(path.join(root, 'dist/biowordle_standalone.html')).equals(base['dist/biowordle_standalone.html']);
  for (const mode of ['raw', 'gzip', 'cache', ...(fileMatches ? ['file'] : [])]) {
    output.samples[mode] = [];
    for (let i = 0; i < runs; i++) output.samples[mode].push(await sample(mode));
  }
  if (!fileMatches) output.file_note = 'Skipped: working standalone differs from baseline commit.';
  console.log(JSON.stringify(output, null, 2));
} finally {
  await browser.close();
  server.closeAllConnections?.();
  await new Promise(r => server.close(r));
}
