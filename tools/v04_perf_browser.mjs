// Isolated mobile-emulation benchmark; accepts folders emitted by v04_perf.py.
import { chromium } from '@playwright/test';
import { createServer } from 'node:http';
import { readFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';
import { pathToFileURL } from 'node:url';
import path from 'node:path';

const [root78, root120, repetitionsRaw = '2'] = process.argv.slice(2);
const repetitions = Number(repetitionsRaw);
if (!root78 || !root120 || !Number.isInteger(repetitions) || repetitions < 1) throw Error('Usage: node tools/v04_perf_browser.mjs <78-folder> <120-folder> [repetitions]');
const roots = { baseline78: root78, preview120: root120 };
const resources = {};
for (const [name, root] of Object.entries(roots)) {
  resources[name] = {};
  for (const f of ['index.html', 'styles.css', 'app.js', 'data/puzzles_data.js', 'dist/biowordle_standalone.html']) resources[name][f] = readFileSync(path.join(root, f));
}
const server = createServer((req, res) => {
  const [name, ...remaining] = new URL(req.url, 'http://localhost').pathname.slice(1).split('/');
  const f = remaining.join('/') || 'index.html';
  if (!resources[name]?.[f]) { res.writeHead(404); res.end(); return; }
  const b = gzipSync(resources[name][f], {level: 9});
  res.writeHead(200, {'Content-Type': f.endsWith('.js') ? 'application/javascript; charset=utf-8' : f.endsWith('.css') ? 'text/css; charset=utf-8' : 'text/html; charset=utf-8', 'Content-Encoding': 'gzip', 'Content-Length': b.length, 'Cache-Control':'no-store', 'Vary': 'Accept-Encoding'});
  res.end(b);
});
await new Promise(r => server.listen(0, '127.0.0.1', r));
const port = server.address().port;
const browser = await chromium.launch({headless: true});
async function run(name, mode) {
  const context = await browser.newContext({viewport:{width:390,height:844},deviceScaleFactor:2,isMobile:true,hasTouch:true});
  const page = await context.newPage();
  const cdp = await context.newCDPSession(page);
  await cdp.send('Emulation.setCPUThrottlingRate', {rate:4});
  if (mode === 'gzip_http') {
    await cdp.send('Network.enable');
    await cdp.send('Network.emulateNetworkConditions', {offline:false, latency:150, downloadThroughput:Math.floor(1.6*1024*1024/8), uploadThroughput:Math.floor(750*1024/8)});
  }
  const errors = [], external = [];
  page.on('pageerror', e => errors.push(e.message));
  page.on('request', r => {if (!r.url().startsWith(`http://127.0.0.1:${port}/`) && !/^(file|data|blob):/.test(r.url())) external.push(r.url());});
  const target = mode === 'gzip_http' ? `http://127.0.0.1:${port}/${name}/index.html?profile=learn` : `${pathToFileURL(path.join(roots[name], 'dist/biowordle_standalone.html')).href}?profile=learn`;
  const begin = performance.now();
  await page.goto(target, {waitUntil:'load',timeout:60000});
  await page.locator('#chemle-search-input').waitFor({state:'visible',timeout:60000});
  const ready_ms = performance.now()-begin;
  const heap = (await cdp.send('Runtime.getHeapUsage')).usedSize;
  const nav = await page.evaluate(() => [...performance.getEntriesByType('navigation'),...performance.getEntriesByType('resource')].map(x=>({transfer:x.transferSize,encoded:x.encodedBodySize})));
  const t = performance.now();
  await page.locator('#chemle-search-input').fill('Xylene');
  await page.locator('#chemle-dropdown .chemle-opt[data-mol-id="m_xylene"]').waitFor({state:'visible'});
  const search_ms = performance.now()-t;
  const q = performance.now();
  await page.locator('#chemle-dropdown .chemle-opt[data-mol-id="m_xylene"]').click();
  await page.locator('#chemle-guesses .chemle-guess-card').first().waitFor({state:'visible'});
  const guess_ms = performance.now()-q;
  const badge = await page.locator('#badge-chemle-count').innerText();
  await context.close();
  return {ready_ms,search_ms,guess_ms,js_heap_used_bytes:heap,transfer_bytes:mode==='gzip_http'?nav.reduce((s,x)=>s+x.transfer,0):null,badge,page_errors:errors,external_requests:external};
}
try {
  const results={conditions:'Chromium 390x844 mobile, 4x CPU, Fast 3G 1.6 Mbps/150 ms for gzip HTTP, unthrottled file',repetitions,samples:{}};
  for (const name of Object.keys(roots)) for (const mode of ['gzip_http','file']) {
    const key=`${name}_${mode}`;results.samples[key]=[];
    for (let i=0;i<repetitions;i++) results.samples[key].push(await run(name,mode));
  }
  process.stdout.write(JSON.stringify(results,null,2)+'\n');
} finally {
  await browser.close();server.closeAllConnections?.();await new Promise(r=>server.close(r));
}
