#!/usr/bin/env node
// Render a show (viewer/shows/*.json) frame by frame to a dome-master video or a PNG sequence.
//
//   node viewer/capture_show.js GLOBE_DIR SHOW.json OUT.mp4|OUT_DIR [--size 2048] [--fps 30] [--from S] [--to S] [--crf 14]
//
// GLOBE_DIR holds a packed globe (pack_globe.py's files). An OUT ending in .mp4/.mov/.mkv is encoded with ffmpeg
// (H.264, 4:2:0); anything else is taken as a directory for frame_000000.png ..., which is what most planetarium
// encoders want. Needs playwright (Chromium) and, for video, ffmpeg. Frames are drawn at show time, not wall-clock
// time, so the result is the same however slow the machine is.
const http = require("http"), fs = require("fs"), path = require("path"), { spawn } = require("child_process");
const { chromium } = require("playwright");

const args = process.argv.slice(2), opt = {}, pos = [];
for (let i = 0; i < args.length; i++) args[i].startsWith("--") ? (opt[args[i].slice(2)] = args[++i]) : pos.push(args[i]);
if (pos.length !== 3) { console.error("usage: capture_show.js GLOBE_DIR SHOW.json OUT.mp4|OUT_DIR [--size N] [--fps N] [--from S] [--to S] [--crf N]"); process.exit(2); }
const [globeDir, showFile, out] = pos, viewer = path.join(__dirname, "globe.html");

const server = http.createServer((req, res) => {
  const u = decodeURIComponent(req.url.split("?")[0]);
  const f = u === "/" ? viewer : u === "/__show.json" ? showFile : path.join(globeDir, path.basename(u));
  fs.readFile(f, (e, b) => {
    if (e) { res.writeHead(404); res.end(); return; }
    res.writeHead(200, { "content-type": f.endsWith(".html") ? "text/html; charset=utf-8" : f.endsWith(".json") ? "application/json" : "text/plain" });
    res.end(b);
  });
});

(async () => {
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  const show = JSON.parse(fs.readFileSync(showFile, "utf8"));
  const size = +(opt.size || show.size || 2048);
  const browser = await chromium.launch({ args: ["--ignore-gpu-blocklist", "--enable-unsafe-swiftshader"] });
  const page = await browser.newPage({ viewport: { width: size, height: size }, deviceScaleFactor: 1 });
  page.on("console", (m) => { if (m.type() === "warning" || m.type() === "error") console.error("page:", m.text()); });
  await page.goto(`http://127.0.0.1:${server.address().port}/?show=__show.json&capture`);
  const info = await page.evaluate(() => window.earthdna.ready());
  const fps = +(opt.fps || info.fps), t0 = +(opt.from || 0), t1 = +(opt.to || info.duration), n = Math.round((t1 - t0) * fps);
  const video = /\.(mp4|mov|mkv)$/i.test(out);
  let ff = null;
  if (video) {
    ff = spawn("ffmpeg", ["-y", "-loglevel", "error", "-f", "image2pipe", "-framerate", String(fps), "-c:v", "png", "-i", "-",
      "-c:v", "libx264", "-preset", "slow", "-crf", String(opt.crf || 14), "-pix_fmt", "yuv420p", "-movflags", "+faststart", out], { stdio: ["pipe", "inherit", "inherit"] });
  } else fs.mkdirSync(out, { recursive: true });
  const started = Date.now();
  for (let i = 0; i < n; i++) {
    await page.evaluate((t) => window.earthdna.seek(t), t0 + i / fps);
    const png = await page.screenshot({ type: "png" });
    if (ff) { if (!ff.stdin.write(png)) await new Promise((r) => ff.stdin.once("drain", r)); }
    else fs.writeFileSync(path.join(out, `frame_${String(i).padStart(6, "0")}.png`), png);
    if (i % fps === 0) process.stderr.write(`\r${i}/${n} frames, ${((Date.now() - started) / 1000 / (i + 1)).toFixed(2)} s each `);
  }
  if (ff) { ff.stdin.end(); await new Promise((r) => ff.on("close", r)); }
  process.stderr.write(`\r${n} frames (${size} x ${size}, ${fps} fps) -> ${out}\n`);
  await browser.close(); server.close();
})().catch((e) => { console.error(e); process.exit(1); });
