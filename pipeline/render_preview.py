#!/usr/bin/env python3
"""Still preview of a packed globe (pack_globe.py output, .bin or .b64.txt files) for picker cards and READMEs.

Orthographic view of one layout, points coloured by group with the viewer's palette (microscape-app globe.js:
the 24 PALETTE colours, then golden-angle hues), nearer points drawn over farther ones, limb darkened,
transparent outside the globe. Writes a PNG or WebP (supersampled, then downscaled).

usage: render_preview.py GLOBE_DIR out.png [--layout sculpt] [--lon L --lat B] [--size 640]
  without --lon/--lat the view faces the points' centre of mass
"""
import argparse, base64, colorsys, json, os
import numpy as np
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("globe"); ap.add_argument("out")
ap.add_argument("--layout", default="sculpt"); ap.add_argument("--lon", type=float)
ap.add_argument("--lat", type=float); ap.add_argument("--size", type=int, default=640)
ap.add_argument("--ss", type=int, default=2, help="supersampling factor")
a = ap.parse_args()

def load(name):
    p = os.path.join(a.globe, name)
    if os.path.exists(p + ".bin"): return np.fromfile(p + ".bin", np.uint8)
    return np.frombuffer(base64.b64decode(open(p + ".b64.txt", "rb").read()), np.uint8)

meta = json.load(open(os.path.join(a.globe, "globe.json")))
N = meta["n"]
pos = load("pos")
grp = pos[4 * N:5 * N]
ll = pos[:4 * N].view("<u2") if a.layout == "similarity" else load(f"layout_{a.layout}").view("<u2")
lon = ll[:N] / 65535 * 2 * np.pi - np.pi
lat = ll[N:2 * N] / 65535 * np.pi - np.pi / 2

PALETTE = ['#4e9bf0', '#f2994a', '#58c78a', '#e8595f', '#b28df2', '#f2d35a', '#41c8c8', '#f07ab8', '#9ccf4f', '#e3a37a',
           '#7f8cf2', '#c25ad6', '#5fb0e8', '#f5b8d0', '#a6e0b2', '#d9774a', '#8fd3f0', '#ccb25a', '#9b7bd1', '#6fd6a0',
           '#e07f99', '#bfc8d6', '#e8c06f', '#7ab8a8']
def hsl(h, s, l): return colorsys.hls_to_rgb(h / 360, l, s)
ng = len(meta["groups"])
cols = np.array([[int(c[i:i + 2], 16) / 255 for i in (1, 3, 5)] for c in PALETTE] +
                [hsl((g * 137.5) % 360, 0.32, 0.56) if g < ng else (0.5, 0.55, 0.62) for g in range(len(PALETTE), 256)])

# rotate so (lon0, lat0) faces the viewer: x right, y up, z toward the viewer
x, y, z = np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)
m = np.array([x.mean(), y.mean(), z.mean()])
l0 = np.radians(a.lon) if a.lon is not None else np.arctan2(m[1], m[0])
b0 = np.radians(a.lat) if a.lat is not None else np.arcsin(m[2] / np.linalg.norm(m))
xr, yr = x * np.cos(l0) + y * np.sin(l0), -x * np.sin(l0) + y * np.cos(l0)       # spin about the pole
depth, up = xr * np.cos(b0) + z * np.sin(b0), -xr * np.sin(b0) + z * np.cos(b0)   # tilt toward the viewer
right = yr
S = a.size * a.ss
R = 0.46 * S
vis = depth > 0
order = np.argsort(depth[vis])                                   # far first, near drawn last
px = (S / 2 + R * right[vis][order]).astype(int)
py = (S / 2 - R * up[vis][order]).astype(int)
shade = (0.45 + 0.55 * depth[vis][order])[:, None]                # limb darkening
c = cols[grp[vis][order]] * shade

bg = np.array([2, 6, 23]) / 255                                    # viewer --space #020617
img = np.tile(bg, (S, S, 1)).astype(np.float32)
# the sphere itself, dark blue, so gaps between points read as ocean
yy, xx = np.mgrid[0:S, 0:S]
rr = np.hypot(xx - S / 2, yy - S / 2) / R
ocean = rr < 1
img[ocean] = (np.array([9, 22, 48]) / 255) * (0.55 + 0.45 * np.sqrt(np.clip(1 - rr[ocean] ** 2, 0, 1)))[:, None]
r = max(1, a.ss)                                                   # point radius in supersampled pixels
for dx in range(-r + 1, r):
    for dy in range(-r + 1, r):
        if dx * dx + dy * dy >= r * r: continue
        qx, qy = np.clip(px + dx, 0, S - 1), np.clip(py + dy, 0, S - 1)
        img[qy, qx] = c                                            # later (nearer) points win
# faint atmosphere ring
ring = np.clip(1 - np.abs(rr - 1.0) * 28, 0, 1) * (rr > 0.97)
img += ring[..., None] * np.array([0.22, 0.45, 0.8]) * 0.5
# transparent outside the globe and its halo, so the image sits on any background
alpha = np.clip(np.maximum(ocean.astype(np.float32), ring * 1.6), 0, 1)
rgba = np.dstack([np.clip(img, 0, 1), alpha])
out = Image.fromarray((rgba * 255).astype(np.uint8), "RGBA").resize((a.size, a.size), Image.LANCZOS)
out.save(a.out, optimize=True)
print(f"{vis.sum():,} of {N:,} points on the near side -> {a.out}")
