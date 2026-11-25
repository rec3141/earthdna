#!/usr/bin/env python3
"""Serve a built globe to every screen of a multi-screen wall and keep the screens in step.

    python3 lg_server.py work/globe [--port 8765] [--viewer viewer/globe.html] [--nav auto|PATH|off]
    python3 lg_server.py --relay https://microscape.app/earthdna/sync/ROOM [--nav auto|PATH]

  GET  /...           files from the globe directory (/ is its index.html, or --viewer's file)
  GET  /sync/events   server-sent events: the leader's latest state (default event), 3D mouse axes ("nav")
  POST /sync/state    the leading screen's state, passed on to every screen and kept for screens that join later
  POST /sync/log      a screen's error message, written to this server's log

Each screen opens  http://SERVER:PORT/?sync&wall=COLS&screen=I  (one with &lead), see globe.html.

With --relay the screens sync through a hosted copy of the viewer instead (microscape.app's /earthdna pages with
?sync=ROOM): nothing is served here, and the 3D mouse's axes are posted to the relay's /nav.

The 3D mouse (3Dconnexion SpaceNavigator and relatives) is read from its Linux input device: the first device whose
name mentions 3Dconnexion, SpaceNavigator or SpaceMouse (or /dev/input/spacenavigator, which Liquid Galaxy's udev rule
makes), or the path given with --nav. Reading it needs read permission on that /dev/input/event* file. Its axes go out about 30 times a second while
it is pushed, normalised to roughly -1..1 (below DEAD counts as zero: the device rests slightly off zero), with one
all-zero message when it is let go.

Standard library only, Python 3.7+.
"""
import argparse
import json
import os
import queue
import re
import struct
import sys
import threading
import time
import http.client
import urllib.parse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

clients = set()            # one queue per connected screen
clients_lock = threading.Lock()
last_state = None          # the leader's latest state, bytes


def publish(event, data):
    msg = (("event: %s\n" % event) if event else "") + "data: %s\n\n" % data
    with clients_lock:
        for q in list(clients):
            try:
                q.put_nowait(msg)
            except queue.Full:  # a screen that stopped reading: drop it, its handler ends on the next write
                clients.discard(q)


class Handler(SimpleHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    viewer = None

    def log_message(self, fmt, *args):
        if "/sync/" not in self.path:
            super().log_message(fmt, *args)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store" if self.path.split("?")[0] in ("/", "/index.html") else "max-age=600")
        super().end_headers()

    def do_GET(self):
        path = self.path.split("?")[0]
        if path == "/sync/events":
            return self.events()
        if self.viewer and path in ("/", "/index.html"):
            body = open(self.viewer, "rb").read()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()

    def do_POST(self):
        global last_state
        path = self.path.split("?")[0]
        n = int(self.headers.get("Content-Length") or 0)
        body = self.rfile.read(n)
        if path == "/sync/log":  # a screen's error report
            print("screen error: %s" % body.decode("utf-8", "replace")[:2000], file=sys.stderr)
            self.send_response(204)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        if path != "/sync/state":
            self.send_error(404)
            return
        try:
            json.loads(body)
        except ValueError:
            self.send_error(400)
            return
        last_state = body
        publish(None, body.decode())
        self.send_response(204)
        self.send_header("Content-Length", "0")
        self.end_headers()

    def events(self):
        q = queue.Queue(maxsize=600)
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Connection", "close")
        self.end_headers()
        self.close_connection = True
        if last_state:
            q.put("data: %s\n\n" % last_state.decode())
        with clients_lock:
            clients.add(q)
        try:
            while True:
                try:
                    msg = q.get(timeout=15)
                except queue.Empty:
                    msg = ": ping\n\n"  # keeps proxies and idle sockets from closing the stream
                self.wfile.write(msg.encode())
                self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError, OSError):
            pass
        finally:
            with clients_lock:
                clients.discard(q)


# ---------- 3D mouse
EV_KEY, EV_REL, EV_ABS = 1, 2, 3
EVENT = struct.Struct("llHHi")  # struct input_event on 64-bit Linux
FULL = 350.0                    # SpaceNavigator axes run about -350..350
DEAD = 0.06                     # a released SpaceNavigator rests a few percent off zero: smaller values count as zero


def find_nav():
    if os.path.exists("/dev/input/spacenavigator"):
        return "/dev/input/spacenavigator"
    try:
        text = open("/proc/bus/input/devices").read()
    except OSError:
        return None
    for block in text.split("\n\n"):
        if re.search(r"3Dconnexion|SpaceNavigator|SpaceMouse|SpacePilot|SpaceExplorer", block, re.I):
            m = re.search(r"H: Handlers=.*?\b(event\d+)", block)
            if m:
                return "/dev/input/" + m.group(1)
    return None


class Relay:
    """posts 3D mouse axes to a hosted room's /nav over one kept-alive connection; a failed post is dropped"""

    def __init__(self, url):
        u = urllib.parse.urlsplit(url.rstrip("/"))
        self.https, self.host, self.path = u.scheme == "https", u.netloc, u.path + "/nav"
        self.conn, self.warned = None, False

    def __call__(self, event, data):
        if event != "nav":
            return
        for attempt in (0, 1):  # once more on a fresh connection if the kept one has gone
            try:
                if self.conn is None:
                    cls = http.client.HTTPSConnection if self.https else http.client.HTTPConnection
                    self.conn = cls(self.host, timeout=5)
                self.conn.request("POST", self.path, body=data, headers={"Content-Type": "application/json"})
                r = self.conn.getresponse()
                r.read()
                if r.status >= 400 and not self.warned:
                    print("relay: %s answered %d" % (self.host, r.status), file=sys.stderr)
                    self.warned = True
                return
            except (OSError, http.client.HTTPException) as e:
                self.conn = None
                if attempt and not self.warned:
                    print("relay: %s (axes dropped until it answers)" % e, file=sys.stderr)
                    self.warned = True


def read_nav(path, sink=publish):
    axes, buttons = [0.0] * 6, 0
    lock = threading.Lock()

    def sender():  # about 30 Hz while pushed; at rest, only a change (the release, a button)
        last = None
        while True:
            time.sleep(1 / 30)
            with lock:
                a, b = list(axes), buttons
            a = [round(v, 3) if abs(v) >= DEAD else 0.0 for v in a]
            if any(a) or (a, b) != last:
                sink("nav", json.dumps({"a": a, "b": b}))
                last = (a, b)

    threading.Thread(target=sender, daemon=True).start()
    while True:
        try:
            with open(path, "rb") as f:
                print("3D mouse: reading %s" % path, file=sys.stderr)
                while True:
                    data = f.read(EVENT.size)
                    if len(data) < EVENT.size:
                        break
                    _, _, typ, code, value = EVENT.unpack(data)
                    with lock:
                        if typ in (EV_REL, EV_ABS) and code < 6:
                            axes[code] = max(-1.0, min(1.0, value / FULL))
                        elif typ == EV_KEY and 256 <= code < 264:  # BTN_0..
                            bit = 1 << (code - 256)
                            buttons = (buttons | bit) if value else (buttons & ~bit)
        except PermissionError:
            print("3D mouse: no permission to read %s (a udev rule or the 'input' group can grant it); carrying on without it" % path, file=sys.stderr)
            return
        except OSError as e:
            print("3D mouse: %s; retrying in 5 s" % e, file=sys.stderr)
        with lock:
            axes[:] = [0.0] * 6
        time.sleep(5)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root", nargs="?", help="built globe directory (pipeline/build.sh's work/globe); not with --relay")
    ap.add_argument("--relay", help="a hosted room (https://HOST/earthdna/sync/ROOM): post the 3D mouse there, serve nothing")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--bind", default="0.0.0.0")
    ap.add_argument("--viewer", help="serve this file as / instead of ROOT/index.html (for working on the viewer)")
    ap.add_argument("--nav", default="auto", help="3D mouse input device: auto, a /dev/input path, or off")
    args = ap.parse_args()
    sys.stderr.reconfigure(line_buffering=True)  # a log file sees each line as it happens

    if bool(args.root) == bool(args.relay):
        ap.error("give either a globe directory or --relay")
    sink = Relay(args.relay) if args.relay else publish
    if args.nav != "off":
        path = find_nav() if args.nav == "auto" else args.nav
        if path:
            threading.Thread(target=read_nav, args=(path, sink), daemon=True).start()
        else:
            print("3D mouse: none found", file=sys.stderr)
    if args.relay:
        print("relaying the 3D mouse to %s" % args.relay, file=sys.stderr)
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            return
    Handler.viewer = os.path.abspath(args.viewer) if args.viewer else None
    server = ThreadingHTTPServer((args.bind, args.port), partial(Handler, directory=os.path.abspath(args.root)))
    server.daemon_threads = True
    print("serving %s on http://%s:%d/" % (args.root, args.bind, args.port), file=sys.stderr)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
