#!/bin/bash
# Put a built globe on a Liquid Galaxy wall: the sync server (lg_server.py) on this machine, and a full-screen
# Firefox on every screen showing its slice of the picture.
#
#   lg_launch.sh GLOBE_DIR       start (restarts anything this script started before)
#   lg_launch.sh stop            close the browsers and the server
#
# Settings (environment):
#   LG_ORDER   screens left to right, as ssh host names       default "lg5 lg6 lg7 lg8 lg1 lg2 lg3 lg4"
#   LG_LEAD    the screen that leads at start, shows the controls and takes the 3D mouse    default: this host
#   LG_SERVER  address the screens reach this machine on      default 10.42.42.1
#   LG_PORT    default 8765          LG_BEZEL   css pixels hidden behind each bezel, default 0
#   LG_EXTRA   more query parameters for every screen, e.g. "dome=30&navspeed=1.5"
#
# Screens that are switched off are skipped; the others still show their own slice.
# Google Earth keeps running underneath; the browsers cover it.
set -u
here=$(cd "$(dirname "$0")" && pwd)
ORDER=${LG_ORDER:-"lg5 lg6 lg7 lg8 lg1 lg2 lg3 lg4"}
LEAD=${LG_LEAD:-$(hostname -s)}
SERVER=${LG_SERVER:-10.42.42.1}
PORT=${LG_PORT:-8765}
BEZEL=${LG_BEZEL:-0}
EXTRA=${LG_EXTRA:-}
PROFILE='$HOME/.earthdna-wall'  # expanded on each screen's machine
SSH="ssh -o BatchMode=yes -o ConnectTimeout=4"

# run a command on a screen's machine (locally for this one)
on() { local host=$1; shift; if [ "$host" = "$(hostname -s)" ]; then bash -c "$*"; else $SSH "$host" "$*"; fi; }

# the bracket keeps pkill from matching the shell that runs it
STOP_BROWSER="pkill -f '[e]arthdna-wall' ; true"

stop_all() {
  for h in $ORDER; do on "$h" "$STOP_BROWSER" 2>/dev/null & done; wait
  pkill -f "[l]g_server.py.*--port $PORT" ; true
}

if [ "${1:-}" = stop ]; then stop_all; echo "stopped"; exit 0; fi
[ $# -ge 1 ] && [ -f "$1/globe.json" ] || { echo "usage: $0 GLOBE_DIR | stop   (GLOBE_DIR holds globe.json)" >&2; exit 2; }
globe=$(cd "$1" && pwd)

stop_all; sleep 1
nohup python3 "$here/lg_server.py" "$globe" --port "$PORT" --viewer "$here/globe.html" > "$HOME/.earthdna-wall-server.log" 2>&1 &
for _ in $(seq 20); do curl -s -o /dev/null "http://127.0.0.1:$PORT/globe.json" && break; sleep 0.5; done
curl -s -o /dev/null "http://127.0.0.1:$PORT/globe.json" || { echo "the server did not start; see ~/.earthdna-wall-server.log" >&2; exit 1; }

# Firefox profile: WebGL on whatever the card is, no first-run pages, no restore prompts
PREFS='user_pref("webgl.force-enabled", true);
user_pref("browser.shell.checkDefaultBrowser", false);
user_pref("browser.sessionstore.resume_from_crash", false);
user_pref("browser.aboutwelcome.enabled", false);
user_pref("datareporting.policy.dataSubmissionPolicyBypassNotification", true);
user_pref("toolkit.telemetry.reportingpolicy.firstRun", false);
user_pref("browser.startup.homepage_override.mstone", "ignore");'

cols=$(wc -w <<< "$ORDER"); i=0
for h in $ORDER; do
  url="http://$SERVER:$PORT/?sync&wall=$cols&screen=$i&bezel=$BEZEL"
  [ "$h" = "$LEAD" ] && url="$url&lead&ui=1"
  [ -n "$EXTRA" ] && url="$url&$EXTRA"
  if on "$h" "mkdir -p $PROFILE && printf '%s\n' '$PREFS' > $PROFILE/user.js && (DISPLAY=:0 nohup firefox --no-remote --profile $PROFILE --kiosk '$url' >/dev/null 2>&1 &)" 2>/dev/null; then
    echo "screen $i: $h"
  else
    echo "screen $i: $h is not answering, skipped"
  fi
  i=$((i + 1))
done
echo "server log: ~/.earthdna-wall-server.log    stop with: $0 stop"
