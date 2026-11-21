#!/bin/bash
# Put a built globe on a Liquid Galaxy wall: the sync server (lg_server.py) on this machine, and a full-screen
# Firefox on every screen showing its slice of the picture.
#
#   lg_launch.sh GLOBE_DIR       start (restarts anything this script started before)
#   lg_launch.sh stop            close the browsers and the server
#
# Settings (environment):
#   LG_ORDER   screens left to right, as ssh host names       default "lg6 lg7 lg8 lg1 lg2 lg3 lg4 lg5"
#   LG_LEAD    the screen that leads at start, shows the controls and takes the 3D mouse    default: this host
#   LG_SERVER  address the screens reach this machine on, or "tunnel" (default): each screen reaches the server
#              through an ssh forward from this machine, for when a firewall here blocks the port
#   LG_PORT    default 8765          LG_BEZEL   css pixels hidden behind each bezel, default 0
#   LG_EXTRA   more query parameters for every screen, e.g. "dome=30&navspeed=1.5"
#   LG_EXHIBIT the full-screen videos and pictures some screens open at login (~/.config/autostart: vlc-video.desktop,
#              fullscreen-image.desktop) sit on top of the browser, so start closes them and stop reopens them
#              from those same entries; "keep" leaves them alone
#
# Screens that are switched off are skipped; the others still show their own slice.
# Google Earth keeps running underneath; the browsers cover it.
set -u
here=$(cd "$(dirname "$0")" && pwd)
ORDER=${LG_ORDER:-"lg6 lg7 lg8 lg1 lg2 lg3 lg4 lg5"}
LEAD=${LG_LEAD:-$(hostname -s)}
SERVER=${LG_SERVER:-tunnel}
PORT=${LG_PORT:-8765}
BEZEL=${LG_BEZEL:-0}
EXTRA=${LG_EXTRA:-}
PROFILE='$HOME/.earthdna-wall'  # expanded on each screen's machine
SSH="ssh -o BatchMode=yes -o ConnectTimeout=4"

# run a command on a screen's machine (locally for this one)
on() { local host=$1; shift; if [ "$host" = "$(hostname -s)" ]; then bash -c "$*"; else $SSH "$host" "$*"; fi; }

# only Firefox started with this profile; the bracket keeps pkill from matching the shell that runs it
STOP_BROWSER="pkill -f '[f]irefox .*--profile .*earthdna-wall' ; true"

EXHIBIT=${LG_EXHIBIT:-}
EXHIBIT_APPS="vlc-video.desktop fullscreen-image.desktop"
EXHIBIT_STOP="pkill -x vlc ; pkill -x lximage-qt ; true"
# reopen each exhibit entry that has nothing running for it
EXHIBIT_START="cd ~/.config/autostart 2>/dev/null || exit 0
for f in $EXHIBIT_APPS; do
  [ -f \$f ] || continue
  cmd=\$(sed -n 's/^Exec=//p' \$f | head -1); bin=\$(basename \$(echo \$cmd | cut -d' ' -f1))
  pgrep -x \$bin >/dev/null || (DISPLAY=:0 nohup sh -c \"\$cmd\" >/dev/null 2>&1 &)
done"

stop_all() {
  for h in $ORDER; do on "$h" "$STOP_BROWSER" 2>/dev/null & done; wait
  pkill -f "[s]sh .*-R $PORT:127.0.0.1:$PORT" ; true
  pkill -f "[l]g_server.py.*--port $PORT" ; true
}

if [ "${1:-}" = stop ]; then
  stop_all
  [ "$EXHIBIT" = keep ] || for h in $ORDER; do on "$h" "$EXHIBIT_START" 2>/dev/null & done; wait
  echo "stopped"; exit 0
fi
[ $# -ge 1 ] && [ -f "$1/globe.json" ] || { echo "usage: $0 GLOBE_DIR | stop   (GLOBE_DIR holds globe.json)" >&2; exit 2; }
globe=$(cd "$1" && pwd)

stop_all
[ "$EXHIBIT" = keep ] || for h in $ORDER; do on "$h" "$EXHIBIT_STOP" 2>/dev/null & done; wait
sleep 1
nohup python3 "$here/lg_server.py" "$globe" --port "$PORT" --viewer "$here/globe.html" > "$HOME/.earthdna-wall-server.log" 2>&1 &
up() { python3 -c "import urllib.request as u; u.urlopen('http://127.0.0.1:$PORT/globe.json', timeout=2).read()" 2>/dev/null; }
for _ in $(seq 20); do up && break; sleep 0.5; done
up || { echo "the server did not start; see ~/.earthdna-wall-server.log" >&2; exit 1; }

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
  host=$SERVER
  if [ "$SERVER" = tunnel ]; then
    host=127.0.0.1
    if [ "$h" != "$(hostname -s)" ]; then
      $SSH -f -N -o ExitOnForwardFailure=yes -o ServerAliveInterval=15 -R "$PORT:127.0.0.1:$PORT" "$h" 2>/dev/null \
        || { echo "screen $i: $h is not answering, skipped"; i=$((i + 1)); continue; }
    fi
  fi
  url="http://$host:$PORT/?sync&wall=$cols&screen=$i&bezel=$BEZEL"
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
