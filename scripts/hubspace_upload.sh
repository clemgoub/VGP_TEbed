#!/usr/bin/env bash
# Upload the built hub to UCSC hubSpace -- storage co-located with the public
# Genome Browser (10 GB per account), loadable there by URL without making the
# hub public: it is not listed anywhere, only people given the URL see it.
#
# One-time setup (browser, logged in to your UCSC account):
#   genome.ucsc.edu -> My Data -> Track Hubs -> Track Development -> API key
#   echo 'apiKey=PASTE_KEY_HERE' > ~/.hubtools.conf     # NO quotes around the key
#
# Then:  ./scripts/hubspace_upload.sh
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
HUBNAME=vgpRepeatConsensus

[ -f ~/.hubtools.conf ] || {
    echo "no ~/.hubtools.conf -- create the API key first (see header)" >&2; exit 1; }

# hubtools' conf parser does line.split("=",1) and never strips quotes, yet its
# own usage text shows apiKey="xxxx" -- following it sends literal quotes in
# the apiKey field and the TUS endpoint answers 500. Strip them defensively.
if grep -q 'apiKey=["'"'"']' ~/.hubtools.conf; then
    sed -i '' 's/^apiKey=["'"'"']*\([^"'"'"']*\)["'"'"']*[[:space:]]*$/apiKey=\1/' ~/.hubtools.conf
    echo "note: removed quotes around apiKey in ~/.hubtools.conf (hubtools sends them verbatim)"
fi

SIZE=$(du -sh "$REPO/hub" | cut -f1)
echo "uploading hub/ ($SIZE) as '$HUBNAME' -- 1.4G over ~10-60 min depending on uplink"
cd "$REPO"

# hubtools needs tuspy (tus resumable-upload client) and tries to pip-install
# it into whatever python3 is on PATH -- which a Homebrew python refuses
# (PEP 668 externally-managed-environment). Give it a private venv instead.
VENV="$REPO/.hubtools-venv"
if ! "$VENV/bin/python" -c 'import tusclient' 2>/dev/null; then
    echo "one-time setup: creating $VENV with tuspy..."
    python3 -m venv "$VENV"
    "$VENV/bin/pip" -q install tuspy requests
fi
"$VENV/bin/python" ./bin/hubtools up "$HUBNAME" -i hub

echo
echo "done. Connect it on the public browser:"
echo "  genome.ucsc.edu -> My Data -> Track Hubs -> Hub Upload tab -> $HUBNAME"
echo "or share the hub.txt URL shown there with the group:"
echo "  https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=<that URL>&genome=GCA_951799975.1&position=OX637609.1:15038664-15045525"