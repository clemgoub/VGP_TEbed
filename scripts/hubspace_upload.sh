#!/usr/bin/env bash
# Upload the built hub to UCSC hubSpace -- storage co-located with the public
# Genome Browser (10 GB per account), loadable there by URL without making the
# hub public: it is not listed anywhere, only people given the URL see it.
#
# One-time setup (browser, logged in to your UCSC account):
#   genome.ucsc.edu -> My Data -> Track Hubs -> Track Development -> API key
#   echo 'apiKey="PASTE_KEY_HERE"' > ~/.hubtools.conf
#
# Then:  ./scripts/hubspace_upload.sh
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
HUBNAME=vgpRepeatConsensus

[ -f ~/.hubtools.conf ] || {
    echo "no ~/.hubtools.conf -- create the API key first (see header)" >&2; exit 1; }

SIZE=$(du -sh "$REPO/hub" | cut -f1)
echo "uploading hub/ ($SIZE) as '$HUBNAME' -- 1.4G over ~10-60 min depending on uplink"
cd "$REPO"
./bin/hubtools up "$HUBNAME" -i hub

echo
echo "done. Connect it on the public browser:"
echo "  genome.ucsc.edu -> My Data -> Track Hubs -> Hub Upload tab -> $HUBNAME"
echo "or share the hub.txt URL shown there with the group:"
echo "  https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=<that URL>&genome=GCA_951799975.1&position=OX637609.1:15038664-15045525"