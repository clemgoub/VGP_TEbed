#!/usr/bin/env bash
# Run a local UCSC Genome Browser (official Docker image) with this repo's hub
# mounted into its docroot, and print the URL that loads the hub -- the real
# rendering test that hubCheck cannot do.
#
# Requires: Docker Desktop running. First build downloads ~1-2 GB (GBiC
# installer inside the official Dockerfile) and takes ~10-15 min; subsequent
# runs start in seconds.
#
#   ./scripts/local_browser_test.sh          # build if needed, start, print URL
#   ./scripts/local_browser_test.sh stop     # stop and remove the container
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE=vgp-ucsc-browser
CONTAINER=vgp-browser
PORT=8080

if [ "${1:-}" = "stop" ]; then
    docker rm -f $CONTAINER 2>/dev/null && echo "stopped $CONTAINER" || echo "not running"
    exit 0
fi

if ! docker info >/dev/null 2>&1; then
    echo "Docker daemon not reachable -- start Docker Desktop first." >&2
    exit 1
fi

if ! docker image inspect $IMAGE >/dev/null 2>&1; then
    echo "building $IMAGE from the official UCSC Dockerfile (one-time, ~10-15 min)..."
    TMP=$(mktemp -d)
    curl -sSf -o "$TMP/Dockerfile" \
        https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/product/installer/docker/Dockerfile
    docker build -t $IMAGE "$TMP"
    rm -rf "$TMP"
fi

docker rm -f $CONTAINER 2>/dev/null || true
# Mount the hub read-only into Apache's docroot. GBiC's docroot is
# /usr/local/apache/htdocs; the hub becomes http://localhost:8080/vgphub/.
docker run -d --name $CONTAINER -p $PORT:80 \
    -v "$REPO/hub":/usr/local/apache/htdocs/vgphub:ro \
    $IMAGE >/dev/null

echo "waiting for the browser to come up..."
for i in $(seq 1 60); do
    if curl -sf -o /dev/null "http://localhost:$PORT/cgi-bin/hgGateway"; then break; fi
    sleep 2
done
curl -sf -o /dev/null "http://localhost:$PORT/cgi-bin/hgGateway" \
    || { echo "browser did not come up; docker logs $CONTAINER" >&2; exit 1; }

HUB="http://localhost:$PORT/vgphub/hub.txt"
echo
echo "hub served inside the container:  $HUB"
echo
echo "open this URL to load the hub at the corroborated-agreement demo locus:"
echo
echo "  http://localhost:$PORT/cgi-bin/hgTracks?hubUrl=$HUB&genome=GCA_951799975.1&position=OX637609.1:15038664-15045525"
echo
echo "sanity checks the browser adds over hubCheck: track rendering, mouseovers,"
echo "the trackDb hierarchy in the track-config page, hgHubConnect validation."
echo "stop with: ./scripts/local_browser_test.sh stop"