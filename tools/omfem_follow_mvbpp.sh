#!/bin/bash
# Follow MVB++ and re-run the OMFEM corpus on every new version, never against a stale binary.
#
# Alf, 2026-09-20: "make sure all the designs start running in the OMFEM agent as soon as a new
# MVB++ version is available ... as many iterations over night as possible, and make sure they
# don't use a stale binary."
#
# THE STALENESS RULE, and why this script exists at all. OMFEM links MVB++/build/libmvb++.so and
# shares MAS/MKF types with it, so an OMFEM binary older than that library is not merely out of
# date -- it can read a differently-laid-out object. OMFEM already refuses to run in that state
# ("ERROR: OMFEM binary is STALE against MVB++"), which is correct and must not be worked around.
# The only safe order is therefore: MVB++ builds -> OMFEM rebuilds -> corpus runs. This script is
# that order in a loop, and it never starts a corpus it has not just rebuilt for.
#
# WHAT IT WATCHES: the MVB++ commit at HEAD *and* the mtime of libmvb++.so. Both, because a commit
# with no rebuild is not yet in the library, and a rebuild with no commit (someone iterating in the
# tree) still changes what OMFEM would link.
#
# RUN IT FROM THE OMFEM SIDE. It builds OMFEM, which is OMFEM's build directory to write; this
# script does not build MVB++ and must not. If MVB++'s library is newer than its own sources'
# last build, it waits rather than guessing -- an unbuilt MVB++ change is the one case where
# rebuilding OMFEM would bake in a half-finished state.
#
#   usage: tools/omfem_follow_mvbpp.sh [--once] [--interval SECONDS] [--corpus "<command>"]
#
set -u

MVBPP=/home/alf/OpenMagnetics/MVB++
OMFEM=/home/alf/OpenMagnetics/OMFEM
LIB=$MVBPP/build/libmvb++.so
STATE=${OMFEM_FOLLOW_STATE:-$OMFEM/.omfem_follow_state}
INTERVAL=300
ONCE=0
# The corpus command. Override with --corpus to run a subset, a different lane, or a dry run.
CORPUS=${OMFEM_CORPUS_CMD:-"python3 tools/omfem_corpus_real.py"}

while [ $# -gt 0 ]; do
  case "$1" in
    --once) ONCE=1; shift ;;
    --interval) INTERVAL=$2; shift 2 ;;
    --corpus) CORPUS=$2; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

log() { echo "[follow $(date +%H:%M:%S)] $*"; }

# The version of MVB++ that OMFEM would be linking right now: commit + library mtime.
mvbpp_version() {
  local sha lib_mtime
  sha=$(git -C "$MVBPP" rev-parse --short HEAD 2>/dev/null || echo nogit)
  lib_mtime=$(stat -c %Y "$LIB" 2>/dev/null || echo 0)
  echo "$sha:$lib_mtime"
}

# Is the library actually built from the current sources? If a source file is newer than the
# library, MVB++ is mid-edit: wait. Rebuilding OMFEM against a library that is about to change
# just spends an hour to arrive stale.
mvbpp_settled() {
  local newest
  newest=$(find "$MVBPP/src" "$MVBPP/include" -name '*.cpp' -o -name '*.h' 2>/dev/null \
           | xargs -r stat -c %Y 2>/dev/null | sort -rn | head -1)
  [ -z "$newest" ] && return 0
  local lib_mtime
  lib_mtime=$(stat -c %Y "$LIB" 2>/dev/null || echo 0)
  [ "$lib_mtime" -ge "$newest" ]
}

# Do not start a build or a corpus while the box is short of memory: several sessions share it.
wait_for_memory() {
  while [ "$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)" -lt 8 ]; do
    log "waiting for memory (need 8 GB available)"
    sleep 60
  done
}

run_once() {
  local version
  version=$(mvbpp_version)
  local last=""
  [ -f "$STATE" ] && last=$(cat "$STATE")

  if [ "$version" = "$last" ]; then
    return 1   # nothing new
  fi
  if ! mvbpp_settled; then
    log "MVB++ sources are newer than libmvb++.so -- waiting for its build to finish"
    return 1
  fi

  log "new MVB++ version $version (was ${last:-none}) -- rebuilding OMFEM"
  wait_for_memory
  if ! ( cd "$OMFEM" && CMAKE_BUILD_PARALLEL_LEVEL=3 cmake --build build -j3 ); then
    log "OMFEM build FAILED against $version -- not running the corpus, and not recording the"
    log "version, so the next tick retries rather than silently skipping this MVB++ change"
    return 1
  fi

  # Record the version ONLY after a successful rebuild: the state file means "OMFEM is built
  # against this", not "this was seen".
  echo "$version" > "$STATE"

  log "OMFEM rebuilt -- starting the corpus against MVB++ $version"
  wait_for_memory
  ( cd "$OMFEM" && eval "$CORPUS" )
  log "corpus finished for $version (exit $?)"
  return 0
}

if [ "$ONCE" = 1 ]; then
  run_once
  exit $?
fi

log "following MVB++ every ${INTERVAL}s; state in $STATE"
while true; do
  run_once
  sleep "$INTERVAL"
done
