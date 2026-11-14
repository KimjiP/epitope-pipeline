#!/usr/bin/env bash
# Replay every locked example and check the pipeline still produces what it did.
#
#   ./reproduce.sh            replay all examples
#   ./reproduce.sh her2       replay one
#
# Pass 1 replays each example from its inputs.lock.json. It downloads the locked
# structures once (about 110 MB, into .cache/), checks every file against its
# recorded SHA-256, and compares the result with the committed summary.json.
#
# Pass 2 replays the same locks with the network blocked. A lock carries a
# snapshot of everything else the run reads, so the numbers must not change.
# Skip it with EPITOPE_OFFLINE_CHECK=0.
#
# Any difference is printed and the script exits 1.
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
CACHE="${EPITOPE_CACHE:-$HERE/.cache}"
WORK="${EPITOPE_WORK:-$HERE/.replay}"
TARGETS=("$@")
if [ ${#TARGETS[@]} -eq 0 ]; then
  TARGETS=()
  for d in "$HERE"/examples/*/; do
    [ -f "$d/inputs.lock.json" ] && TARGETS+=("$(basename "$d")")
  done
fi

mkdir -p "$CACHE" "$WORK"
fail=0

# compare WANT GOT : print one line per headline number, return 1 on any difference
compare() {
  python3 -I - "$1" "$2" <<'PY'
import json, sys
want, got = (json.load(open(p)) for p in sys.argv[1:3])
keys = ['structures_used', 'complexes', 'units', 'distinct_binders',
        'reference', 'reference_offset', 'most_targeted_patch']
bad = 0
for k in keys:
    a, b = want.get(k), got.get(k)
    if k == 'most_targeted_patch':
        a = a and (a['first_uniprot'], a['last_uniprot'], a['distinct_binders'])
        b = b and (b['first_uniprot'], b['last_uniprot'], b['distinct_binders'])
    ok = a == b
    bad += not ok
    print(f"  {'ok ' if ok else 'BAD'} {k:<20} committed {a!r}" + ('' if ok else f"   got {b!r}"))
sys.exit(1 if bad else 0)
PY
}

# ---------------------------------------------------------------- pass 1 ----
for name in "${TARGETS[@]}"; do
  lock="$HERE/examples/$name/inputs.lock.json"
  want="$HERE/examples/$name/summary.json"
  if [ ! -f "$lock" ]; then
    echo "no lock for $name"; fail=1; continue
  fi
  label=$(python3 -I -c "import json,sys;print(json.load(open(sys.argv[1]))['target'])" "$lock")
  acc=$(python3 -I -c "import json,sys;print(json.load(open(sys.argv[1]))['uniprot'])" "$lock")

  echo
  echo "=== $label ($acc) ==="
  out="$WORK/$name"
  rm -rf "${out:?}"
  if ! python3 -I "$HERE/pipeline/run_pipeline.py" \
        --from-lock "$lock" --out "$out" --cache "$CACHE" > "$WORK/$name.log" 2>&1; then
    echo "  the run failed; last lines:"
    tail -5 "$WORK/$name.log" | sed 's/^/    /'
    fail=1; continue
  fi
  grep -E "matches the hash|re-deposited|locked reference|has no snapshot" "$WORK/$name.log" | sed 's/^ */  /'
  compare "$want" "$out/summary.json" || fail=1
done

# ---------------------------------------------------------------- pass 2 ----
if [ "${EPITOPE_OFFLINE_CHECK:-1}" = "1" ] && [ "$fail" -eq 0 ]; then
  echo
  echo "--- same replays, network blocked ---"
  DEAD="http://127.0.0.1:9"          # nothing listens here: any network call fails at once
  for name in "${TARGETS[@]}"; do
    lock="$HERE/examples/$name/inputs.lock.json"
    want="$HERE/examples/$name/summary.json"
    if ! grep -q '"snapshot"' "$lock"; then
      echo "  $name: the lock has no snapshot, skipped"; continue
    fi
    out="$WORK/$name-offline"
    rm -rf "${out:?}"
    if HTTPS_PROXY="$DEAD" HTTP_PROXY="$DEAD" https_proxy="$DEAD" http_proxy="$DEAD" \
         python3 -I "$HERE/pipeline/run_pipeline.py" \
           --from-lock "$lock" --out "$out" --cache "$CACHE" > "$WORK/$name-offline.log" 2>&1; then
      echo "  $name ran with no network"
      compare "$want" "$out/summary.json" || fail=1
    else
      echo "  $name failed with the network blocked; last lines:"
      tail -5 "$WORK/$name-offline.log" | sed 's/^/    /'
      fail=1
    fi
  done
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "every example reproduced"
else
  echo "something did not reproduce"
fi
exit "$fail"
