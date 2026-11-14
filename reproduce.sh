#!/usr/bin/env bash
# Replay every locked example and check the pipeline still produces what it did.
#
#   ./reproduce.sh            replay all examples
#   ./reproduce.sh her2       replay one
#
# Each replay uses the exact PDB entries in that example's inputs.lock.json, checks
# every file against its recorded SHA-256, and compares the result with the
# committed summary.json. Any difference is printed and the script exits 1.
#
# Downloads about 120 MB the first time, into .cache/, and reuses it after.
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

for name in "${TARGETS[@]}"; do
  lock="$HERE/examples/$name/inputs.lock.json"
  want="$HERE/examples/$name/summary.json"
  if [ ! -f "$lock" ]; then
    echo "no lock for $name"; fail=1; continue
  fi
  acc=$(python3 -I -c "import json,sys;print(json.load(open(sys.argv[1]))['uniprot'])" "$lock")
  label=$(python3 -I -c "import json,sys;print(json.load(open(sys.argv[1]))['target'])" "$lock")

  echo
  echo "=== $label ($acc) ==="
  out="$WORK/$name"
  rm -rf "$out"
  if ! python3 -I "$HERE/pipeline/run_pipeline.py" \
        --uniprot "$acc" --name "$label" --from-lock "$lock" \
        --out "$out" --cache "$CACHE" > "$WORK/$name.log" 2>&1; then
    echo "  the run failed; last lines:"
    tail -5 "$WORK/$name.log" | sed 's/^/    /'
    fail=1; continue
  fi
  grep -E "matches the hash|re-deposited|locked reference" "$WORK/$name.log" | sed 's/^ */  /'

  python3 -I - "$want" "$out/summary.json" <<'PY' || fail=1
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
done

echo
if [ "$fail" -eq 0 ]; then
  echo "every example reproduced"
else
  echo "something did not reproduce"
fi
exit "$fail"
