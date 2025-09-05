#!/usr/bin/env bash
set -euo pipefail

# Lint all source YAML files by loading them with Ruby/Psych.
# Exits non-zero if any file fails to parse.

shopt -s nullglob
cd "$(dirname "$0")/.."  # repo summaries dir root

files=( $(ls *.yaml | sort) )

# Exclude generated artifacts if present
files=("${files[@]}")
filtered=()
for f in "${files[@]}"; do
  [[ "$f" == "combined.yaml" ]] && continue
  filtered+=("$f")
done

pass=0
fail=0
echo "Linting ${#filtered[@]} YAML files..."
for f in "${filtered[@]}"; do
  if ruby -ryaml -e 'YAML.load_file(ARGV[0])' "$f" 2>/dev/null; then
    printf "[OK]   %s\n" "$f"
    ((pass++))
  else
    printf "[FAIL] %s\n" "$f"
    ruby -ryaml -e 'YAML.load_file(ARGV[0])' "$f" 2>&1 | sed 's/^/       /'
    ((fail++))
  fi
done

echo "Summary: $pass ok, $fail fail"
if (( fail > 0 )); then
  exit 1
fi

