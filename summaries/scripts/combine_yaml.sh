#!/usr/bin/env bash
set -euo pipefail

# Combine YAML files into a single YAML stream with a header doc.
# Usage: scripts/combine_yaml.sh OUTPUT.yaml file1.yaml file2.yaml ...

out=${1:-combined.yaml}
shift || true

shopt -s nullglob
files=("$@")
# Ensure we never include the output file itself if passed
out_base=$(basename "$out")
filtered=()
for f in "${files[@]}"; do
  [[ $(basename "$f") == "$out_base" ]] && continue
  filtered+=("$f")
done
files=("${filtered[@]}")
if [[ ${#files[@]} -eq 0 ]]; then
  # fallback to local directory
  mapfile -t files < <(printf '%s\n' *.yaml | sort | grep -v "^$(printf '%q' "$out_base")$")
fi

ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
tmp="$out.tmp"
{
  echo '---'
  echo 'collection_meta:'
  echo "  generated_at: $ts"
  echo "  source_dir: $(pwd)"
  echo "  file_count: ${#files[@]}"
  echo 'table_of_contents:'
  for f in "${files[@]}"; do
    id=${f%.yaml}
    # Extract a few helpful meta fields if present (best-effort)
    title=$(sed -n 's/^[[:space:]]*title_long:[[:space:]]*\(.*\)$/\1/p' "$f" | head -n1)
    if [[ -z "$title" ]]; then
      title=$(sed -n 's/^[[:space:]]*title_short:[[:space:]]*\(.*\)$/\1/p' "$f" | head -n1)
    fi
    authors=$(sed -n 's/^[[:space:]]*authors:[[:space:]]*\(.*\)$/\1/p' "$f" | head -n1)
    year=$(sed -n 's/^[[:space:]]*year:[[:space:]]*\(.*\)$/\1/p' "$f" | head -n1)
    area=$(sed -n 's/^[[:space:]]*area:[[:space:]]*\(.*\)$/\1/p' "$f" | head -n1)
    echo "  - doc_id: $id"
    echo "    file: $f"
    [[ -n "$title" ]]   && echo "    title: $title"
    [[ -n "$authors" ]] && echo "    authors: $authors"
    [[ -n "$year" ]]    && echo "    year: $year"
    [[ -n "$area" ]]    && echo "    area: $area"
  done

  for f in "${files[@]}"; do
    id=${f%.yaml}
    echo '---'
    echo "doc_id: $id"
    echo "source: $f"
    echo "ingested_at: $ts"
    cat "$f"
  done
} > "$tmp"

mv "$tmp" "$out"
