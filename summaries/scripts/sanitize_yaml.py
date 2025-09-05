#!/usr/bin/env python3
"""
Lightweight YAML sanitiser for our summary files.

Fixes common issues that Ruby/Psych is strict about:
- Doubles literal backslashes inside double-quoted scalars.
- Quotes bare mapping values that contain ':' (both block and flow style).
- Quotes mapping keys that contain '*' (e.g., symbols like 𝔞*), which YAML
  would otherwise treat as an alias.

Usage: scripts/sanitize_yaml.py file1.yaml [file2.yaml ...]
Modifies files in-place (UTF-8).
"""
import re
import sys
from pathlib import Path


def fix_backslashes_in_dq(line: str) -> str:
    out = []
    i = 0
    in_dq = False
    while i < len(line):
        ch = line[i]
        if ch == '"':
            # toggle double-quote state unless escaped
            # if it's escaped (preceded by a backslash), leave state unchanged
            # naive but fine for our docs (we rarely escape quotes)
            if i == 0 or line[i - 1] != '\\':
                in_dq = not in_dq
            out.append(ch)
            i += 1
            continue
        if in_dq and ch == '\\':
            # ensure a literal backslash is represented as \\
            if i + 1 < len(line) and line[i + 1] == '\\':
                # already escaped; keep as is
                out.append('\\\\')
                i += 2
            else:
                out.append('\\\\')
                i += 1
            continue
        out.append(ch)
        i += 1
    return ''.join(out)


_key_with_star_re = re.compile(r'^(\s*)([^"\'\[{][^:]*\*[^:]*):')


def quote_keys_with_star(line: str) -> str:
    # Quote keys containing '*' to avoid alias parsing
    m = _key_with_star_re.match(line)
    if m:
        indent, key = m.groups()
        return f'{indent}"{key}":' + line[m.end():]
    return line


_flow_kv_re = re.compile(r'(\{[^}]*\})')
_flow_value_needs_quote = re.compile(r'(:\s*)([^,}\"\']*:[^,}\"\']*)')


def quote_colon_values_in_flow(line: str) -> str:
    # Inside {...}, quote values after k: that contain ':' and are unquoted
    def _fix_segment(seg: str) -> str:
        return _flow_value_needs_quote.sub(lambda m: m.group(1) + '"' + m.group(2) + '"', seg)

    parts = []
    last = 0
    for m in _flow_kv_re.finditer(line):
        # append text before {...}
        parts.append(line[last:m.start()])
        seg = m.group(1)
        parts.append(_fix_segment(seg))
        last = m.end()
    parts.append(line[last:])
    return ''.join(parts)


_block_colon_value_re = re.compile(r'^(\s*[^:#\n]+:\s*)([^\'\"\[{][^#\n]*)$')


def quote_colon_values_in_block(line: str) -> str:
    # For lines like: src: def:test-class
    # If the value contains a ':' and is unquoted/unstyled, wrap in quotes.
    m = _block_colon_value_re.match(line)
    if m:
        prefix, val = m.groups()
        if ':' in val:
            return f'{prefix}"{val.strip()}"\n'
    return line


def sanitize_file(path: Path) -> bool:
    text = path.read_text(encoding='utf-8')
    changed = False
    new_lines = []
    for line in text.splitlines(keepends=True):
        orig = line
        line = fix_backslashes_in_dq(line)
        line = quote_keys_with_star(line)
        line = quote_colon_values_in_flow(line)
        line = quote_colon_values_in_block(line)
        if line != orig:
            changed = True
        new_lines.append(line)
    if changed:
        path.write_text(''.join(new_lines), encoding='utf-8')
    return changed


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    rc = 0
    for arg in argv[1:]:
        p = Path(arg)
        if not p.exists():
            print(f"warn: {p} missing")
            rc = 1
            continue
        if sanitize_file(p):
            print(f"fixed: {p}")
        else:
            print(f"ok:    {p}")
    return rc


if __name__ == '__main__':
    sys.exit(main(sys.argv))

