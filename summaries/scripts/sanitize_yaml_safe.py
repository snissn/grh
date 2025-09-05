#!/usr/bin/env python3
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
            # toggle quote state if not escaped
            if i == 0 or line[i-1] != '\\':
                in_dq = not in_dq
            out.append(ch)
            i += 1
            continue
        if in_dq and ch == '\\':
            # if already escaped backslash, keep pair
            if i+1 < len(line) and line[i+1] == '\\':
                out.append('\\\\')
                i += 2
            else:
                out.append('\\\\')
                i += 1
            continue
        out.append(ch)
        i += 1
    return ''.join(out)


_key_star = re.compile(r'^(\s*)([^"\'\[{][^:]*\*[^:]*):')


def quote_keys_with_star(line: str) -> str:
    m = _key_star.match(line)
    if m:
        return f'{m.group(1)}"{m.group(2)}":' + line[m.end():]
    return line


def quote_src_in_flow(line: str) -> str:
    # Quote src: value inside {...} if value is unquoted and contains ':'
    out = []
    i = 0
    while i < len(line):
        if line[i] == '{':
            # find matching '}' on this line (naive; braces are single-line here)
            j = line.find('}', i)
            if j == -1:
                out.append(line[i:])
                break
            seg = line[i:j]
            # replace occurrences of src: <bare-with-colon>
            seg = re.sub(r'(src:\s*)([^,\"\'}][^,}]*)',
                         lambda m: m.group(1) + ('"'+m.group(2)+'"' if ':' in m.group(2) else m.group(2)),
                         seg)
            out.append(seg)
            out.append('}')
            i = j + 1
        else:
            out.append(line[i])
            i += 1
    return ''.join(out)


def quote_src_block(line: str) -> str:
    # Lines like: src: def:test-class
    m = re.match(r'^(\s*src:\s*)([^\'\"\[{].*)$', line)
    if m:
        val = m.group(2)
        if ':' in val:
            return f'{m.group(1)}"{val.strip()}"\n'
    return line


def sanitize(path: Path) -> bool:
    s = path.read_text(encoding='utf-8')
    out_lines = []
    changed = False
    for line in s.splitlines(keepends=True):
        orig = line
        line = fix_backslashes_in_dq(line)
        line = quote_keys_with_star(line)
        line = quote_src_in_flow(line)
        line = quote_src_block(line)
        if line != orig:
            changed = True
        out_lines.append(line)
    if changed:
        path.write_text(''.join(out_lines), encoding='utf-8')
    return changed


def main(argv):
    if len(argv) < 2:
        print('usage: sanitize_yaml_safe.py file1.yaml [file2.yaml ...]')
        return 2
    for fn in argv[1:]:
        p = Path(fn)
        if sanitize(p):
            print('fixed', fn)
        else:
            print('ok   ', fn)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))

