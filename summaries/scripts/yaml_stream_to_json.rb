#!/usr/bin/env ruby
# Convert a YAML multi-document stream into a pretty JSON array.
# Robust to per-document YAML syntax issues: on parse error, embeds raw text.
# Usage: scripts/yaml_stream_to_json.rb INPUT.yaml [OUTPUT.json]

require 'yaml'
require 'json'

input = ARGV[0] || '-'
output = ARGV[1]

content = if input == '-' || input.nil?
  STDIN.read
else
  File.read(input)
end

# Split into docs on lines that are exactly '---' (YAML doc separator)
docs_raw = []
current = []
content.each_line do |line|
  if line.strip == '---'
    unless current.empty?
      docs_raw << current.join
      current = []
    end
  else
    current << line
  end
end
docs_raw << current.join unless current.empty?

docs_json = []
docs_raw.each do |doc_s|
  next if doc_s.strip.empty?
  begin
    obj = YAML.load(doc_s)
    docs_json << obj
  rescue Psych::SyntaxError => e
    fallback = { 'raw' => doc_s, 'parse_error' => e.message }
    # try to capture a doc_id if present
    if (m = doc_s.lines.find { |l| l =~ /^\s*doc_id:\s*(.+)\s*$/ })
      fallback['doc_id'] = m.strip.split(':', 2)[1].strip
    end
    docs_json << fallback
  end
end

json = JSON.pretty_generate(docs_json)

if output
  File.open(output, 'w') { |f| f.write(json) }
else
  puts json
end
