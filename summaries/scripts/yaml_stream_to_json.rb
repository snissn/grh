#!/usr/bin/env ruby
# Convert a YAML multi-document stream into a pretty JSON array.
# Strict mode: any YAML parse error causes a non-zero exit (Make fails).
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
parse_errors = []
docs_raw.each_with_index do |doc_s, idx|
  next if doc_s.strip.empty?
  begin
    obj = YAML.load(doc_s)
    docs_json << obj
  rescue Psych::SyntaxError => e
    # capture a doc_id if present for better diagnostics
    doc_id = nil
    if (m = doc_s.lines.find { |l| l =~ /^\s*doc_id:\s*(.+)\s*$/ })
      doc_id = m.strip.split(':', 2)[1].strip
    end
    parse_errors << { index: idx, doc_id: doc_id, message: e.message }
  end
end

if parse_errors.any?
  STDERR.puts "YAML parse failed for #{parse_errors.length} document(s):"
  parse_errors.each do |err|
    where = "doc ##{err[:index]}"
    where += " (doc_id: #{err[:doc_id]})" if err[:doc_id]
    STDERR.puts " - #{where}: #{err[:message]}"
  end
  exit 1
end

json = JSON.pretty_generate(docs_json)

if output
  File.open(output, 'w') { |f| f.write(json) }
else
  puts json
end
