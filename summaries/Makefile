# Combine all YAML summaries into a single AI-friendly YAML stream.
#
# Outputs
# - combined.yaml: YAML stream
#   - Doc 1: collection_meta + table_of_contents
#   - Subsequent docs: one per input file, as a mapping with
#     doc_id, source, ingested_at, then the file's original content
#
# Usage
#   make            # builds combined.yaml
#   make clean      # removes generated files

SHELL := /bin/bash

# Input YAML sources (exclude generated combined.yaml)
SOURCES := $(filter-out combined.yaml,$(sort $(wildcard *.yaml)))
COMBINED := combined.yaml
COMBINED_JSON := combined.json
DATE := $(shell date -u +"%Y-%m-%dT%H:%M:%SZ")


# Phony targets
.PHONY: all clean show show-json pack json yaml lint

# Default: build both YAML stream and JSON array
all: $(COMBINED_JSON)

# Convenience targets
yaml: $(COMBINED)
json: $(COMBINED_JSON)

lint:
	@bash scripts/lint_yaml.sh

$(COMBINED): $(SOURCES) scripts/combine_yaml.sh
	@bash scripts/combine_yaml.sh "$(COMBINED)" $(SOURCES)

$(COMBINED_JSON): $(COMBINED) scripts/yaml_stream_to_json.rb
	@ruby scripts/yaml_stream_to_json.rb "$(COMBINED)" "$(COMBINED_JSON)"

show: $(COMBINED)
	@echo "Built $(COMBINED) with $$(grep -c '^---' $(COMBINED)) documents."; \
	head -n 30 $(COMBINED)

show-json: $(COMBINED_JSON)
	@echo "Built $(COMBINED_JSON) (pretty JSON array of documents)."; \
	jq -r 'length as $n | "docs: \($n)"' $(COMBINED_JSON) 2>/dev/null || head -n 30 $(COMBINED_JSON)

pack: $(COMBINED).gz $(COMBINED_JSON).gz

$(COMBINED).gz: $(COMBINED)
	@gzip -c "$<" > "$@"

$(COMBINED_JSON).gz: $(COMBINED_JSON)
	@gzip -c "$<" > "$@"

clean:
	rm -f $(COMBINED) $(COMBINED_JSON) $(COMBINED).gz $(COMBINED_JSON).gz
