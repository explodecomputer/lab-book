#!/bin/sh
set -eu
report_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
stage_dir=$(mktemp -d "${TMPDIR:-/tmp}/chd-quarto.XXXXXX")
trap 'rm -rf "$stage_dir"' EXIT HUP INT TERM
cp "$report_dir/chd_followup.qmd" "$stage_dir/"
ln -s "$report_dir/offline" "$stage_dir/offline"
mkdir -p "$report_dir/chd-followup-results"
ln -s "$report_dir/chd-followup-results" "$stage_dir/chd-followup-results"
export http_proxy=http://127.0.0.1:9 https_proxy=http://127.0.0.1:9 ALL_PROXY=http://127.0.0.1:9 NO_PROXY=
quarto render "$stage_dir/chd_followup.qmd" --execute-dir "$report_dir" --no-cache
sed 's/[[:blank:]]*$//' "$stage_dir/chd_followup.html" > "$report_dir/chd_followup.html"
