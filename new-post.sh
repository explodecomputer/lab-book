#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<-EOF
	Create a new post for the lab book.

	Usage: $0 <title> [options]

	Options:
	  -n, --notebook    Create a .ipynb instead of .qmd
	  --slug <slug>     Custom slug (default: auto-derived from title)
	  -h, --help        Show this help

	Examples:
	  $0 "Post title"
	  $0 "Post title" --notebook
	  $0 "Post title" --slug custom-slug
	EOF
    exit 0
}

TITLE=""
NOTEBOOK=false
CUSTOM_SLUG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--notebook) NOTEBOOK=true; shift ;;
        --slug) CUSTOM_SLUG="$2"; shift 2 ;;
        -h|--help) usage ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *) TITLE="$1"; shift ;;
    esac
done

if [[ -z "$TITLE" ]]; then
    echo "Error: title is required"
    usage
fi

DATE=$(date +%Y-%m-%d)

if [[ -n "$CUSTOM_SLUG" ]]; then
    SLUG="$CUSTOM_SLUG"
else
    SLUG=$(echo "$TITLE" \
        | tr '[:upper:]' '[:lower:]' \
        | sed 's/[^a-z0-9 ]//g' \
        | sed 's/  */ /g' \
        | sed 's/ /-/g' \
        | sed 's/^-//;s/-$//' \
        | cut -c1-80 \
        | sed 's/-$//')
fi

FOLDER="${DATE}-${SLUG}"
DIR="posts/${FOLDER}"

COUNTER=2
while [[ -d "$DIR" ]]; do
    DIR="posts/${FOLDER}-${COUNTER}"
    COUNTER=$((COUNTER + 1))
done

mkdir -p "$DIR"

if $NOTEBOOK; then
    python3 /dev/stdin "$TITLE" "$DATE" > "$DIR/index.ipynb" <<-'PYEOF'
import json, sys

title = sys.argv[1]
date = sys.argv[2]

nb = {
    "cells": [
        {
            "cell_type": "markdown",
            "id": "yaml-front-matter",
            "metadata": {},
            "source": [
                f"---\n",
                f"title: \"{title}\"\n",
                f"author: Gibran Hemani\n",
                f"date: {date}\n",
                f"categories: []\n",
                f"---"
            ]
        },
        {
            "cell_type": "markdown",
            "id": "background",
            "metadata": {},
            "source": [
                "## Background\n"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "id": "code-cell",
            "metadata": {},
            "outputs": [],
            "source": []
        }
    ],
    "metadata": {
        "language_info": {
            "name": "r"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 5
}

json.dump(nb, sys.stdout, indent=1)
print()
PYEOF
    echo "Created $DIR/index.ipynb"
else
    cat > "$DIR/index.qmd" <<-EOF
---
title: "$TITLE"
author: Gibran Hemani
date: "$DATE"
categories: []
---

## Background

---

\`\`\`{r}
sessionInfo()
\`\`\`
EOF
    echo "Created $DIR/index.qmd"
fi
