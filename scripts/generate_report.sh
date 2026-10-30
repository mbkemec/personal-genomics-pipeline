
#!/usr/bin/env bash
set -euo pipefail

# Find project root regardless of where the script is executed.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPORT_DIR="$PROJECT_ROOT/reports"

# Check dependencies.
for cmd in pandoc xelatex; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd is not installed." >&2
        exit 1
    fi
done

# Ask for the Markdown filename.
read -r -p "Markdown report filename: " REPORT_NAME

if [[ -z "$REPORT_NAME" ]]; then
    echo "ERROR: Filename cannot be empty." >&2
    exit 1
fi

# Allow filenames with or without the .md extension.
REPORT_NAME="${REPORT_NAME%.md}"

INPUT="$REPORT_DIR/$REPORT_NAME.md"
OUTPUT="$REPORT_DIR/$REPORT_NAME.pdf"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: File not found: $INPUT" >&2
    exit 1
fi

echo "Generating PDF: $REPORT_NAME.pdf"

# Run from the reports directory so relative image paths work.
cd "$REPORT_DIR"


pandoc "$INPUT" \
    --from markdown+pipe_tables+link_attributes \
    --to pdf \
    --pdf-engine=xelatex \
    --resource-path="$REPORT_DIR:$PROJECT_ROOT" \
    --standalone \
    --toc \
    --toc-depth=2 \
    -V documentclass=article \
    -V mainfont="DejaVu Sans" \
    -V geometry:margin=2.2cm \
    -V fontsize=10pt \
    -V colorlinks=true \
    -V linkcolor=black \
    -V urlcolor=blue \
    -o "$OUTPUT"


echo "PDF generated successfully: $OUTPUT"
