#!/bin/bash
# Post-edit syntax check for Python files

INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty')

# Only check Python files
if [[ "$FILE_PATH" != *.py ]]; then
  exit 0
fi

# Run py_compile
OUTPUT=$(python -m py_compile "$FILE_PATH" 2>&1)
if [ $? -eq 0 ]; then
  exit 0
else
  echo "Syntax error in $FILE_PATH:" >&2
  echo "$OUTPUT" >&2
  exit 2
fi
