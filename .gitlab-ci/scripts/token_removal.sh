#!/bin/bash

declare -a REPLACEMENTS

# Value pattern replacements
REPLACEMENTS+=(-e 's/glcbt-.\{23\}/[MASKED]/g')     # GitLab ... Token
REPLACEMENTS+=(-e 's/ghp_.\{36\}/[MASKED]/g')       # GitHub Personal access Token

# Variable value replacements
if [ ! -z ${SITE_ID_TOKEN+x} ]; then
	REPLACEMENTS+=(-e "s/$SITE_ID_TOKEN/[MASKED]/g")
fi
if [ ! -z ${CI_RUNNER_SHORT_TOKEN+x} ]; then
	REPLACEMENTS+=(-e 's/$CI_RUNNER_SHORT_TOKEN/[MASKED]/g')
fi
if [ ! -z ${CI_DEPENDENCY_PROXY_PASSWORD+x} ]; then
	REPLACEMENTS+=(-e 's/$CI_DEPENDENCY_PROXY_PASSWORD/[MASKED]/g')
fi

# Find all files and MASK SECRETS IN-PLACE
find "${1:-.}" -type f -exec sed -i ${REPLACEMENTS[@]} '{}' \;
