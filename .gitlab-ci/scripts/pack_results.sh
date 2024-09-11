#!/bin/bash
#
# Pack each job in $OUTPATH/000000 into a UUID-named tgz file.
#
# If a variable
#
set -euo pipefail
set +x

if [ -x "${SCRIPTS_DIR:-.}/token_removal.sh" ]; then
	echo "Removing sensitive data from $OUTPATH"
	"${SCRIPTS_DIR:-.}/token_removal.sh" "$OUTPATH"
	echo "cleaning done. – Remember to always check carefully if any secrets are forgotten before publishing data packets!"
fi
cd "$OUTPATH/000000";

for job in 0*; do
    ID=$(uuidgen);
    if tar -czf "${ID}.tgz" "$job" *.xml; then
        echo "$job stored as ${ID}.tgz";
    else
        echo "$job packing failed";
    fi

    # export result blobs to datadump if available.
    if [ -z ${CB_DATADUMP+x} ]; then  # if CB_DATADUMP is not defined continue with next file.
        continue;
    fi
    if [ -d "${CB_DATADUMP}" -a -w "${CB_DATADUMP}" ] ; then
        if [ -e "${CB_DATADUMP}/${ID}.tgz" ]; then
            echo "ERROR: UUID COLLISION! This is highly unlikely. Probably something else is horribly wrong! PLEASE CHECK CAREFULLY!";
            exit 1;
        fi
        cp -v "${ID}.tgz" "${CB_DATADUMP}/";  # choose 'mv' instead of 'cp' if the results should not be added to the artifact.
    else
        echo "CB_DATADUMP=\"$CB_DATADUMP\" is not a usable directory:"
        ls -lisad "$CB_DATADUMP"
    fi;
done;

if [ ! -d "$CB_DATADUMP" ] ; then
    echo "WARNING: ${CB_DATADUMP} is not a directory! Results are only available in this ephemeral artifact!";
fi;

