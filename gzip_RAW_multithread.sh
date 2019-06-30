#!/bin/bash
# zip .raw folder from any raw data format
# Used for proteome data submission to PRIDE database (ProteomeXchange)
# pigz is used here for parallel zipping - http://zlib.net/pigz/

echo 'Please give a output dir:'
read targetDir

targetDir="$(echo -e "${targetDir}" | sed -e 's/[[:space:]]*$//')"

if [[ "$targetDir" = *"/" ]] ; then
    targetDir=${targetDir%"/"}
    # remove trailing /
fi

for RawDir in */ ; do
    if [[ "$RawDir" = *".raw/" ]] ; then
        RawDir=${RawDir%"/"}
        # remove trailing /
        targetGzip="${targetDir}/${RawDir}.tgz"
        if [ -e $targetGzip ] ; then
            echo "${RawDir}.gzip exists"
        else
            echo "Zipping: ${RawDir}"
            tar cfv - ${RawDir} | pigz > ${targetGzip}
        fi
    fi
done
