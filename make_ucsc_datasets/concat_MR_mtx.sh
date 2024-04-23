#!/usr/bin/bash
# Merge all of the temporary mtx files
pconf -t DEVTRAJ -s multiRegion

UCSCDIR=$SDBDIR/ucsc_datasets/
ADDIR=${UCSCDIR}/ad-multi-region/

cd $ADDIR
if [[ ! -s ${ADDIR}/matrix.mtx.gz ]]; then
    cat ${ADDIR}/tmp_matrix_header.mtx > ${ADDIR}/matrix.mtx
    while read -r file; do
        echo "-- $file"
        cat ${file} >> ${ADDIR}/matrix.mtx
    done < <( ls tmp_matrix_0*mtx | sort -n )
    gzip ${ADDIR}/matrix.mtx

    gzip -t ${ADDIR}/matrix.mtx.gz

    rm ${ADDIR}/tmp_matrix_0*mtx
fi

