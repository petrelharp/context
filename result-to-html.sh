#!/bin/sh

for f in $@; do
    base=$(basename $f .RData)
    test $f = $base || { # If it indeed has an .RData extension then process.
        Rscript $(dirname $0)/RData-to-md.R -o $base.md $f
        pandoc \
            -f markdown-yaml_metadata_block \
            -o $base.html \
            -c http://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/css/bootstrap.min.css \
            $base.md
    }
done
