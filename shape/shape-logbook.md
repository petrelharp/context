
Setup
=====

Make generator matrices, for 5, 7, 8, and 9 (doing the length-9 took about 26 hours):
```{.sh}

makegenmat () {
    LONGWIN=$1
    MODEL=$2
    Rscript ../make-genmat.R -c ${MODEL} -w ${LONGWIN} -o genmatrices/${MODEL%.json}-genmatrix-${LONGWIN}.RData"
}

for MODEL in base-model.json cpg-model.json
do
    for LONGWIN in 5 6 7 8 9
    do
        makegenmat $LONGWIN $MODEL
    done
done

for MODEL in shape-model-MGW-no-CpG shape-model-all-variables-no-CpG
do
    for LONGWIN in 5 6 7
    do
        makegenmat $LONGWIN $MODEL
    done
done

LONGWIN=8; MODEL="shape-model-MGW-no-CpG"; echo "cd $PWD; source ~/cmb/bin/R-setup-usc.sh; Rscript ../make-genmat.R -c ${MODEL}.json -w ${LONGWIN} -o genmatrices/${MODEL}-genmatrix-${LONGWIN}.RData" | qsub -q cmb -l nodes=1:dl165:ppn=24 -l walltime=120:00:00 -l mem=48000mb -l vmem=48000mb -l pmem=2000mb
LONGWIN=8; MODEL="shape-model-all-variables-no-CpG"; echo "cd $PWD; source ~/cmb/bin/R-setup-usc.sh; Rscript ../make-genmat.R -c ${MODEL}.json -w ${LONGWIN} -o genmatrices/${MODEL}-genmatrix-${LONGWIN}.RData" | qsub -q cmb -l nodes=1:dl165:ppn=24 -l walltime=120:00:00 -l mem=48000mb -l vmem=48000mb -l pmem=2000mb
LONGWIN=9; MODEL="shape-model-MGW-no-CpG"; echo "cd $PWD; source ~/cmb/bin/R-setup-usc.sh; Rscript ../make-genmat.R -c ${MODEL}.json -w ${LONGWIN} -o genmatrices/${MODEL}-genmatrix-${LONGWIN}.RData" | qsub -q cmb -l nodes=1:dl165:ppn=24 -l walltime=120:00:00 -l mem=48000mb -l vmem=48000mb -l pmem=2000mb
LONGWIN=9; MODEL="shape-model-all-variables-no-CpG"; echo "cd $PWD; source ~/cmb/bin/R-setup-usc.sh; Rscript ../make-genmat.R -c ${MODEL}.json -w ${LONGWIN} -o genmatrices/${MODEL}-genmatrix-${LONGWIN}.RData" | qsub -q cmb -l nodes=1:dl165:ppn=24 -l walltime=120:00:00 -l mem=48000mb -l vmem=48000mb -l pmem=2000mb
```

Look at the likelihood surface(s)
=================================

```{.sh}
for MODEL in base-model.json cpg-model.json shape-model-MGW-no-CpG shape-model-all-variables-no-CpG
do
    MODELFILE="${MODEL}.json"
    for TRIPLET in "5 5 0" "7 5 1" "8 6 1" "9 5 2"
    do
        read -a trip <<< "$TRIPLET"
        OUTFILE="${MODEL}-likelihood-${trip[0]}-${trip[1]}-${trip[2]}.html"
        likelihood-surface.sh $MODELFILE $OUTFILE longwin=${trip[0]} shortwin=${trip[1]} leftwin=${trip[2]} tlen=0.02 ncounts=10000000
    done
done
```


Model fitting
=============

I will start with the short ones:
```{.sh}
SCRIPT_PREFIX="cd /home/rcf-40/pralph/panfs/context/shape; source ~/cmb/bin/R-setup-usc.sh;"
# QSUB="qsub -q cmb -l nodes=1:sl230s:ppn=16 -l walltime=24:00:00 -l mem=64000mb -l vmem=64000mb -l pmem=4000mb"
QSUB="qsub -q cmb -l nodes=1:dl165:ppn=12 -l walltime=72:00:00 -l mem=24000mb -l vmem=24000mb -l pmem=1000mb"  # half a dl165
QSUB="qsub -q cmb -l nodes=1:dl165:ppn=24 -l walltime=72:00:00 -l mem=48000mb -l vmem=48000mb -l pmem=2000mb"  # a whole dl165

for COUNTDIR in $(find RegulatoryFeature-regions-from-axt -mindepth 2 -type 'd' -name "*5-5-0")
do
    for MODEL in base-model.json cpg-model.json shape-model-MGW-no-CpG shape-model-all-variables-no-CpG
    do
        MODELFILE="${MODEL}.json"
        COUNTFILE="${COUNTDIR}/total.counts.gz"
        FITFILE="${COUNTDIR}/${MODEL}-fit.RData"
        LONGWIN=$(echo $COUNTDIR | sed -e 's/.*\([0-9]\+\)-\([0-9]\+\)-\([0-9]\+\)/\1/')
        GENMAT="genmatrices/${MODEL}-genmatrix-${LONGWIN}.RData"
        ls $COUNTFILE $MODELFILE $GENMAT && \
            (echo $SCRIPT_PREFIX; echo "Rscript ../fit-model.R -i $COUNTFILE -t 1.0 --maxit 500 -c $MODELFILE -m $GENMAT -o $FITFILE")  | $QSUB
    done
done

for COUNTDIR in $(find RegulatoryFeature-regions-from-axt -mindepth 2 -type 'd' -name "*7-5-1")
do
    for MODEL in shape-model-MGW-no-CpG shape-model-all-variables-no-CpG
    do
        MODELFILE="${MODEL}.json"
        COUNTFILE="${COUNTDIR}/total.counts.gz"
        FITFILE="${COUNTDIR}/${MODEL}-fit.RData"
        LONGWIN=$(echo $COUNTDIR | sed -e 's/.*\([0-9]\+\)-\([0-9]\+\)-\([0-9]\+\)/\1/')
        GENMAT="genmatrices/${MODEL}-genmatrix-${LONGWIN}.RData"
        ls $COUNTFILE $MODELFILE $GENMAT && \
            (echo $SCRIPT_PREFIX; echo "Rscript ../fit-model.R -i $COUNTFILE -t 1.0 --maxit 500 -c $MODELFILE -m $GENMAT -o $FITFILE")  | $QSUB
    done
done

```
