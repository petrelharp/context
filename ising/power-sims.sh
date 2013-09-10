TLEN="-t 0.2"
BETA="-b 1"
GAMMA="-m .5"
SEQLEN="1e4"
WIN="-w 3"
LWIN="-l 3"

# common reference
SIMDIR="reference-sims"
for k in $(seq 10)
do
    Rscript sim-ising.R -d $SIMDIR $TLEN $BETA $GAMMA $SEQLEN
done

SIMDIR="seqlen-sims"
mkdir -p $SIMDIR
for slen in 1e2 1e3 1e5
do
    for k in $(seq 10)
    do
        Rscript sim-ising.R -d $SIMDIR $TLEN $BETA $GAMMA -s $slen
    done
done

SIMDIR="tlen-sims"
mkdir -p $SIMDIR
for tlen in .05 .1 .4
do
    for k in $(seq 10)
    do
        Rscript sim-ising.R -d $SIMDIR -t $tlen $BETA $GAMMA $SEQLEN 
    done
done

SIMDIR="gamma-sims"
mkdir -p $SIMDIR
for gamma in .05 .1 .4
do
    for k in $(seq 10)
    do
        Rscript sim-ising.R -d $SIMDIR $TLEN $BETA -m $gamma $SEQLEN 
    done
done
