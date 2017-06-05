Notes:

1. Polarization can be a problem -- for instance, for a two-taxon tree although in principle
    a mutation occurring on one branch is distinguishable from the *reverse* mutation occurring on the other branch
    if ancestral frequencies are assumed known, this requires in practice a large number of mutations.

2. Suppose there are mutation rates for `C->T` and also for `CA->TA`, `CC->TC`, and `CG->TG`.
    Then the `C->T` rate effectively fills in for `CT->TT`;
    if this is the highest of the four dimer rates, then the others will be negative.
    If the `CT->TT` rate is included as well, then the model is overparameterized. 

3. Mutational signatures -- if we extend the model to allow signatures (like `selfactors`) with a single coefficient,
    this could allow multistage modeling of, say, fixing the number of signatures then alternating updating the signatures
    and their relative intensities.

4. We ought to allow negative mutation rates.

5. Also, include gaps.
