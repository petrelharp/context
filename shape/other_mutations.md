# list of processes

Generally: interaction between of *source of damage* and *repair method*.

Very good review: Roberts & Gordenin 2014

    - Another mutation signature, [AT -> AC], detected in HV-C positive hepatocellular carcinomas was suggested to originate from error prone synthesis by one of the TLS polymerases based on overexpression of POL zeta and POL iota in this type of cancer

    - pol iota does A->G 3:1 more often than A->A; subsequently extended by zeta or eta
        * and "GT/TT tandem mismatches occur slightly more often than the AA/TT correct matches,"
            * leading to AA -> GT
        * "increased probability of fixing the favored G:T mispair as a mutation when T is the 5′-template base."
            * leading to AA -> AG

Review of CpG by Pfeifer 2006 

- spontaneous deamination (C->U) at CpG

- "When adjacent to another pyrimidine, 5-methylcytosine preferentially undergoes sunlight-induced pyrimidine dimer formation." and "The vast majority of base changes in skin lesions are C to T or CC to TT mutations at dipyrimidine sequences."
    * TCG -> TTG or CCG -> CTG
    * but this is maybe combination of dipyramidine C->T and CpG

- "Certain polycyclic aromatic hydrocarbons form guanine adducts and induce G to T transversion mutations with high selectivity at mCpG sequences."
    * CG -> CT

From https://www.nature.com/articles/ncomms9866 which did NMF on some cancer:

- AID likes to deaminate Cs in RC so AC to AT or GC to GT

- deamination (C->U) due to AID resolved by pol eta leads to WA to WC (W=A or T)

- AID followed by normal replication: C to T or G at WRCY (R is purine, A or G, Y is pyramidine, T or C)
    * from https://www.ncbi.nlm.nih.gov/pubmed/17328676
    * WRCY = AACT / AACC / AGCT / AGCC / TACT / TACC / TGCT / TGCC

in Alexandrov 2013 Nature, massive NMF on lots of cancers:

- "Signature 2" - probably "APOBEC family of cytidine deaminases, which convert cytidine to uracil" causing TC to TT or TG
    * actually maybe TCW to TTW or TGW (Roberts & Gordenin review)

- UV-induced pyramidine dimers: mostly CC to TT

- suggested that pol eta responsible for AT to AG and TT to TG

Rogovin et al looked for Ig hotspots and correlated these to in vitro noisy pol spectra:

- pol eta: WA -> WG (gives relative proportions: 49% ->WG, 32% ->WT, 19%->WC

Other notes:

- strand bias caused by "transcription-coupled nucleotide excision repair (NER) that operates predominantly on the transcribed strand of genes and is recruited by RNA polymerase II when it encounters bulky DNA helix-distorting lesions"
- microsattelite instability (many small indels) caused by errors in mismatch repair
- kataegis is regions with hella mutations (from "thunder")
- " It is possible that both result from cytidine to uracil deamination by an APOBEC family member, but that the different signatures are sequelae of different repair mechanisms following the deamination step. C>T transitions may simply result from DNA replication across uracil. However, if uracil is excised by uracil-DNA glycosylase (UNG) as part of base excision repair (BER), an abasic site is generated (Wilson and Bohr, 2007). The partiality for C>G transversions in Signature E may reflect preferential insertion of cytosine opposite such an UNG-mediated abasic site. The propensity to introduce cytosine opposite an abasic site is characteristic of REV1 translesion polymerase"
- "uracil DNA glycosylase (UNG), a base-excision repair enzyme that removes uracil from DNA."
- "The mainstream model for SHM holds that base excision repair (BER) broadens the range of mutations at C:G basepairs, and mismatch repair (MMR) extends the mutation process to A:T basepairs"

# arbeithuber2015crossovers: https://www.ncbi.nlm.nih.gov/pubmed/25646453

- "Observed de novo mutations changed strong (S) CG into weak (W) TA base pairs and they all occurred mainly at CpG sites. As it is shown that GC base pairs are preferentially transferred during crossover, the authors suggested that GC-biased gene conversion (gBGC) is the dominant force shaping the nucleotide composition at hotspots and potentially in other recombination products, which might explain the high GC content associated with recombination. It is possible that gBGC is an adaptation to reduce the opposing mutational load of recombination, knowing that mutation favors weak over strong nucleotides. "


# matthieson and reich:

- possibly constrained by choice of 3-1-1 Tmers
- signature 1: "TCC>T, ACC>T, CCC>T and TCT>T (possibly also including CCG>T)"
- "Signature 3 is characterized by GT>GG mutations, particularly GTG>GGG." but is probably an artifact



# Additional signatures after adding more CpG mutations:

These have z-scores of 20--40

high!   .... low!!   also high
GG  AG  ....    CC  CT
GT  AT  ....    AC  AT
TA  CA  ....           TA  TG
GC  AC  .... GC  GT
AC  GC  ....           GT  GC
GT  CT  ....           AC  AG
GA  AA  ....    TC  TT
GT  TT  ....    AC  AA
--> doubles!
GT  AC  ....           AC  GT
TT  GA  ....           AA  TC

These have z-scores of 10--20

high!    ... also high 
CTAT  TG ... ATAG  CA   
CTGT  TC ... ACAG  GA
CTAC  TG ... GTAG  CA
ATTG  CT ... CAAT  AG
--> double!            z
GTTG  GA ... CAAC  TC  7
--> but single overrepresented
GTTG  CT ... CAAC  AG  5
GTTG  CT ... CAAC  AG  5
GTTG  GT ... CAAC  AC  3



# Pol eta

Good recent review: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3630503/

- pol eta is goot at inserting AA opposite damaged TT
- ... but people back through history with deficient pol eta will have excess mutations
- pol iota MOSTLY puts a G opposite T instead of an A !?!?

From Myron's review:

  - "DNA polymerase eta catalyzes error-free replication across a TT cis-syn
photodimer by incorporating two A nucleotides" 
    - "responsible for error-prone incorporation within TAA motifs that generates mutations in the
variable region of immunoglobulin genes" 
    - in yeast "responsible for suppressing UV mutations by copying TC (6 – 4) and CC (6 – 4) photoproducts accurately" 
    - has "misincorporation frequency 0.1 to 0.001"

from "Somatic mutation hotspots correlate with DNA polymerase eta error spectrum"

- pol-eta hotspots in IG somatic hypermutability: mutation at G in RGYW  and middle base in TAA:
- R = A/G # Y = C/T # W = A/T

so this is 
```
AGCA / TGCT
AGCT / AGCT
AGTA / TACT
AGTT / AACT
GGCA / TGCC
GGCT / AGCC
GGTA / TACC
GGTT / AACC
.^.. / ..^.
```

# Proposal:

Add parameters for all mutations near CPDs (see below);
this mimcs errors in repair at UV-induced damage.
This is in base-model-plus-cpg-and-CPD.json .


# Terminiology

Cyclobutane pyrimidine dimer
    (CPD). The most common ultraviolet light photoproduct, in which two adjacent pyrimidines are joined together by a cyclobutane ring between their 5 and 6 positions. T-T dimers are the most abundant CPDs, although C-T, T-C and C-C dimers are also possible. Several isomers of a CPD can exist. The most common isomer is a cis-syn CPD.
    --> so, this includes with reverse complements, TT, AA, TC, GA, CT, AG, CC, and GG.

Activation-induced deaminase
    (AID). An enzyme that hydrolytically deaminates cytosines in nucleic acids, resulting in the substitution of cytosines with uracils.
    "Although AID can only deaminate dC to dU, its action gives rise to mutations at all four bases in a series of reactions that crucially depend on the Y-family polymerases114. The dU formed by the action of AID is removed by uracil DNA glycosylase (UNG), resulting in an abasic site. Direct replication of this abasic site involves REV1 and generates mutations at dG-dC base pairs52,54,55. Recognition of dU can also result in the formation of a single-strand gap, and the filling of these gaps by Pol η results in mutations at dA–dT base pairs."


# Residuals


Top 4-2-1 resids:
```
gTCg    GA           4.00            0.0345         3.96553     10.679123
cAAc    TC           7.75            0.1196         7.63036     11.029893
cGGg    CT           9.25            0.1400         9.11003     12.175151

# other
caGt    aA         821.50          421.8388       399.66118      9.729459
acGg    cT         140.25           31.0960       109.15404      9.787186

# G->C | C->G
cGtc    Ct         135.00           29.1996       105.80044      9.789697
aaCg    aG         160.00           32.9841       127.01587     11.057973
ccGt    cC         169.75           33.1254       136.62461     11.869116
aCgg    Gg         175.25           33.6337       141.61635     12.209457
cGtg    Ct         176.75           45.8013       130.94867      9.674579
gaCg    aG         162.50           29.2444       133.25565     12.320679
tcGt    cC         152.25           26.0346       126.21542     12.368218
aCga    Gg         158.75           26.5486       132.20141     12.828779

# A->G | T->C
tAtg    Gt         682.00          332.2221       349.77791      9.595070
caTa    aC         746.75          332.3065       414.44352     11.367526
caAt    aG         893.50          368.5698       524.93020     13.671359
aTtg    Ct         919.25          368.1682       551.08182     14.360282
ctAt    tG         816.75          294.5638       522.18623     15.212681
aTag    Ca         867.00          296.5920       570.40803     16.560595

tTaa    Aa         364.50          111.3741       253.12587     11.992628
ttAa    tT         375.00          111.3741       263.62587     12.490099
aAAt    TT          10.25            0.2348        10.01515     10.333207
tTTt    AA          20.50            0.2953        20.20472     18.590995
aAAa    TT          21.75            0.2964        21.45358     19.702190
tTTa    AA          30.75            0.1874        30.56259     35.299033
tAAa    TT          33.50            0.1879        33.31206     38.420136
```

Top 3-1-1 resids:
```
ACT     T         2111    1441.7          668.958        10.17198
TGT     A         2632    1867.6          764.014        10.20698
AGT     A         2147    1438.5          708.175        10.78033
ACA     T         2677    1865.7          811.435        10.84609
ATG     C         3681    2567.3         1113.360        12.68644
CAT     G         3688    2566.9         1121.060        12.77498
TCG     G          328     102.6          224.913        12.82053
CGA     C          329     103.3          226.095        12.84461
CGG     C          427     150.4          276.712        13.02641
CCG     G          433     150.6          282.790        13.30497
CGC     T          390     125.7          264.082        13.60086
GCG     G          434     142.9          290.934        14.04928
GCG     A          400     125.2          275.042        14.19135
CGC     C          446     143.5          302.533        14.58265
TAT     G         2529    1507.3         1021.800        15.19503
ACG     A          403     112.9          289.679        15.73694
ATA     C         2599    1506.1         1092.521        16.25330
CGT     T          419     111.9          306.958        16.75210
CGG     T          515     146.7          368.175        17.55004
CCG     A          539     146.9          391.757        18.66352
CGT     C          564     127.8          436.341        22.28569
ACG     G          573     129.0          444.416        22.59451
```

Top 5-1-2 resids:
```
AGCGT   G          43.000           9.120         33.88010       5.017237
CCGCT   C          50.333          11.697         38.63666       5.052234
CTGTC   C         132.333          50.757         81.57669       5.120762
ACGTG   C          56.333          13.734         42.59918       5.140618
GCCGA   A          38.000           7.100         30.89957       5.185910
CAGTT   A         238.333         113.971        124.36274       5.209654
TCGTC   T          31.000           4.955         26.04462       5.232323
GATTG   C         161.333          65.844         95.48914       5.262720
ACGGA   T          42.000           8.232         33.76760       5.263226
ACATG   G         322.333         168.813        153.52059       5.284197
CAATA   G         204.333          91.360        112.97307       5.285807
CCATG   G         331.667         172.176        159.49114       5.435830
CTGTA   A         221.000          99.291        121.70898       5.462389
TACAG   T         220.000          98.540        121.45978       5.471933
TCGTA   C          31.000           4.621         26.37946       5.488254
AACGG   G          41.667           7.557         34.11011       5.549273
CTATT   G         229.667         103.280        126.38640       5.561688
CCGTA   C          32.000           4.742         27.25781       5.597789
CATAC   C         155.000          58.701         96.29935       5.621043
CCGCA   C          54.000          11.346         42.65364       5.662952
CATGT   C         336.667         169.674        166.99264       5.733300
CAGTG   A         303.333         146.624        156.70931       5.787716
CATGG   C         339.333         170.250        169.08366       5.795267
CCGTT   C          43.000           7.337         35.66338       5.888292
AATAG   C         239.333         103.748        135.58508       5.953008
TCGTC   C          37.333           5.658         31.67540       5.955365
CCCGG   G          71.333          16.495         54.83840       6.038437
TGCGG   G          56.333          11.150         45.18354       6.051484
CCGGG   C          73.000          16.628         56.37222       6.182483
CCGTC   C          48.333           8.284         40.04983       6.223126
CTATC   G         157.667          54.442        103.22444       6.256474
ACCGT   A          43.000           6.707         36.29278       6.267065
TCGTT   C          49.000           7.436         41.56412       6.816590
AACGA   G          51.000           7.772         43.22838       6.934706
CTATA   G         195.667          66.170        129.49689       7.119418
GACGG   G          58.333           8.439         49.89412       7.680925
GACGA   G          47.333           5.719         41.61410       7.781913
TATAG   C         210.667          67.527        143.13950       7.789960
CTATG   G         242.667          74.647        168.01961       8.696983
CAATG   G         329.667         104.972        224.69427       9.807751
CATTG   C         352.667         104.822        247.84513      10.826053
CATAG   C         296.667          75.434        221.23265      11.391492
TTAAA   T         277.333          45.169        232.16442      15.448661
TTTAA   A         279.000          45.217        233.78258      15.547990
```

