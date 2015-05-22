Questions to demonstrate with simulation
----------------------------------------

1. How does accuracy increase with increasing tuple width?
2. How long a time can we get at (as a function of tuple width)?


Simulation plan
---------------

1. Correctness of posterior distribution : TASEP
    Goal : boxplots of posterior distributions and percent within 95% credible intervals

    Parameters : 

        - sequence length : 1e4 

        - repications per value : 100 

        - scaled time : [ 0.25, 1.5 ]

        - tuple widths :

            - [ 4, 2 ]
            - [ 6, 2 ]


2. Ising model: tuple width and total amount of time.
    Goal: boxplots of posterior distributions across different parameter values, showing decreasing precision with increasing time and decreasing window length

    - Sequence length : 1e4 

    - replications per value : 1
    
    - Parameters :
        
        - inverse temperature : strong
        - magnetization : moderate
        - scaled time : [ 0.125, 0.25, 1, 4 ]

    - Tuple widths [ long, short ] : 

        - [ 2, 2 ]
        - [ 4, 2 ]
        - [ 5, 2 ]
        - [ 7, 2 ]
        - [ 3, 3 ]
        - [ 5, 3 ]
        - [ 7, 3 ]
        - [ 9, 3 ]

3. CpG+[CG-bias] model : how much divergence is feasible to estimate with in a nucleotide model?
    Goal: boxplots of posterior distribution of CpG and [CG] parameters across different total times

    - Sequence length : 1e4

    - replications per value : 1

    - Parameters :

        - total time : [ 0.125, 0.25, 0.5, 1, 2 ]
        - single-nuc transitions : all equal to 1.0
        - CpG : 3.0
        - CG bias (2Ns) : 0.01

    - Tuple widths [ long, short ] :

        - [ 3, 1 ]
        - [ 5, 3 ]

