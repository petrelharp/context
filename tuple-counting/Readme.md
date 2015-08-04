---
title: Tools for counting paired tuples
author: Peter Ralph
date: 3 August 2015
---

*General-purpose:*

codons.py

: List of codons.

count-paired-tuples.py

: Counts paired tuples from a paired fasta or an axt file.

count-tuples.sh

: Script to run count-paired-tuples.py on everything in a directory.

get-regions-from-alignment.py

: Write out the subset of a .axt file
given by an input file of start and end positions (two whitespace-separated columns).
Positions are in the "primary" organism, i.e. the one listed first.

get-positions.py

: Write out a collection of files giving the locations of certain sequences.

is-coding.py

: Split up position files into coding, noncoding etc.

plrutils.py

: Utilities.

sum-counts.py

: Combine many 'count' tables into one.

*Primate data:*

explore-tuple-counts.R

: Data exploration of human-chimp-gorilla tuples.

*Drosophila data:*

makefile

: Stuff for processing Drosophila data.

get-all-positions.sh

: runs get-positions.py on the Drosophila data

is-coding.sh

: Runs is-coding.py on the Drosophila data.

extract-coding.sh

: Extracts coding sequence from Drosophila data
