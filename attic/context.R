#!/usr/bin/R
# input data

bases <- c("A","C","G","T")
nbases <- length(bases)

codons <- data.frame( rbind(
        c("TTT", "F", "Phe", "TT[TC]"),
        c("TTC", "F", "Phe", "TT[TC]"),
        c("TTA", "L", "Leu", "TT[AG]"),
        c("TTG", "L", "Leu", "TT[AG]"),
        c("CTT", "L", "Leu", "CT[ACGT]"),
        c("CTC", "L", "Leu", "CT[ACGT]"),
        c("CTA", "L", "Leu", "CT[ACGT]"),
        c("CTG", "L", "Leu", "CT[ACGT]"),
        c("ATT", "I", "Ile", "AT[TCA]"),
        c("ATC", "I", "Ile", "AT[TCA]"),
        c("ATA", "I", "Ile", "AT[TCA]"),
        c("ATG", "M", "Met", "ATG"),
        c("GTT", "V", "Val", "GT[TCAG]"),
        c("GTC", "V", "Val", "GT[TCAG]"),
        c("GTA", "V", "Val", "GT[TCAG]"),
        c("GTG", "V", "Val", "GT[TCAG]"),
        c("TCT", "S", "Ser", "TC[TCAG]"),
        c("TCC", "S", "Ser", "TC[TCAG]"),
        c("TCA", "S", "Ser", "TC[TCAG]"),
        c("TCG", "S", "Ser", "TC[TCAG]"),
        c("CCT", "P", "Pro", "CC[TCAG]"),
        c("CCC", "P", "Pro", "CC[TCAG]"),
        c("CCA", "P", "Pro", "CC[TCAG]"),
        c("CCG", "P", "Pro", "CC[TCAG]"),
        c("ACT", "T", "Thr", "AC[TCAG]"),
        c("ACC", "T", "Thr", "AC[TCAG]"),
        c("ACA", "T", "Thr", "AC[TCAG]"),
        c("ACG", "T", "Thr", "AC[TCAG]"),
        c("GCT", "A", "Ala", "GC[TCAG]"),
        c("GCC", "A", "Ala", "GC[TCAG]"),
        c("GCA", "A", "Ala", "GC[TCAG]"),
        c("GCG", "A", "Ala", "GC[TCAG]"),
        c("TAT", "Y", "Tyr", "TA[TC]"),
        c("TAC", "Y", "Tyr", "TA[TC]"),
        c("TAA", "stop", "Ter", "TA[AG]"),
        c("TAG", "stop", "Ter", "TA[AG]"),
        c("CAT", "H", "His", "CA[TC]"),
        c("CAC", "H", "His", "CA[TC]"),
        c("CAA", "Q", "Gln", "CA[AG]"),
        c("CAG", "Q", "Gln", "CA[AG]"),
        c("AAT", "N", "Asn", "AA[TC]"),
        c("AAC", "N", "Asn", "AA[TC]"),
        c("AAA", "K", "Lys", "AA[AG]"),
        c("AAG", "K", "Lys", "AA[AG]"),
        c("GAT", "D", "Asp", "GA[TC]"),
        c("GAC", "D", "Asp", "GA[TC]"),
        c("GAA", "E", "Glu", "GA[AG]"),
        c("GAG", "E", "Glu", "GA[AG]"),
        c("TGT", "C", "Cys", "TG[TC]"),
        c("TGC", "C", "Cys", "TG[TC]"),
        c("TGA", "stop", "Ter", "TGA"),
        c("TGG", "W", "Trp", "TGG"),
        c("CGT", "R", "Arg", "CG[TCAG]"),
        c("CGC", "R", "Arg", "CG[TCAG]"),
        c("CGA", "R", "Arg", "CG[TCAG]"),
        c("CGG", "R", "Arg", "CG[TCAG]"),
        c("AGT", "S", "Ser", "AG[TC]"),
        c("AGC", "S", "Ser", "AG[TC]"),
        c("AGA", "R", "Arg", "AG[AG]"),
        c("AGG", "R", "Arg", "AG[AG]"),
        c("GGT", "G", "Gly", "GG[TCAG]"),
        c("GGC", "G", "Gly", "GG[TCAG]"),
        c("GGA", "G", "Gly", "GG[TCAG]"),
        c("GGG", "G", "Gly", "GG[TCAG]")
    ) )
names(codons) <- c("codon", "aacode", "aa", "regexp")

synons <- setdiff( levels(codons$aa), c("Met","Trp","Ter") )

