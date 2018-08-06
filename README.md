# Prediction of overexpressed miRNA in given transcriptome data

I coded this program for Introduction to Computer Science for Biologists course project (Spring 2016 semester, Prof. Daehyun Baek).

This program finds an overexpressed miRNA in a given transcriptome data.
The program inputs gene expression fold change values and outputs the overexpressed miRNA.

For example, the program should output miR-1 if inputting gene expression fold change values of miR-1 transfected HeLa cells.

The program first finds miRNA motifs that are enriched in down-regulated genes in miR-1 transfected Hela cells.
Then, the program reports miRNAs having that miRNA motifs, so in this case, it should report miR-1.
