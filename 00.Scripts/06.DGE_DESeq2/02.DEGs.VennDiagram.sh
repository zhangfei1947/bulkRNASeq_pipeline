#!/bin/bash

module restore R431
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

out_path="/scratch/group/lilab/Phil/20250107_isogrp/06.DGE_DESeq2"

sig_file_1=body_d1_iso.grp.logFC1.padj0.05.csv
sig_file_2=body_d7_iso.grp.logFC1.padj0.05.csv
sig_file_3=body_grp_d7.d1.logFC1.padj0.05.csv
sig_file_4=body_iso_d7.d1.logFC1.padj0.05.csv

sig_file_5=head_d1_iso.grp.logFC1.padj0.05.csv
sig_file_6=head_d7_iso.grp.logFC1.padj0.05.csv
sig_file_7=head_grp_d7.d1.logFC1.padj0.05.csv
sig_file_8=head_iso_d7.d1.logFC1.padj0.05.csv


Rscript 02.DEGs.VennDiagram.R  ${sig_file_1},${sig_file_2}  body_d1_iso.grp,body_d7_iso.grp $out_path

Rscript 02.DEGs.VennDiagram.R  ${sig_file_3},${sig_file_4}  body_grp_d7.d1,body_iso_d7.d1 $out_path


Rscript 02.DEGs.VennDiagram.R  ${sig_file_5},${sig_file_6}  head_d1_iso.grp,head_d7_iso.grp $out_path

Rscript 02.DEGs.VennDiagram.R  ${sig_file_7},${sig_file_8}  head_grp_d7.d1,head_iso_d7.d1 $out_path

