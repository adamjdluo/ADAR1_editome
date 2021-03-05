# ADAR1_editome
The scripts are applied to the study of ADAR1 editome, published on PNAS[link]().
- Align_and_varCall.R: follow the step of GATK best practice, variants were annotated with annovar
- rawCount_AF_AC.R: calculate depth (DP), alleleic fraquency (AF), and alternative alleleic counts (AC) for each variant
- gather_and_proecss_varInfo.R: integrate information of each variant into single table. filtering and grouping variants
- report_rmd.Rmd: Make visualized reports
