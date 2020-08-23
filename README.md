# SLiM_scripts

SLiM scripts for looking at how synergistic epistasis and admixture can affect patterns of LD

run_slim.sh
bash file for running many parrallel instances of SLiM script


slim_admixture_3pops.eidos
SliM script simulating admixture between 3 subpopulations, timing of admixture can be varied, outputs vcf of a sample of individuals

slim_synergistic_epistasis.eidos
SLiM script simulating negative synergistic epistasis between deleterious mutations
Future could be combined with admixture to test joint effect on patterns of LD at neutral vs. selected sites

LD_calc_ref_adj_PCA.R
R script for calculating LD from list of SLiM outputed vcf's
Calculated LD separately for deleterious and neutral mutations
Has filtering step to remove recent migrants ny performing genomic PCA (neccessary for admixture simulations)
Also mimics LD calculation from real world data by setting alleles with over 50% frequency as the reference allele by default
