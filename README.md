# Prioritizing genes associated with brain disorders by leveraging enhancer-promoter interactions in diverse neural cells and tissues
These scirpts are for the reconstructing enhancer-promoter interactions (EPIs) in 958 FANTOM5 samples, downstream analysis of the cell- and tissue-specificity of the EPIs in human neural cell and tissue types, and prioritizing the critical genes associated with diverse brain disorders and behavioral-cognitive phenotypes based on the EPIs.

___
Please install the dependency for these scripts

- Python v2.7
- Perl v5.16
- R v4.1.0 
- HOMER v4.10
- LDSC (LD SCore) v1.0.1
- MAGMA v1.08

___
Here, we organized the custom codes of the computational analyses in our study into 6 parts as shown below.

- Part 1. <b>Predicting_EPIs</b>.
We predicted the set of active EPIs in each of the 439 human cell and tissue types, and quantified the strength of each active EPI in each cell and tissue type. All the results can be obtained via https://soulnature.github.io/brainepl/EPIs.html.

- Part 2. <b>EPI_validation</b>.
We estimated the accuracy (AUPR scores) of the reconstructed EPIs in the cell and tissue types using GTEx eQTL and pcHi-C datasets from the matched cell or tissue types.

- Part 3. <b>Specificity_analysis</b>.
We performed differential expression analysis on the enhancers and promoters in neural cell and tissue types to determine their cell- and tissue-specificity. We also predicted enriched TF binding on the enhancer regions with different cell- and tissue-specificity.

- Part 4. <b>Partitioned_LDSC</b>.
We estimated the partitioned heritability of 17 brain disorders and behavioral-cognitive phenotypes in the active regulatory elements (i.e., enhancers and promoters) from each of the neural cell and tissue type.

- Part 5. <b>Prioritizing_genes</b>.
We prioritized a set of critical genes associated with each of the brain disorders and behavioral-cognitive phenotypes for each neural cell and tissue type by leveraging the reconstructed EPIs. The prioritized genes of 17 brain disorders and behavioral-cognitive phenotypes can be obtained via https://soulnature.github.io/brainepl/prioritized_genes.html.

- Part 6. <b>Gene_analysis</b>.
Finally, we analyzed the pleiotropy and gene expression patterns of the prioritized genes of each brain disorder and behavioral-cognitive phenotype.

___
Supplementary website: https://soulnature.github.io/brainepl

___
If you have any questions, please contact Drs. Xing-Ming Zhao (xmzhao AT fudan.edu.cn) or Yucheng T. Yang (yangyy AT fudan.edu.cn).

___
If you use the results in your study, please cited
> Zhao et al. Prioritizing genes associated with brain disorders by leveraging enhancer-promoter interactions in diverse neural cells and tissues.

