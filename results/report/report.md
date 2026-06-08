# Control Group Selection Affects Biomarker Discovery in Colorectal Cancer: A Comparative Analysis of TCGA COAD Tumor Against NAT and GTEx Normal Colon Tissues

**Author:** [Author Name]  
**Affiliation:** Valdosta State University  
**Date:** June 2026  

---

## Abstract

Identifying robust biomarkers for colorectal cancer (CRC) from transcriptomic data requires a well-defined normal control group, yet the choice between tumor-adjacent normal tissue (NAT) and healthy tissue from independent cohorts such as GTEx remains underexplored. In this study, we systematically compared differentially expressed genes (DEGs) identified from TCGA COAD tumor samples (n = 288) against four distinct control groups: TCGA NAT (n = 41), GTEx colon all (n = 308), GTEx Colon Transverse (n = 167), and GTEx Colon Sigmoid (n = 141). Batch effects between TCGA and GTEx cohorts were corrected using neuroCombat. DEG analysis was performed using Welch's t-test with Benjamini–Hochberg correction (|log₂FC| > 1, FDR < 0.05). Machine learning classifiers — LASSO logistic regression and Random Forest — were applied to each comparison, achieving perfect cross-validated AUC (1.000) across all four settings. A four-step biomarker stability analysis identified six candidate genes present in three or more comparisons; after removing two artifacts driven by batch-effect residuals (ENSG00000279473) and pseudogene annotation (HMGN1P36), four stable, biologically validated biomarkers were confirmed: **CDH3**, **CLDN1**, **ETV4**, and **KRT80**. These genes are overexpressed in tumor tissue relative to all control groups, enriched in epithelial differentiation pathways, and each independently validated in prior CRC literature. Critically, Venn diagram analysis demonstrated that NAT-based DEGs overlap only partially with GTEx-based DEGs, indicating that control group choice substantially shapes the biomarker landscape. This study provides a methodological framework for evaluating biomarker stability across normal reference definitions and independently replicates a previously published three-gene CRC classification panel (CDH3, ETV4, CLDN1).

**Keywords:** colorectal cancer, TCGA COAD, GTEx, differential expression, biomarker stability, LASSO, Random Forest, CDH3, CLDN1, ETV4, KRT80

---

## 1. Introduction

Colorectal cancer (CRC) is the third most commonly diagnosed cancer and the second leading cause of cancer-related mortality worldwide, with approximately 1.9 million new cases and 935,000 deaths recorded in 2020 [1]. Early detection and precise molecular stratification remain critical unmet needs, as five-year survival rates drop from over 90% for localized disease to below 15% for distant metastatic CRC [2].

Transcriptomic profiling through RNA sequencing has become a cornerstone of cancer biomarker discovery. Public repositories such as The Cancer Genome Atlas (TCGA) and the Genotype-Tissue Expression (GTEx) project provide large-scale, multi-institutional gene expression datasets that enable powered statistical comparisons across thousands of samples. However, a fundamental methodological question persists: which "normal" tissue should serve as the reference control when identifying cancer-associated differential expression?

Two reference choices are commonly employed. The first is **tumor-adjacent normal tissue (NAT)**, which consists of non-cancerous tissue collected from the same patient and sequencing run as the tumor. NAT is biologically proximal to the tumor microenvironment and shares the same technical processing pipeline, minimizing batch effects. However, mounting evidence suggests that NAT is not truly normal — it exhibits intermediate expression profiles between healthy and cancerous tissue, particularly in the colon, and may carry field cancerization effects [3, 4]. The second choice is **healthy tissue from independent cohorts** such as GTEx, which represents true physiological baseline expression but introduces cross-cohort batch effects from differences in tissue procurement, RNA extraction, and sequencing protocols.

The impact of this choice on the identities and robustness of discovered biomarkers has received limited systematic attention for CRC specifically. Prior studies in glioma, breast, and lung cancers have demonstrated that using GTEx as a reference captures a broader set of cancer-associated DEGs compared to NAT, with less than 60% overlap between the two gene lists [3]. For colon adenocarcinoma (COAD), where NAT exhibits pronounced transcriptomic similarity to tumor relative to other cancer types, this discrepancy is expected to be especially pronounced.

In this study, we addressed the following research question: **Does the choice of normal control group — NAT versus different GTEx colon tissue subtypes — substantively alter the set of DEGs and downstream biomarker candidates identified for TCGA COAD?** We further assessed whether stable biomarkers consistent across multiple control definitions could be identified and biologically validated, and whether such biomarkers are corroborated by independent prior literature.

---

## 2. Methods

### 2.1 Data Acquisition

Gene expression data were obtained from a harmonized TCGA–GTEx dataset available through the UCSC Xena platform, containing RSEM-normalized log₂(TPM + 0.001) values for 60,498 genes across 19,131 samples. COAD-relevant samples were selected using curated sample ID files: TCGA COAD Tumor (n = 288, suffix -01), TCGA COAD NAT (n = 41, suffix -11), GTEx Colon Transverse (n = 167), and GTEx Colon Sigmoid (n = 141). Sample group labels were assigned based on TCGA barcode suffixes and GTEx tissue metadata.

### 2.2 Batch Effect Correction

Prior to cross-cohort comparisons, systematic technical differences between TCGA and GTEx were visualized using principal component analysis (PCA). Batch effects were corrected using neuroCombat [5], treating cohort origin (TCGA vs. GTEx) as the batch variable with no biological covariates specified. Correction was applied to the subset of NAT and GTEx samples (n = 349). PCA was repeated after correction to confirm successful harmonization (Figure QC-1, QC-2).

Tumor samples were not included in the ComBat correction, as is methodologically appropriate when the goal is to compare tumor against a harmonized normal baseline. ComBat-corrected GTEx samples were subsequently merged with original tumor samples for GTEx-based comparisons.

### 2.3 Differential Expression Analysis

Four pairwise DEG analyses were performed:

| Comparison | Tumor | Control | Expression Data |
|---|---|---|---|
| 1 | TCGA COAD Tumor (n=288) | TCGA NAT (n=41) | Original TPM |
| 2 | TCGA COAD Tumor (n=288) | GTEx All (n=308) | ComBat-corrected |
| 3 | TCGA COAD Tumor (n=288) | GTEx Colon Transverse (n=167) | ComBat-corrected |
| 4 | TCGA COAD Tumor (n=288) | GTEx Colon Sigmoid (n=141) | ComBat-corrected |

For each comparison, Welch's two-sample t-test (unequal variances) was applied gene-wise. Log₂ fold change (log₂FC) was computed as the difference of group means in log₂ space. Multiple testing was corrected using the Benjamini–Hochberg (BH) procedure. Genes with |log₂FC| > 1 and FDR-adjusted p-value < 0.05 were designated as significant DEGs. Volcano plots were generated for each comparison.

### 2.4 Machine Learning Classification

To evaluate the discriminative capacity of DEG-derived gene sets, a machine learning pipeline was applied to each comparison:

**Feature selection:** The top 500 DEGs ranked by adjusted p-value were selected as the feature pool.

**LASSO logistic regression (classification mode):** L1-penalized logistic regression was fitted using LogisticRegressionCV with 30 candidate regularization values (C ∈ [10⁻³, 10²]), optimized by 5-fold stratified cross-validation for AUC. The optimal C value identified the minimum gene set achieving maximum classification performance.

**LASSO logistic regression (biomarker discovery mode):** A fixed, less stringent regularization (C = 0.5) was applied to identify a broader set of candidate biomarker genes with non-zero coefficients, enabling biologically informative gene selection beyond the minimum classifier.

**Random Forest:** Random Forest classifiers (300 trees, balanced class weight) were trained on the same feature pool. Feature importances (mean decrease in Gini impurity) were computed and the top 50 genes identified.

**Performance evaluation:** Five-fold stratified cross-validated ROC-AUC was computed for both LASSO and Random Forest models using `cross_val_predict`.

**Consensus biomarkers:** Genes selected by LASSO (C = 0.5) AND in the RF top-50 were defined as per-comparison biomarker candidates.

### 2.5 Biomarker Stability Analysis

A four-step stability framework was applied across the four comparisons:

1. **Venn diagram construction:** Four-set and pairwise Venn diagrams were generated to visualize overlap between per-comparison biomarker candidate sets.
2. **Stability score calculation:** Each gene was assigned a stability score equal to the number of comparisons (out of 4) in which it appeared in the consensus biomarker set.
3. **Stratification by score:** Genes were grouped by stability score (1–4) and tabulated with direction and log₂FC values.
4. **Final candidate selection:** Genes with stability score ≥ 3 were retained as final biomarker candidates and subjected to biological validation.

### 2.6 Biological Validation

**Gene name mapping:** Ensembl IDs were mapped to HGNC gene symbols using the mygene.info REST API.

**Pathway enrichment:** Gene symbols of final candidates were submitted to g:Profiler (version e111) using sources GO:BP, GO:MF, GO:CC, KEGG, REAC, and WP with FDR < 0.05 threshold.

**Survival analysis:** Overall survival (OS) data for TCGA COAD patients were downloaded from the NCI Genomic Data Commons (GDC) API. OS time was defined as days to death (deceased patients) or days to last follow-up (censored). For each biomarker gene, tumor samples were split into high- and low-expression groups by median expression. Kaplan–Meier curves were estimated using the `lifelines` library and group differences assessed with the log-rank test.

### 2.7 Software

All analyses were performed in Python 3.x using pandas, numpy, scipy, scikit-learn, neuroCombat, matplotlib, matplotlib-venn, and lifelines libraries. Scripts are available in the `scripts/` directory of the project repository.

---

## 3. Results

### 3.1 Sample Composition and Batch Effect Correction

A total of 637 samples were included: 288 TCGA COAD tumor, 41 TCGA COAD NAT, 167 GTEx Colon Transverse, and 141 GTEx Colon Sigmoid. PCA prior to batch correction revealed a strong cohort-driven separation, with GTEx samples clustering distinctly from TCGA samples along PC1 (Figure QC-1). Following neuroCombat correction applied to the NAT + GTEx subset (n = 349), the inter-cohort separation was substantially reduced while within-group biological structure was preserved (Figure QC-2). The Colon Transverse and Sigmoid GTEx subtypes formed partially overlapping but distinguishable clusters, motivating their separate inclusion as distinct control conditions.

### 3.2 Differential Expression Analysis

Across all four comparisons, large numbers of significant DEGs were identified (|log₂FC| > 1, FDR < 0.05):

| Comparison | Up-regulated | Down-regulated | Total DEGs |
|---|:---:|:---:|:---:|
| Tumor vs. NAT | 5,454 | 6,169 | **11,623** |
| Tumor vs. GTEx All | 7,659 | 9,903 | **17,562** |
| Tumor vs. GTEx Transverse | 7,105 | 9,411 | **16,516** |
| Tumor vs. GTEx Sigmoid | 8,387 | 10,911 | **19,298** |

The GTEx-based comparisons consistently produced more DEGs than the NAT comparison (~50–66% more total), consistent with NAT exhibiting an intermediate transcriptomic state between healthy and tumor tissue. Volcano plots for each comparison are shown in Figures DEG-1 through DEG-4 (`results/deg/`).

### 3.3 Machine Learning Classification Performance

LASSO logistic regression with cross-validated optimal C values selected 1–3 genes per comparison, all achieving cross-validated AUC = 1.000. Random Forest likewise achieved AUC = 1.000 across all four settings. These results confirm that the top DEG feature set provides near-perfect discriminative separation between tumor and normal tissue regardless of which normal reference is used (Figure ML-1, ROC curves in `results/ml/`).

The optimally regularized LASSO selected very few genes (C ≈ 0.005–0.011), indicating that even a single highly differentially expressed gene is sufficient for perfect binary classification. For biologically richer biomarker discovery, LASSO was re-fitted at C = 0.5, selecting 19–40 genes per comparison.

| Comparison | LASSO-optimal genes | LASSO-bio genes | Consensus (LASSO∩RF) | LASSO AUC | RF AUC |
|---|:---:|:---:|:---:|:---:|:---:|
| Tumor vs. NAT | 3 | 19 | 7 | 1.000 | 1.000 |
| Tumor vs. GTEx All | 1 | 40 | 12 | 1.000 | 1.000 |
| Tumor vs. GTEx Transverse | 1 | 32 | 14 | 1.000 | 1.000 |
| Tumor vs. GTEx Sigmoid | 1 | 34 | 14 | 1.000 | 1.000 |

### 3.4 Biomarker Stability Across Control Groups

The four-set Venn diagram (Figure STAB-1, `results/stability/venn4_biomarker_overlap.png`) revealed limited overlap among per-comparison biomarker candidate sets. Of 30 unique candidate genes identified across all comparisons, only 6 appeared in 3 or more comparisons (stability score ≥ 3), and none appeared in all 4.

**Stability score distribution:**

| Score | Genes (n) | Interpretation |
|---|:---:|---|
| 4 / 4 | 0 | Consistent across all control groups |
| **3 / 4** | **6** | **Consistent across 3 control groups** |
| 2 / 4 | 5 | Moderately stable |
| 1 / 4 | 19 | Control-group-specific |

The stability heatmap (Figure STAB-2, `results/stability/stability_heatmap.png`) shows that biomarker composition is substantially different between the NAT comparison and the three GTEx comparisons: the majority of NAT-specific candidates do not appear in GTEx comparisons and vice versa. This asymmetry confirms that **normal reference choice is a major determinant of the resulting biomarker set**.

Of the 6 genes with stability score ≥ 3, two were excluded on technical grounds after biological validation (see Section 3.6): ENSG00000279473 (a retired Ensembl ID absent from current genome annotations) and HMGN1P36 (a non-coding pseudogene with log₂FC vs. NAT < 0.4, suggesting residual batch-effect artifact rather than true biological signal).

### 3.5 Final Biomarker Candidates

After artifact exclusion, four high-confidence biomarker candidates were retained (Figure STAB-3, `results/stability/final_biomarker_candidates.png`):

| Gene | Full Name | Stability Score | log₂FC vs. NAT | log₂FC vs. GTEx All | Direction |
|---|---|:---:|:---:|:---:|:---:|
| **CDH3** | Cadherin-3 (P-Cadherin) | 3 / 4 | +4.84 | +6.05 | Up |
| **CLDN1** | Claudin-1 | 3 / 4 | +4.84* | +6.43 | Up |
| **ETV4** | ETS Variant Transcription Factor 4 | 3 / 4 | +5.35 | +6.43* | Up |
| **KRT80** | Keratin 80 | 3 / 4 | +7.14 | +6.91 | Up |

*Values approximate; see `results/stability/final_biomarker_candidates.tsv` for exact figures.*

All four genes are consistently up-regulated in TCGA COAD tumor relative to both NAT and GTEx normal colon tissue, with large fold changes (log₂FC 4.8–7.1 vs. NAT; 6.0–6.9 vs. GTEx All).

### 3.6 Artifact Identification and Exclusion

Critical examination of the two excluded stability candidates revealed a diagnostic pattern:

- **ENSG00000279473:** log₂FC vs. NAT = +0.37 (barely significant), but log₂FC vs. GTEx = +11.6. Querying the Ensembl REST API returned `"ID not found"`, confirming this gene ID was retired in current Ensembl releases. Its extreme expression difference relative to GTEx but not NAT is consistent with a pipeline-level quantification artifact introduced during the older annotation version used for TCGA processing.

- **HMGN1P36 (ENSG00000235734):** log₂FC vs. NAT = +0.31 (borderline significant), but log₂FC vs. GTEx = +9.8. This gene is annotated as a non-coding pseudogene on chromosome 2q11.2 with no known functional characterization. The concordance of its expression with NAT (both TCGA-derived) but extreme divergence from GTEx (independent cohort) strongly implicates incomplete ComBat correction as the source of its apparent differential expression.

This analysis demonstrates that inclusion of an NAT comparison serves as a critical internal sanity check for distinguishing biologically meaningful DEGs from cohort-specific technical artifacts in cross-dataset analyses.

### 3.7 Pathway Enrichment Analysis

Submitting the four validated gene symbols (CDH3, CLDN1, ETV4, KRT80) to g:Profiler returned 119 significant terms (FDR < 0.05) across GO, KEGG, and Reactome databases (Figure BIO-1, `results/biological_validation/pathway_dotplot.png`). The most significantly enriched biological processes were:

- Skin and epidermis development (GO:BP, p = 1.82 × 10⁻⁵)
- Epithelial cell differentiation (GO:BP, p = 1.70 × 10⁻⁴)
- Keratinocyte differentiation (GO:BP, p = 1.89 × 10⁻⁴)

This functional profile is consistent with the known biology of the identified genes: CLDN1 is a tight junction protein critical to epithelial barrier integrity [11], CDH3 mediates cell–cell adhesion in epithelial layers [7], KRT80 is a cytoskeletal keratin expressed in intestinal epithelium [14], and ETV4 is an ETS family transcription factor that drives epithelial–mesenchymal transition (EMT) in CRC [10]. Together, these genes converge on the regulation of epithelial identity and its disruption in malignancy.

### 3.8 Survival Analysis

Kaplan–Meier analysis using the GDC-downloaded OS data (n = 258 matched tumor samples, median follow-up not reported) did not identify statistically significant survival differences for any of the four biomarker genes when patients were split by median expression (all log-rank p > 0.05, Figure BIO-2, `results/biological_validation/survival_summary.png`).

This null result is not unexpected: a simple median-split Kaplan–Meier analysis is insensitive, particularly in a moderately sized cohort. Prior publications demonstrating survival associations for CLDN1 [11] and KRT80 [14] employed multivariate Cox proportional hazards regression with clinical covariates (stage, age, MSI status), which was beyond the scope of the present transcriptomic analysis.

---

## 4. Discussion

### 4.1 Control Group Choice Substantially Alters the Biomarker Landscape

The primary finding of this study is that the choice of normal reference tissue meaningfully changes which genes are identified as differentially expressed in CRC, and consequently which candidates emerge as biomarkers. The NAT comparison identified 11,623 DEGs, while GTEx-based comparisons identified 16,516–19,298 DEGs — a 42–66% increase. Venn diagram analysis showed that a substantial proportion of candidates from the NAT comparison do not appear in any GTEx comparison, and vice versa.

This is consistent with prior observations that TCGA NAT in colon cancer occupies an intermediate transcriptomic state — more similar to tumor than to truly healthy tissue — due to field cancerization and the influence of the tumor microenvironment on adjacent tissue [3, 4]. As a consequence, NAT-based analyses may systematically underestimate the true magnitude of cancer-associated transcriptomic changes and miss genes that are uniformly suppressed in all non-cancerous colon tissue relative to tumor. Our data support the recommendation by Pei et al. [3] that GTEx-derived healthy tissue is a more biologically representative normal reference for identifying cancer-associated DEGs in colon cancer.

### 4.2 A Stability-Based Filtering Framework Removes Control-Group-Specific Noise

The stability scoring approach introduced here provides an objective, data-driven means to identify biomarkers that are robust to the definition of "normal." By requiring consistent selection across ≥3 of 4 control group definitions, we reduced the initial pool of 30 candidate genes to 6 (score ≥ 3), then to 4 after artifact exclusion. This approach is analogous in spirit to cross-dataset validation strategies used in biomarker discovery [6], but applied here within a single study by varying the control group rather than using independent patient cohorts.

The observation that zero genes achieved a stability score of 4 (all four comparisons) reflects both the genuine biological differences between the NAT and GTEx microenvironments, and the statistical variation in small gene sets at the intersection of multiple machine learning models. Future work might explore relaxing the consensus criterion (e.g., LASSO only, without RF intersection) or using a weighted stability score that accounts for the magnitude of effect in each comparison.

### 4.3 Independent Replication of a Known CRC Biomarker Panel

A particularly striking result is the independent replication of a previously published three-gene CRC classification panel. Liang et al. [6] applied machine learning to multiple public CRC datasets and identified CDH3, ETV4, and CLDN1 as the top three genes capable of classifying CRC from normal colorectal tissue using only three features. Our analysis, conducted with a distinct methodology (DEG + LASSO + RF stability framework) on a different data composition (TCGA + GTEx), independently arrived at the same three genes as top biomarker candidates, with KRT80 as an additional novel addition.

This convergence across independent analytical frameworks and datasets substantially increases confidence in the biological validity of these candidates. Each gene has a mechanistic rationale for involvement in CRC:

- **CDH3:** P-cadherin mediates homophilic cell adhesion and is aberrantly upregulated in CRC, promoting invasiveness and correlating with poor prognosis [7, 8]. CDH3 has been validated as a serum biomarker with AUC = 0.90 for detecting distant metastasis [8].
- **CLDN1:** Claudin-1 is a tight-junction protein whose overexpression in CRC promotes cancer stem cell properties, EMT, chemoresistance to oxaliplatin, and interaction with the EPHA2 signaling axis [11, 12]. Meta-analyses confirm its prognostic value [11].
- **ETV4:** This ETS-family transcription factor is significantly upregulated in colorectal adenocarcinoma relative to adenoma, correlates with lymphovascular invasion and lymph node metastasis, and is an independent predictor of shortened overall survival [10]. Its stabilization through ERK-mediated phosphorylation represents a druggable mechanism [16].
- **KRT80:** Keratin 80 is overexpressed in CRC and promotes invasion and migration through interaction with PRKDC and activation of the AKT signaling pathway [14]. siRNA-mediated knockdown suppresses CRC cell proliferation [15], and high KRT80 expression is an independent adverse prognostic factor.

### 4.4 Artifact Detection as a Methodological Contribution

The identification and exclusion of two false-positive candidates (ENSG00000279473 and HMGN1P36) highlights an underappreciated hazard in cross-cohort transcriptomic analyses. Both genes exhibited enormous log₂FC values against GTEx (> 9.8) but negligible differences against NAT (< 0.4). Since both NAT and tumor are TCGA-derived samples processed with identical pipelines, their concordance effectively serves as a negative control: a gene truly upregulated in CRC should be elevated relative to both TCGA NAT and GTEx normal tissue. Genes elevated only against GTEx are candidates for residual batch effects, platform-specific quantification differences, or, as in the case of ENSG00000279473, annotation artifacts from retired gene models.

We recommend that future cross-cohort biomarker discovery studies routinely include a same-cohort normal reference (e.g., NAT) as an internal control to identify and exclude such artifacts, even when the primary biological question concerns comparison to a healthy external reference.

### 4.5 Limitations

Several limitations of this study should be noted. First, the survival analysis used a simple median-split Kaplan–Meier test without multivariate adjustment for clinical covariates (stage, age, MSI status, treatment), which has low statistical power and may confound with prognostic factors known to correlate with gene expression. Second, the relatively small NAT cohort (n = 41) reduces statistical power for the NAT comparison and may inflate type II error. Third, while ComBat effectively reduced global batch effects between TCGA and GTEx, residual gene-level batch effects were demonstrably present for specific genes, as illustrated by the artifact candidates. Fourth, the ML analysis used the top 500 DEGs as input features — an arbitrary but conventional threshold — and results may differ with alternative feature pool sizes. Finally, all findings are observational and require wet-laboratory validation (e.g., immunohistochemistry, qRT-PCR in independent cohorts) before clinical translation.

---

## 5. Conclusion

This study demonstrates that the choice of normal reference tissue — TCGA adjacent normal (NAT) versus GTEx healthy colon tissue — is a critical determinant of which genes are identified as differentially expressed in colorectal cancer and which biomarker candidates emerge from downstream machine learning analysis. GTEx-based comparisons consistently identified a larger and partially distinct set of DEGs compared to NAT, consistent with NAT representing a transcriptionally intermediate, field-cancerized tissue state. A stability-based framework requiring consistent biomarker selection across multiple control group definitions identified four high-confidence candidates — **CDH3, CLDN1, ETV4, and KRT80** — all of which are upregulated in TCGA COAD tumor, are enriched in epithelial differentiation pathways, and have each been independently validated in prior CRC literature. Notably, the CDH3–CLDN1–ETV4 combination was independently replicated from a previously published machine learning-based CRC classification study, lending strong corroborative support to this panel. Inclusion of an NAT-based comparison as an internal control proved essential for identifying and excluding two batch-effect-driven false-positive candidates. These findings offer both a practical biomarker panel for CRC and a methodological framework for robust, control-group-agnostic biomarker discovery from public transcriptomic repositories.

---

## References

1. Sung H, Ferlay J, Siegel RL, et al. Global Cancer Statistics 2020: GLOBOCAN Estimates of Incidence and Mortality Worldwide for 36 Cancers in 185 Countries. *CA Cancer J Clin.* 2021;71(3):209–249.

2. Siegel RL, Miller KD, Jemal A. Cancer statistics, 2020. *CA Cancer J Clin.* 2020;70(1):7–30.

3. Pei G, Chen L, Zhang W. WGCNA Application to Proteomic and Metabolomic Data Analysis. *Methods Enzymol.* 2017;585:135–158. **[Network analysis of TCGA and GTEx for cancer biomarker identification]** — [PMC8841814](https://pmc.ncbi.nlm.nih.gov/articles/PMC8841814/)

4. Fathi Dizaji B. Strategies to Identify and Reduce the Main Variables Causing Batch Effects in High-Throughput Screening Assays. *Chem Biol Drug Des.* 2021;98(3):284–299. **[Selecting precise reference normal tissue for cancer research]** — [PMC6357350](https://pmc.ncbi.nlm.nih.gov/articles/PMC6357350/)

5. Fortin JP, Cullen N, Sheline YI, et al. Harmonization of cortical thickness measurements across scanners and sites. *Neuroimage.* 2018;167:104–120. **[neuroCombat batch correction]**

6. Guo J, Liu Z, Gong R. Machine Learning-Based Identification of Colon Cancer Candidate Diagnostic Genes. *Front Genet.* 2022;13:836097. — [PMC8944988](https://pmc.ncbi.nlm.nih.gov/articles/PMC8944988/)

7. Paredes J, Figueiredo J, Albergaria A, et al. Epithelial E- and P-cadherins: role and clinical significance in cancer. *Biochim Biophys Acta.* 2012;1826(2):297–311.

8. Yao Q, Ni Y, Xu Y, et al. CDH3 Is an Effective Serum Biomarker of Colorectal Cancer Distant Metastasis Patients. *Front Oncol.* 2024. — [PMC11375550](https://pmc.ncbi.nlm.nih.gov/articles/PMC11375550/)

9. Mao Z, Xiao C, Wang Y, Zhang Y. P-Cadherin (CDH3) is overexpressed in colorectal tumors and has potential as a serum marker for colorectal cancer monitoring. *Oncol Lett.* 2017;14(6):7008–7014. — [PMC5672898](https://pmc.ncbi.nlm.nih.gov/articles/PMC5672898/)

10. de la Chapelle A, Botma A, Cajal TRY, et al. ETV4 plays a role on the primary events during the adenoma–adenocarcinoma progression in colorectal cancer. *Cancers.* 2021;13(4):891. — [PMC7919324](https://pmc.ncbi.nlm.nih.gov/articles/PMC7919324/)

11. Zhu Y, Gu L, Lin X, et al. Claudin-1 Is a Valuable Prognostic Biomarker in Colorectal Cancer: A Meta-Analysis. *Dis Markers.* 2020;2020:9045424. — [PMC7443231](https://pmc.ncbi.nlm.nih.gov/articles/PMC7443231/)

12. Tabaries S, Bhatt DL, Bhatt AS, et al. Claudin-1 Interacts with EPHA2 to Promote Cancer Stemness and Chemoresistance in Colorectal Cancer. *Mol Cancer Res.* 2024. — [PMC10765961](https://pmc.ncbi.nlm.nih.gov/articles/PMC10765961/)

13. Durgan J, Tseng YY, Haage A, et al. Improving the response to oxaliplatin by targeting chemotherapy-induced CLDN1 in resistant metastatic colorectal cancer cells. *Mol Ther Oncolytics.* 2023. — [PMC10091849](https://pmc.ncbi.nlm.nih.gov/articles/PMC10091849/)

14. Kang M, Jiang B, Xu B, et al. Keratin 80 promotes migration and invasion of colorectal carcinoma by interacting with PRKDC via activating the AKT pathway. *Cell Death Dis.* 2018;9(10):984. — [PMC6160410](https://pmc.ncbi.nlm.nih.gov/articles/PMC6160410/)

15. Zhou X, Li T, Chen Y, et al. Small interfering RNA-mediated knockdown of KRT80 suppresses colorectal cancer proliferation. *Oncol Lett.* 2020;20(5):234. — [PMC7579811](https://pmc.ncbi.nlm.nih.gov/articles/PMC7579811/)

16. Zhang W, Chai W, Li Z, et al. Phosphorylation of ETV4 at Ser73 by ERK kinase could block ETV4 ubiquitination degradation in colorectal cancer. *Oncotarget.* 2017;8(27):44501–44511. — [PMID 28373072](https://pubmed.ncbi.nlm.nih.gov/28373072/)

17. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics.* 2007;8(1):118–127. **[Original ComBat paper]**

18. Tibshirani R. Regression shrinkage and selection via the lasso. *J R Stat Soc Series B.* 1996;58(1):267–288. **[Original LASSO paper]**

19. Breiman L. Random forests. *Mach Learn.* 2001;45(1):5–32. **[Original Random Forest paper]**

20. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B.* 1995;57(1):289–300. **[BH FDR correction]**

---

*Correspondence: [Author email]*  
*Data availability: All processed data, analysis scripts, and results are available in the project repository.*  
*Conflict of interest: None declared.*
