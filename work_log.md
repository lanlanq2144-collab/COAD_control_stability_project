# Work Log

> All dates and times are recorded in Korea Standard Time (KST).
>
> This file documents the time I spend learning and building this project, including Python, GitHub, Kaggle, statistics, Claude Code, ChatGPT, and research-related biomedical data analysis.

---

## May 31, 2026 — Research Topic Selection

**Time spent:** 5 hours

### What I worked on

I wanted to find a research topic I could do independently using only public data, without a lab or experiments. My starting question was whether the definition of "normal" affects what an AI model identifies as cancer-related genes.

Most of this session was spent trying to understand basic terminology. Words like gene expression, tumor-adjacent normal tissue, RNA-seq, and biomarker were all new to me. I read through five related papers, not because I fully understood them, but to check whether my question was something researchers actually study. I confirmed that the comparison between NAT and GTEx as reference controls is a real methodological debate.

### ✅ What I actually understood

- **Core research question:** Does changing the definition of "normal" change which genes an AI model selects as cancer-related?
- **Basic structure of the comparison:** Fix the tumor group, change the normal group, see what happens
- Public gene expression data exists and can be downloaded from UCSC Xena

### ❌ What I did not yet understand

- What TCGA and GTEx actually were — I assumed they were general category names for cancer and healthy samples, not names of separate research projects
- What the data actually looked like or how to work with it

### 📝 Honest note

I did not write any code today. The goal was just to find a direction worth pursuing.

---

## June 1, 2026 — First Data Download

**Time spent:** 3.5 hours

### What I worked on

Claude was not yet available to me, so I used ChatGPT as a coding assistant. I did not write code myself. I asked ChatGPT what to do and followed the instructions step by step.

I downloaded data from UCSC Xena and ran Python code to extract COAD-related samples and generate a basic scatter plot. The plot showed points grouped in different clusters, but I could not explain what it meant. I could see the groups were separated, but I did not know the term batch effect or understand why that separation was happening.

### ✅ What I actually understood

- Combining data from different institutions might cause problems — I did not know the word for it yet, but something looked off
- The COAD-relevant subset was much smaller than the full dataset

### ❌ What I did not yet understand

- What TCGA and GTEx were — I still thought they were category names, not separate research projects
- What batch effect was or why it mattered
- What the plot was actually showing

### 📝 Honest note

Most of this session was following instructions I did not fully understand. I had a result but could not explain it.

---

## June 3, 2026 — Full Analysis Pipeline Completed

**Time spent:** 7 hours (6 hours project, 1 hour Kaggle)

### What I worked on

I subscribed to Claude Pro and used it for the first time. A significant portion of this session was spent on concepts I had been missing, not on running new code.

### 🔄 Things I had wrong that I corrected today

| What I believed | What is actually true |
|---|---|
| TCGA = category name for cancer samples, GTEx = category name for healthy samples | Two completely separate large-scale data collection projects run by different organizations |
| The June 1 plot was showing cancer vs normal | The separation was partly due to batch effect — data from different institutions using different equipment |
| I understood my own research | I had been running code without understanding what I was trying to find |

### ✅ Concepts I understood for the first time today

- **PCA:** Compresses tens of thousands of gene expression values into a 2D plot. Each point is one sample, not one gene
- **Volcano plot:** Each point is one gene, not one sample. X-axis shows how different a gene is between cancer and normal. Y-axis shows how statistically certain that difference is
- **DEG:** A list of genes expressed differently between cancer and normal tissue
- **Batch effect:** A systematic difference caused by measurement environment, not biology
- **ComBat:** A statistical method to remove batch effects. Not perfect — residual artifacts can remain
- **Deep learning vs traditional ML:** LASSO and Random Forest are math-based models that run in minutes on a laptop. Not the same as training large language models
- **Python basics (from Kaggle):** print statements and basic arithmetic operations — understood. Functions — covered but still a bit confusing, need to review

### 🔬 Analysis completed today

I installed and used Claude Code for the first time. I did not know that Homebrew → Node → Claude Code were three separate steps required before it would work.

Claude Code ran the full pipeline below. My role was to give instructions and ask questions about the results as they came out.

| Step | Result |
|---|---|
| QC (PCA + ComBat batch correction) | Completed |
| DEG analysis (4 control groups) | Tumor vs NAT: 11,623 / GTEx All: 17,562 / Transverse: 16,516 / Sigmoid: 19,298 |
| Machine learning (LASSO + Random Forest) | Completed |
| Biomarker stability (Venn diagrams) | NAT vs GTEx overlap: 1–2 genes / GTEx groups overlap: 5–8 genes |
| Biological validation | 4 final biomarker candidates confirmed |

**Finding during biological validation:**

Two of the six final candidates showed almost no difference compared to NAT but extreme differences compared to GTEx. I noticed this pattern in the numbers and asked about it. One turned out to be a retired Ensembl gene ID that no longer exists. The other was a pseudogene with no known biological function. Both were removed.

> I did not identify these problems independently. I noticed the pattern and asked a question. Claude confirmed the explanation.

**Final 4 biomarker candidates:**

| Gene | Role | Prior literature |
|---|---|---|
| CDH3 | Cell adhesion | Established CRC biomarker |
| CLDN1 | Epithelial barrier function | Confirmed by meta-analysis |
| ETV4 | Transcription factor | Key role in CRC progression |
| KRT80 | Cytoskeleton | Promotes invasion and metastasis |

### ❌ What I still do not fully understand

- The mathematical principles behind LASSO, Random Forest, and ComBat
- How to interpret pathway analysis results beyond a basic level
- Why survival analysis did not show significant results and what multivariate Cox regression would add
- What AUC means and how to properly interpret the machine learning performance results
- What each column in the DEG output file means (log2FC, p-value, padj, mean_tumor, mean_control)

### 🐛 Technical errors encountered

- Forgot to run `cd` before commands every time I opened a new terminal window
- Pasted Python code directly into terminal instead of saving it to a file first
- pycombat threw an ImportError due to lowercase spelling — required capital C in `Combat` — after fixing that, got a dimension mismatch error — switched to neuroCombat which worked

### 📦 Output

- 6 QC plots
- 4 volcano plots
- ROC curves, feature importance plots
- 2 Venn diagrams
- 4 final biomarker candidates: CDH3, CLDN1, ETV4, KRT80
- Full draft report (Report_first.pages)

### 📝 Honest note

Claude Code ran the analysis. My role today was to give instructions, follow the results, and try to understand what was happening at each step. The pipeline is complete but my understanding is still catching up.

---

## Next Steps

- [ ] Understand what AUC means and how to interpret machine learning performance results
- [ ] Understand what each column in the DEG output file means (log2FC, p-value, padj, mean_tumor, mean_control)
- [ ] Fill in empty tables in the report
- [ ] Add AUC leakage limitation to methods
- [ ] Strengthen ComBat methodology explanation
- [ ] Fix reference formatting
- [ ] Push everything to GitHub
- [ ] Continue Kaggle — currently 2/5 chapters of Intro to Programming complete
- [ ] Review Python functions — covered today but still confusing

---

## June 4, 2026 — Kaggle Python Fundamentals

**Time spent:** 4 hours (Kaggle) + 3 hours (project thinking)

### What I worked on

Continued Kaggle courses and spent time thinking about next steps for the research project.

Completed Intro to Programming chapters 3, 4, and 5 — finished the full course. Started Python course and completed up to step 2 out of 7.

### ✅ Concepts I understood for the first time today

- **Data types:** int, float, and type() — understood the difference between number types
- **Boolean logic:** True functions as 1 and False functions as 0
- **Conditionals:** if, else, and elif — learned how to write conditional logic
- **Indentation rules:** one level of indentation inside an if block, two levels inside a nested block
- **= vs ==:** assignment vs comparison — understood the difference more clearly
- **return:** understood what it actually means in a function and how to distinguish the function definition from its body
- **Lists:** learned how to use lists to organize data, and understood the difference between storing items in a list vs storing them as a plain string
- **List methods:** append and remove
- **help() function:** learned that this can be used to look up how a function works
- **docstring:** learned that this is the description written inside a function to explain what it does
- **round(value, digits):** learned how to use this function and what each argument does

### ❌ What I did not yet understand

- Functions with a variable number of parameters — the concept of parameters in general is still confusing
- Why the number of parameters in a function can change — this is where things started breaking down today

### 📝 Honest note

Functions were starting to make sense earlier in the day but became confusing again as the material got more complex. I plan to come back to the Intro to Programming course if I keep getting stuck in the Python course. The practice exercises in Intro to Programming were also harder than expected — got stuck there too.

Still thinking through how to move the research project forward. No new analysis today.

### ⏭️ Next session goals

- [ ] Review functions and parameters until they feel solid
- [ ] Continue Python course from step 3
- [ ] Start thinking concretely about next steps for the research project

---

## June 8, 2026 — Research Understanding + GitHub

**Time spent:** 7 hours

### What I worked on

Went through the entire analysis pipeline step by step to understand what each stage actually did. Also completed external validation and pushed everything to GitHub.

### ✅ Concepts I understood for the first time today

- **QC purpose:** Confirming whether data can be trusted before analysis. PCA shows group separation; separation between same-type groups signals batch effect
- **Batch effect:** Systematic difference caused by measurement environment, not biology. NAT and GTEx separated in PCA because they came from different institutions
- **Why Tumor was excluded from ComBat:** Including Tumor would erase the real biological signal between cancer and normal
- **DEG filtering logic:** Two conditions must both be met — |log2FC| > 1 (2x difference) and FDR < 0.05 (statistically certain). FDR is stricter than p-value because 60,498 genes are tested simultaneously
- **log2FC direction:** Positive = higher in tumor, negative = higher in normal
- **Volcano plot:** Each point is one gene. X-axis = how different, Y-axis = how certain. Top right = most important candidates
- **Why 500 DEGs for ML input:** Too many features causes noise; top 500 by FDR are the most reliably different genes
- **LASSO logic:** Zeros out unimportant genes, keeps only the most predictive ones. Two modes: strict (1–3 genes) and loose (19–40 genes)
- **Random Forest logic:** 300 trees each vote; genes used most often across trees get high importance scores
- **Why LASSO ∩ RF top-50:** Two methods with different approaches both agreeing on a gene means higher confidence
- **Nested CV:** Fixes feature-selection leakage by doing DEG selection inside each training fold. AUC stayed 1.000 after correction — confirms genuine biological signal
- **Venn diagram key finding:** NAT vs GTEx overlap only 1–2 genes. GTEx groups overlap 5–8 genes. Biomarker candidates change substantially depending on how "normal" is defined — this is the core finding of the study
- **Why two artifact genes were removed:** ENSG00000279473 (retired Ensembl ID) and HMGN1P36 (pseudogene) both showed near-zero log2FC vs NAT but extreme log2FC vs GTEx — pattern consistent with batch effect artifact, not biology
- **Final 4 biomarkers:** CDH3, CLDN1, ETV4, KRT80 — all upregulated in tumor, all confirmed in prior literature, CDH3+ETV4+CLDN1 combination independently replicated by Guo et al. 2022
- **Why survival analysis showed no significant results:** Single gene median split is too simple. Survival is determined by many factors (stage, age, treatment). Cox regression with covariates would be needed
- **Liquid biopsy:** Blood or stool-based cancer detection already exists (Cologuard, Guardant Shield). Main barriers are cost, false positives, and lack of standardization

### 🔬 Analysis completed today

**External validation (GSE156451)**

Independently validated all 4 biomarker candidates in a Chinese CRC cohort (n=72 paired samples) from Wuhan University.

| Gene | log₂FC | p-value |
|---|---|---|
| CDH3 | +3.16 | 1.1×10⁻³⁶ |
| CLDN1 | +3.35 | 7.8×10⁻³⁵ |
| ETV4 | +4.08 | 2.8×10⁻⁴⁴ |
| KRT80 | +2.61 | 1.5×10⁻³¹ |

All 4 genes replicated in a completely independent dataset from a different country, different hospital, different equipment.

**GitHub**

Pushed the full project to GitHub including scripts, results, README, and WORK_LOG.

```
scripts/   → 11 analysis scripts
results/   → 83 files across all analysis stages
README.md  → project overview
WORK_LOG.md → this file
```

### ❌ What I still do not fully understand

- The mathematical principles behind LASSO, Random Forest, and ComBat
- How to interpret pathway analysis results beyond a basic level
- Why survival analysis showed no significant results and what Cox regression would add
- How to write the code myself without Claude Code

### 📝 Honest note

Most of today was understanding what the pipeline I built actually does. Claude Code ran the external validation and handled the GitHub push. My role was to ask questions, follow the logic of each step, and identify where I still had gaps. By the end of the session I could explain the full pipeline from QC to external validation in my own words.

### ⏭️ Next session goals

- [ ] Add external validation results to the report
- [ ] Fill in empty tables in the report
- [ ] Strengthen ComBat methodology explanation
- [ ] Fix reference formatting
- [ ] Continue Kaggle — Python course step 3 onward
- [ ] Review functions and parameters

---

## Project Reflections

### Honest thoughts after wrapping up this first project

**On having no mentor**

Without a supervisor or lab, I spent more time thinking about the project outside of actual work sessions than inside them. Ideas would come to me randomly — something I thought was a breakthrough — and when I looked it up, it was either already done, not feasible, or both. I felt the limits of working alone constantly. A mentor would have saved me from a lot of dead ends.

**On Python and code**

I went in thinking Python was supposed to be easy. It is not easy, at least not at the level I was trying to understand it. I tried to work toward understanding the code Claude Code produced, but I quickly realized that getting from where I am now to that level of understanding is going to take a very long time. I also realized that free resources online are not enough to get to the depth I want. That limitation became very obvious very fast.

**On biology**

I thought I had a reasonable foundation in life sciences going into this. I did not. Concepts like pseudogenes stopped me completely — I had no idea what they were. I also did not understand how much interpretation matters at every step. The method you choose at one stage completely shapes what is possible at the next. Because I do not know enough, I could have been heading in a wrong direction without any way to notice.

**On the project timeline**

I originally planned for three to four months. The further I got in, the more I understood that sustaining a project for that long, with no structure or external accountability, is genuinely difficult. It is a different kind of hard than I expected.

**On Claude Code**

The idea generation was mine. Claude Code executed those ideas quickly. But the token usage was much faster than I anticipated, and I realized I am not using the tool at anywhere near its full capability yet. There is a gap between giving it instructions and actually knowing how to direct it well.

**On whether this approach is right**

I want to know: is this actually how independent research is supposed to work? And in a real research or engineering setting, how do people actually do this kind of work day to day — how structured is it, how do they handle not knowing things, how do they know when a direction is wrong? I genuinely do not know the answer to any of that yet.
