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

**Time spent:** 4 hours

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
