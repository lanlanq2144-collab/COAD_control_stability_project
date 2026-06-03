Work Log
All dates and times are recorded in Korea Standard Time (KST).
This file documents the time I spend learning and building this project, including Python, GitHub, Kaggle, statistics, Claude Code, ChatGPT, and research-related biomedical data analysis.

May 31, 2026
Time spent: 5 hours
What I worked on:
I wanted to find a research topic I could do independently using only public data, without a lab or experiments. My starting question was whether the definition of "normal" affects what an AI model identifies as cancer-related genes.
Most of this session was spent trying to understand basic terminology. Words like gene expression, tumor-adjacent normal tissue, RNA-seq, and biomarker were all new to me. I read through five related papers, not because I fully understood them, but to check whether my question was something researchers actually study. I confirmed that the comparison between NAT and GTEx as reference controls is a real methodological debate.
What I actually understood by the end of this session:
The core research question: does changing the definition of "normal" change which genes an AI model flags as cancer-related?
The basic structure of the comparison: fix the tumor group, change the normal group, see what happens.
That public gene expression data exists and can be downloaded from UCSC Xena.
What I did not yet understand:
What TCGA and GTEx actually were. I assumed TCGA meant cancer samples and GTEx meant healthy samples as general categories. I did not know they were names of separate research projects.
What the data actually looked like or how to work with it.
Honest note:
I did not write any code today. The goal was just to find a direction worth pursuing.
Next step:
Download data from UCSC Xena and try to extract COAD-relevant samples.

June 1, 2026
Time spent: 3.5 hours
What I worked on:
Claude was not yet available to me, so I used ChatGPT as a coding assistant. I did not write code myself. I asked ChatGPT what to do and followed the instructions step by step.
I downloaded data from UCSC Xena and ran Python code to extract COAD-related samples and generate a basic scatter plot. The plot showed points grouped in different clusters, but I could not explain what it meant. I could see the groups were separated, but I did not know the term batch effect or understand why that separation was happening.
What I actually understood by the end of this session:
That combining data from different institutions might cause problems. I did not know the word for it yet, but I noticed something looked off.
That the COAD-relevant subset was much smaller than the full dataset.
What I did not yet understand:
What TCGA and GTEx were. I still thought they were just category names for cancer and healthy samples, not separate research projects.
What batch effect was or why it mattered.
What the plot was actually showing.
Honest note:
Most of this session was following instructions I did not fully understand. I had a result but could not explain it.
Next step:
Figure out what the plot is showing and whether the data can be used.

June 3, 2026
Time spent: 7 hours (6 hours project, 1 hour Kaggle)
What I worked on:
I subscribed to Claude Pro and used it to understand what I had built so far. A significant portion of this session was spent on concepts I had been missing, not on running new code.
Things I had wrong that I corrected today:
I thought TCGA meant cancer samples and GTEx meant healthy samples as general labels. I learned today that they are names of two completely separate large-scale data collection projects run by different organizations. TCGA collected data from cancer patients, GTEx collected data from healthy tissue donors including postmortem samples.
I did not understand what the scatter plot from June 1 was showing. I learned today that the group separation I saw was called a batch effect, meaning the groups were separating not only because of biological differences but because the data came from different institutions using different equipment and protocols.
I did not know my own research question clearly. I had been running code without fully understanding what I was trying to find. Today I was able to say it in one sentence for the first time: does changing the definition of normal change which genes a machine learning model selects as cancer biomarker candidates?
Concepts I understood for the first time today:
PCA: a method to compress tens of thousands of gene expression values into a 2D plot so you can see patterns. Each point is one sample, not one gene.
Volcano plot: each point is one gene, not one sample. The x-axis shows how different a gene is between cancer and normal. The y-axis shows how statistically certain that difference is.
DEG: a list of genes that are expressed differently between cancer and normal tissue. Each row in the output file is one gene, with columns showing the fold change, p-value, and adjusted p-value.
Batch effect: a systematic difference caused by measurement environment, not biology. The separation I saw in my June 1 plot was partly this.
ComBat: a statistical method to remove batch effects. It does not work perfectly and residual artifacts can remain.
The difference between deep learning and traditional machine learning. LASSO and Random Forest are math-based models that run in minutes on a laptop. They are not the same as the kind of AI used to train large language models.
AUC: a measure of how well a model separates two groups. 0.5 means random chance, 1.0 means perfect separation. All four comparisons returned AUC 1.000, which I initially found suspicious. I learned this happens because the input features were already the most extreme genes from DEG analysis. However, there is also a risk called feature-selection leakage that I need to flag as a limitation.
What happened during the analysis today:
Claude Code ran the full pipeline including DEG analysis, machine learning, biomarker stability comparison, and biological validation. I gave instructions and asked questions about the results as they came out.
During the biological validation step, I noticed something suspicious. Two of the six final candidate genes showed almost no difference compared to NAT but extreme differences compared to GTEx. I asked about this discrepancy. One gene turned out to be a retired Ensembl ID that no longer exists. The other was a pseudogene with no known biological function. Both were removed from the final candidate list. The four remaining genes were CDH3, CLDN1, ETV4, and KRT80, all of which have prior literature support in colorectal cancer research.
I did not identify these problems independently. I noticed the pattern in the numbers and asked a question about it. Claude confirmed the explanation.
What I still do not fully understand:
The detailed mathematics behind LASSO, Random Forest, and ComBat.
How to interpret pathway analysis results beyond the basic level.
Why survival analysis did not show significant results and what multivariate Cox regression would add.
Whether the AUC 1.000 results reflect a real biological signal or are inflated by the analysis design.
Technical problems encountered:
Forgot to run cd before commands in terminal every time I opened a new window, causing repeated file not found errors.
Pasted Python code directly into terminal instead of saving it to a file first.
pycombat threw an ImportError due to lowercase spelling. Required capital C in Combat. After fixing that, got a dimension mismatch error. Switched to neuroCombat which worked correctly.
Kaggle:
Completed 2 out of 5 chapters of Intro to Programming.
Output:
6 QC plots, 4 volcano plots, ROC curves, feature importance plots, Venn diagrams, final 4 biomarker candidates, full draft report.
Next steps:
Fill in empty tables in the report. Add limitation around AUC and ComBat. Fix reference formatting. Push to GitHub. Continue building understanding of each analysis step through Kaggle and further reading.
