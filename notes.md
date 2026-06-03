# Notes

This file contains study notes related to Python, GitHub, Kaggle, statistics, and biomedical data analysis for this project.

## GitHub Notes

* A repository is a project folder that stores files, code, documentation, and version history.
* A commit saves a version of the project at a specific point in time.
* `README.md` is the main introduction page for the project.
* `work_log.md` is used to document time spent and project progress.
* `project_plan.md` is used to organize the research question, data sources, methods, and timeline.

## Python Notes

* Python will be used for data cleaning, analysis, visualization, and machine learning.
* `pandas` is commonly used to work with table-like data.
* `NumPy` is commonly used for numerical calculations.
* `matplotlib` is commonly used to create graphs and figures.
* `scikit-learn` is commonly used for machine learning models and feature selection.

## Research Notes

* Colon adenocarcinoma is commonly abbreviated as COAD.
* Tumor-adjacent normal tissue may look normal under a microscope, but it may not be molecularly identical to healthy normal tissue.
* This project focuses on whether biomarker candidates selected by machine learning remain stable when the normal control group changes.
* The planned comparison is between tumor-adjacent normal tissue and healthy GTEx colon tissue.

## Data Source Notes

* TCGA can provide tumor and tumor-adjacent normal samples.
* GTEx can provide healthy normal tissue samples.
* UCSC Xena can be used to access and download TCGA and GTEx gene expression data.

## Questions to Explore

* How different are tumor-adjacent normal tissue and healthy normal tissue at the gene expression level?
* Do machine learning-selected biomarker candidates change when the normal control group changes?
* Which genes remain stable across different control group comparisons?
* How should batch effects between TCGA and GTEx be handled?

