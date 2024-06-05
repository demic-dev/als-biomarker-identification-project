# **ALS Biomarker Identification Project**

Welcome to the ALS Biomarker Identification Project repository. This project, conducted in collaboration with [Aliaksei](https://github.com/aliakseibrown) at Paris Saclay University, aims to identify potential biomarkers associated with Amyotrophic Lateral Sclerosis (ALS) through RNA-Seq data analysis.

## **Objective**

The primary objective is to use advanced computational techniques to identify genes with significant expression differences between ALS patients and a Non-Neurological control group. By identifying these biomarkers, we aim to contribute to early diagnosis and targeted interventions for ALS.

## **Key Highlights**

- **Data Preprocessing**: Parsed RNA-Seq data into DataFrame, split samples, and organized data using an object-oriented approach.
- **Descriptive Analysis**: Analyzed and visualized the distribution of gene samples, highlighting differences between ALS and control groups.
- **Dimensionality Reduction Techniques**: Used Principal Component Analysis (PCA) and t-distributed Stochastic Neighbor Embedding (tSNE) to visualize complex patterns and clusters.
- **Advanced Analysis Techniques**: Conducted univariate analyses and used the PyDESeq2 library to identify genes with significant expression differences.
- **Model Tuning and Generalization**: Applied normalization and ElasticNet model tuning to accurately identify ALS-impacted genes.

## **Getting Started**

To get started with this project:

1. Clone this repository to your local machine.
```shell
git clone https://github.com/demic-dev/als-biomarker-identification-project.git
```
2. Create a conda environment after opening the repository
```
conda create -n bio python=3.9 -y
conda activate bio
```
3. Install the necessary dependencies listed in **`requirements.txt`**.
```shell
pip install -r requirements.txt
```
4. Explore the Jupyter notebooks in the **`analysis`** directory to understand our methodology and findings.

## **Acknowledgments**

We thank Paris-Saclay University for their resources and support. We also appreciate the guidance and expertise of our mentor and collaborators.

## **Sources**
The data was provided by a research Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: Retrotransposon Activation, Oxidative Stress, and Activated Glia

- Data set [add link]


