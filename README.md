# JalNetra: An AI-Powered eDNA Classification Pipeline

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://jalnetra.streamlit.app/) 
*(Note: The live deployment is currently unstable due to challenges with external binary dependencies. Please see the "Deployment Status" section below.)*

JalNetra ("The Eye of the Water") is a fully automated, AI-driven pipeline designed to identify and quantify eukaryotic biodiversity from raw environmental DNA (eDNA) sequences. This project was developed for the Smart India Hackathon to address the critical challenge of classifying unknown species from unique environments like the deep sea, where traditional reference databases fall short.

## The Problem: The "Unassigned" Majority in eDNA Analysis

Traditional bioinformatics pipelines rely on direct alignment of eDNA sequences to reference databases (like BLAST). This approach fails in biodiverse and under-studied environments, leading to:
1.  **Massive Data Loss:** A high percentage of sequences from novel or divergent species are labeled "Unassigned" or "Unknown," causing a significant underestimation of true biodiversity.
2.  **Slow Performance:** Aligning tens of thousands of sequences against massive databases is computationally expensive and can take hours or days for a single sample.
3.  **Misclassification:** Novel species are often incorrectly assigned to their closest, well-known (but biologically different) relatives.

## Our Solution: A Discovery-First, AI-Driven Approach

JalNetra flips the traditional model on its head. Instead of matching first, we **discover first**.

 
*(Suggestion: Take a screenshot of the "Technical Workflow" section of your PPT and upload it to a site like Imgur, then paste the link here.)*

Our pipeline works in two main stages:
1.  **Unsupervised Discovery:** The pipeline first uses `VSEARCH`, a high-performance bioinformatics tool, to perform dereplication and clustering. It groups all similar DNA sequences into **Operational Taxonomic Units (OTUs)** at 97% identity. This step is completely database-free and allows us to discover every potential species signature in the sample.
2.  **AI-Powered Classification:** We then use a `Logistic Regression` model, trained on the comprehensive PR2 database, to classify these discovered OTUs. The model learned the fundamental k-mer patterns of major eukaryotic groups and achieved **99.6% accuracy** during validation. It can provide a high-confidence classification (e.g., `Ascomycota`, `Arthropoda`) even for OTUs that don't have a perfect match in any database.

This entire engine is wrapped in a user-friendly Streamlit web application, making advanced genomic analysis accessible to any researcher.

## Tech Stack

*   **Backend Engine:** Python, Pandas, Biopython, Scikit-learn
*   **Command-Line Tools:** `cutadapt`, `VSEARCH`
*   **Frontend UI:** Streamlit
*   **Environment Management:** Conda

## How the Prototype Works

The repository contains two main components:
1.  **`run_pipeline.py`:** The master command-line script that automates the entire workflow:
    *   **Input:** Paired-end raw FASTQ files.
    *   **Process:** Calls `cutadapt` for cleaning, then `VSEARCH` for clustering, and finally uses the pre-trained `scikit-learn` model to classify the resulting OTUs.
    *   **Output:** A final biodiversity summary in CSV format.
2.  **`app.py`:** The Streamlit web application.
    *   Provides a simple UI for users to upload their FASTQ files.
    *   Calls the `run_pipeline.py` script as a subprocess.
    *   Reads the final CSV report and displays the results as an interactive table and bar chart.

## Deployment Status and Known Challenges

The application has been deployed to Streamlit Community Cloud for demonstration purposes.

**URL:** [https://jalnetra.streamlit.app/](https://jalnetra.streamlit.app/)

**Current Status: Unstable.**

The primary challenge for a stable cloud deployment lies in the pipeline's dependency on external, pre-compiled command-line tools, specifically **`VSEARCH`**. While the Python environment and its dependencies are perfectly managed by `requirements.txt`, ensuring a consistent and memory-efficient execution of the `VSEARCH` binary within the resource-limited free tier of Streamlit Cloud has proven difficult.

The application frequently crashes with a `503 Service Unavailable` error, which indicates the backend process is being terminated by the host, likely due to exceeding memory or CPU limits during the intensive clustering stage.

**Future Work for Deployment:** A more robust deployment would involve:
*   Containerizing the entire application with **Docker** to create a fully self-contained environment.
*   Deploying the container on a more flexible cloud service (like AWS Fargate, Google Cloud Run, or Heroku) that provides more control over computational resources.