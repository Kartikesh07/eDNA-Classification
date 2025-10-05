import argparse
import subprocess
import os
import pandas as pd
import joblib
from Bio import SeqIO

# A helper function to run external commands (like cutadapt, vsearch)
def run_command(command):
    """Runs a command in the shell and checks for errors."""
    print(f"--- Running Command ---\n{command}\n-----------------------")
    try:
        # Using shell=True for simplicity on Windows, check=True to raise an error if the command fails
        subprocess.run(command, shell=True, check=True)
        print("--- Command successful ---\n")
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        # Exit the script if a command fails
        raise

# A helper function for our k-mer feature engineering
def get_kmer_features(sequence, k=6):
    """Converts a DNA sequence into a dictionary of its k-mer counts."""
    kmer_counts = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts

def run_cleaning(args):
    """Stage 2: Runs cutadapt to clean the raw FASTQ files."""
    print(">>> STAGE 2: CLEANING RAW READS <<<")
    
    # Define the output file paths
    args.trimmed_f = os.path.join(args.temp_dir, "reads_1.trimmed.fastq")
    args.trimmed_r = os.path.join(args.temp_dir, "reads_2.trimmed.fastq")
    
    # --- FIX STARTS HERE ---
    # Build the cutadapt command as a clean, single-line string
    command = (
        f"cutadapt -q 20 -m 150 "
        f"-o {args.trimmed_f} -p {args.trimmed_r} "
        f"{args.forward_reads} {args.reverse_reads}"
    )
    # --- FIX ENDS HERE ---
    
    run_command(command)

def run_clustering(args):
    """Stage 3: Runs vsearch to discover OTUs."""
    print(">>> STAGE 3: DISCOVERING OTUs VIA CLUSTERING <<<")
    
    # Define intermediate file paths
    combined_fastq = os.path.join(args.temp_dir, "all.trimmed.fastq")
    unique_fasta = os.path.join(args.temp_dir, "unique_sequences.fasta")
    args.otus_fasta = os.path.join(args.temp_dir, "otus.fasta") # This is the final output of this stage

    # --- FIX STARTS HERE ---
    # 1. Combine trimmed files using Python, avoiding shell commands.
    print("--- Combining trimmed FASTQ files... ---")
    with open(combined_fastq, 'wb') as outfile:
        with open(args.trimmed_f, 'rb') as infile:
            outfile.write(infile.read())
        with open(args.trimmed_r, 'rb') as infile:
            outfile.write(infile.read())
    print("--- File combination successful ---\n")
    # --- FIX ENDS HERE ---

    # 2. Dereplicate (This command is universal and doesn't change)
    derep_command = f'vsearch --fastx_uniques {combined_fastq} --sizeout --fastaout {unique_fasta}'
    run_command(derep_command)

    # 3. Cluster (This command is universal and doesn't change)
    cluster_command = f'vsearch --cluster_size {unique_fasta} --id 0.97 --centroids {args.otus_fasta} --uc {os.path.join(args.temp_dir, "clusters.uc")}'
    run_command(cluster_command)

def run_classification(args):
    """Stage 4 & 5: Loads the AI model and classifies the discovered OTUs."""
    print(">>> STAGE 4 & 5: CLASSIFYING OTUs AND GENERATING REPORT <<<")

    # 1. Load the trained model assets
    print("--- Loading trained model and assets... ---")
    model = joblib.load(os.path.join(args.model_dir, 'tax_classifier.joblib'))
    vectorizer = joblib.load(os.path.join(args.model_dir, 'kmer_vectorizer.joblib'))
    label_encoder = joblib.load(os.path.join(args.model_dir, 'label_encoder.joblib'))
    
    # 2. Load the newly discovered OTUs
    otu_data = []
    for record in SeqIO.parse(args.otus_fasta, "fasta"):
        header_parts = record.id.split(';size=')
        otu_id = header_parts[0]
        abundance = int(header_parts[1])
        sequence = str(record.seq)
        otu_data.append([otu_id, abundance, sequence])
    
    otu_df = pd.DataFrame(otu_data, columns=['otu_id', 'abundance', 'sequence'])
    print(f"--- Loaded {len(otu_df)} OTUs for classification. ---")

    # 3. Featurize the OTUs
    otu_df['kmer_counts'] = otu_df['sequence'].apply(lambda seq: get_kmer_features(seq, k=6))
    X_otus = vectorizer.transform(otu_df['kmer_counts'])

    # 4. Predict and decode
    predictions_encoded = model.predict(X_otus)
    predictions_decoded = label_encoder.inverse_transform(predictions_encoded)
    otu_df['predicted_class'] = predictions_decoded

    # 5. Save the final reports
    # Full report with every OTU
    full_report_path = os.path.join(args.output_dir, f"{args.sample_name}_full_report.csv")
    otu_df.to_csv(full_report_path, index=False)
    
    # Summary report
    summary_report_path = os.path.join(args.output_dir, f"{args.sample_name}_summary_report.csv")
    biodiversity_summary = otu_df.groupby('predicted_class')['abundance'].sum().sort_values(ascending=False)
    biodiversity_summary.to_csv(summary_report_path)

    print(f"--- Final reports saved to {args.output_dir} ---")
    print("\n--- Biodiversity Summary ---")
    print(biodiversity_summary)

def main():
    # This is the main controller of our script.
    
    # Setup argument parser
    parser = argparse.ArgumentParser(description="AI-driven eDNA Classification Pipeline")
    parser.add_argument("--forward_reads", required=True, help="Path to the forward FASTQ file (R1).")
    parser.add_argument("--reverse_reads", required=True, help="Path to the reverse FASTQ file (R2).")
    parser.add_argument("--output_dir", required=True, help="Directory to save the final reports.")
    parser.add_argument("--sample_name", required=True, help="A name for the sample, used for output files.")
    parser.add_argument("--model_dir", default="models", help="Directory where the trained model assets are stored.")
    
    args = parser.parse_args()

    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    args.temp_dir = os.path.join(args.output_dir, "temp_files") # For intermediate files
    os.makedirs(args.temp_dir, exist_ok=True)

    print(f"=== Starting Pipeline for Sample: {args.sample_name} ===")

    # Run the pipeline stages in order
    run_cleaning(args)
    run_clustering(args)
    run_classification(args)

    print(f"\n=== Pipeline for Sample: {args.sample_name} Completed Successfully! ===")


# This standard Python construct makes the script runnable from the command line
if __name__ == "__main__":
    main()