# --- This is the FINAL version of run_pipeline.py ---

import argparse
import subprocess
import os
import pandas as pd
import joblib
from Bio import SeqIO
import sys

def run_command(command):
    print(f"--- Running Command ---\n{command}\n-----------------------")
    try:
        subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        print("--- Command successful ---\n")
    except subprocess.CalledProcessError as e:
        print(f"Error stdout: {e.stdout}")
        print(f"Error stderr: {e.stderr}")
        raise

def get_kmer_features(sequence, k=6):
    kmer_counts = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts

def run_cleaning(args):
    print(">>> STAGE 2: CLEANING RAW READS <<<")
    args.trimmed_f = os.path.join(args.temp_dir, "reads_1.trimmed.fastq")
    args.trimmed_r = os.path.join(args.temp_dir, "reads_2.trimmed.fastq")
    
    cutadapt_command = (
        f"{sys.executable} -m cutadapt -q 20 -m 150 "
        f"-o {args.trimmed_f} -p {args.trimmed_r} "
        f"{args.forward_reads} {args.reverse_reads}"
    )
    run_command(cutadapt_command)

def run_clustering(args):
    print(">>> STAGE 3: DISCOVERING OTUs VIA CLUSTERING <<<")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    vsearch_path = os.path.join(script_dir, "bin", "vsearch")
    
    if os.path.exists(vsearch_path):
        os.chmod(vsearch_path, 0o755)
    else:
        raise FileNotFoundError(f"VSEARCH executable not found at {vsearch_path}")

    combined_fastq = os.path.join(args.temp_dir, "all.trimmed.fastq")
    unique_fasta = os.path.join(args.temp_dir, "unique_sequences.fasta")
    args.otus_fasta = os.path.join(args.temp_dir, "otus.fasta")

    print("--- Combining trimmed FASTQ files... ---")
    with open(combined_fastq, 'wb') as outfile:
        with open(args.trimmed_f, 'rb') as infile:
            outfile.write(infile.read())
        with open(args.trimmed_r, 'rb') as infile:
            outfile.write(infile.read())
    print("--- File combination successful ---\n")

    derep_command = f'{vsearch_path} --fastx_uniques {combined_fastq} --sizeout --fastaout {unique_fasta}'
    run_command(derep_command)

    cluster_command = f'{vsearch_path} --cluster_size {unique_fasta} --id 0.97 --centroids {args.otus_fasta} --uc {os.path.join(args.temp_dir, "clusters.uc")}'
    run_command(cluster_command)

def run_classification(args):
    print(">>> STAGE 4 & 5: CLASSIFYING OTUs AND GENERATING REPORT <<<")
    model = joblib.load(os.path.join(args.model_dir, 'tax_classifier.joblib'))
    vectorizer = joblib.load(os.path.join(args.model_dir, 'kmer_vectorizer.joblib'))
    label_encoder = joblib.load(os.path.join(args.model_dir, 'label_encoder.joblib'))
    
    otu_data = []
    for record in SeqIO.parse(args.otus_fasta, "fasta"):
        header_parts = record.id.split(';size=')
        otu_id = header_parts[0]
        abundance = int(header_parts[1])
        sequence = str(record.seq)
        otu_data.append([otu_id, abundance, sequence])
    
    otu_df = pd.DataFrame(otu_data, columns=['otu_id', 'abundance', 'sequence'])
    print(f"--- Loaded {len(otu_df)} OTUs for classification. ---")

    otu_df['kmer_counts'] = otu_df['sequence'].apply(lambda seq: get_kmer_features(seq, k=6))
    X_otus = vectorizer.transform(otu_df['kmer_counts'])

    predictions_encoded = model.predict(X_otus)
    predictions_decoded = label_encoder.inverse_transform(predictions_encoded)
    otu_df['predicted_class'] = predictions_decoded

    summary_report_path = os.path.join(args.output_dir, f"{args.sample_name}_summary_report.csv")
    biodiversity_summary = otu_df.groupby('predicted_class')['abundance'].sum().sort_values(ascending=False)
    biodiversity_summary.to_csv(summary_report_path)

    print(f"--- Final reports saved to {args.output_dir} ---")

def main():
    parser = argparse.ArgumentParser(description="AI-driven eDNA Classification Pipeline")
    parser.add_argument("--forward_reads", required=True)
    parser.add_argument("--reverse_reads", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--sample_name", required=True)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_model_dir = os.path.join(script_dir, "models")
    parser.add_argument("--model_dir", default=default_model_dir)
    
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    args.temp_dir = os.path.join(args.output_dir, "temp_files")
    os.makedirs(args.temp_dir, exist_ok=True)

    print(f"=== Starting Pipeline for Sample: {args.sample_name} ===")
    run_cleaning(args)
    run_clustering(args)
    run_classification(args)
    print(f"\n=== Pipeline for Sample: {args.sample_name} Completed Successfully! ===")

if __name__ == "__main__":
    main()