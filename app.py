# --- This is the FINAL version of app.py ---

import streamlit as st
import subprocess
import os
import time
import sys
import pandas as pd

st.set_page_config(page_title="JalNetra AI eDNA Classifier", page_icon="🧬", layout="wide")

st.title("🧬 JalNetra: AI-Powered eDNA Classification")
st.write("Upload paired-end 18S rRNA FASTQ files to discover and classify eukaryotic biodiversity.")

st.header("1. Upload Your Data")
uploaded_file_r1 = st.file_uploader("Upload Forward Read FASTQ File (R1)", type=["fastq", "fq", "gz"])
uploaded_file_r2 = st.file_uploader("Upload Reverse Read FASTQ File (R2)", type=["fastq", "fq", "gz"])
sample_name_input = st.text_input("Enter a name for your sample", "my_edna_sample")

st.header("2. Run the Analysis")
run_button = st.button("Start Pipeline", type="primary")

if run_button:
    if uploaded_file_r1 and uploaded_file_r2 and sample_name_input:
        run_timestamp = int(time.time())
        temp_dir = f"temp_run_{run_timestamp}"
        os.makedirs(temp_dir, exist_ok=True)
        
        forward_reads_path = os.path.join(temp_dir, uploaded_file_r1.name)
        reverse_reads_path = os.path.join(temp_dir, uploaded_file_r2.name)
        with open(forward_reads_path, "wb") as f:
            f.write(uploaded_file_r1.getbuffer())
        with open(reverse_reads_path, "wb") as f:
            f.write(uploaded_file_r2.getbuffer())

        output_dir = os.path.join(temp_dir, "results")
        
        st.write("Pipeline started... This may take several minutes.")
        
        with st.spinner('Processing... Logs will appear below.'):
            log_box = st.empty()
            command = [
                sys.executable, "run_pipeline.py",
                "--forward_reads", forward_reads_path,
                "--reverse_reads", reverse_reads_path,
                "--output_dir", output_dir,
                "--sample_name", sample_name_input
            ]
            
            try:
                process = subprocess.run(command, capture_output=True, text=True, check=True, encoding='utf-8')
                log_box.text_area("Pipeline Log", process.stdout + process.stderr, height=300)
                st.success("Pipeline Completed Successfully!")
                
                st.header("3. Results")
                summary_file_path = os.path.join(output_dir, f"{sample_name_input}_summary_report.csv")
                
                if os.path.exists(summary_file_path):
                    # Only import pandas here, when we absolutely need it to show results
                    summary_df = pd.read_csv(summary_file_path)
                    st.subheader("Biodiversity Summary")
                    st.dataframe(summary_df)
                    
                    st.subheader("Community Composition Chart")
                    st.bar_chart(summary_df.set_index('predicted_class'))
                else:
                    st.error("Result file not found. Check the pipeline log for errors.")

            except subprocess.CalledProcessError as e:
                st.error("The pipeline failed! See the log for details.")
                log_box.text_area("Error Log", e.stdout + e.stderr, height=300)

    else:
        st.warning("Please upload both FASTQ files and provide a sample name.")