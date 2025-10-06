import sys
import streamlit as st
import subprocess
import os
import pandas as pd
import time

# --- Page Configuration ---
st.set_page_config(
    page_title="AI eDNA Classifier",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 AI-Powered eDNA Classification Pipeline")
st.write("Upload your paired-end FASTQ files, and this pipeline will automatically process them and predict the eukaryotic biodiversity.")

# --- File Uploader ---
st.header("1. Upload Your Data")
uploaded_file_r1 = st.file_uploader("Upload Forward Read FASTQ File (R1)", type=["fastq", "fq"])
uploaded_file_r2 = st.file_uploader("Upload Reverse Read FASTQ File (R2)", type=["fastq", "fq"])
sample_name_input = st.text_input("Enter a name for your sample", "my_edna_sample")

# --- Run Pipeline Button ---
st.header("2. Run the Analysis")
run_button = st.button("Start Pipeline", type="primary")

if run_button and uploaded_file_r1 and uploaded_file_r2 and sample_name_input:
    # --- Setup temporary directories ---
    # Create a unique temporary directory for this run to avoid conflicts
    run_timestamp = int(time.time())
    temp_dir = f"temp_run_{run_timestamp}"
    os.makedirs(temp_dir, exist_ok=True)
    
    # Save uploaded files to the temporary directory
    forward_reads_path = os.path.join(temp_dir, uploaded_file_r1.name)
    reverse_reads_path = os.path.join(temp_dir, uploaded_file_r2.name)
    with open(forward_reads_path, "wb") as f:
        f.write(uploaded_file_r1.getbuffer())
    with open(reverse_reads_path, "wb") as f:
        f.write(uploaded_file_r2.getbuffer())

    output_dir = os.path.join(temp_dir, "results")
    
    # --- Run the pipeline using our run_pipeline.py script ---
    st.write("Pipeline started... This will take several minutes. Please wait.")
    
    # Display a spinner and a status box while the pipeline runs
    with st.spinner('Processing...'):
        status_box = st.empty()
        
        # Build the command
        command = [
            sys.executable, "run_pipeline.py", # <--- THIS IS THE FIX
            "--forward_reads", forward_reads_path,
            "--reverse_reads", reverse_reads_path,
            "--output_dir", output_dir,
            "--sample_name", sample_name_input
        ]
        
        try:
            # We capture the output to display it in the app
            process = subprocess.run(command, capture_output=True, text=True, check=True)
            status_box.text_area("Pipeline Log", process.stdout + process.stderr, height=300)
            st.success("Pipeline Completed Successfully!")
            
            # --- Display Results ---
            st.header("3. Results")
            summary_file_path = os.path.join(output_dir, f"{sample_name_input}_summary_report.csv")
            
            if os.path.exists(summary_file_path):
                summary_df = pd.read_csv(summary_file_path)
                st.subheader("Biodiversity Summary")
                st.dataframe(summary_df)
                
                st.subheader("Community Composition Chart")
                st.bar_chart(summary_df.set_index('predicted_class'))
            else:
                st.error("Result file not found. Check the pipeline log for errors.")

        except subprocess.CalledProcessError as e:
            st.error("The pipeline failed!")
            status_box.text_area("Error Log", e.stdout + e.stderr, height=300)

else:
    if run_button:
        st.warning("Please upload both FASTQ files and provide a sample name.")