import os
import sys
import subprocess
import argparse
from datetime import datetime

def log(message,logfile="preprocess_bulk.log"):
    # Log messages to a file with a timestamp.
    if not os.path.exists(logfile):
        with open(logfile, "w") as f:
            f.write("Log file created.\n")
    timestamp = datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    with open(logfile, "a") as f:
        f.write(f"{timestamp}: {message}\n")

def parse_args():
    # Parse command line arguments for the script.    
    parser = argparse.ArgumentParser(description="Bulk RNA-seq preprocessing pipeline.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sra", type=str, help="SRA ID to process (e.g. SRR1234567)")
    group.add_argument("--list", type=str, help="Path to file with one SRA ID per line")
    parser.add_argument("--input_dir", type=str, default="./data", help="Directory to store raw files")
    parser.add_argument("--output_dir", type=str, default="./output", help="Directory to store FastQC output")
    parser.add_argument("--GTF", type=str, help="Gene annotation of the reference genome in GTF format")
    parser.add_argument("--genome", type=str, help="Reference genome in FASTA format")
    parser.add_argument("--logfile", type=str, default="preprocess_bulk.log", help="Log file to record operations")
    return parser.parse_args()


def execute(cmd, outputfile = None, quitOnError = True):
    # Execute a shell command and log the output.
    if outputfile and os.path.exists(outputfile):
        print("Output file already exists. Skipping function.\n")
    else:
        print(f"Running {cmd} \n")
        log(f"{cmd}")
        res = subprocess.run(cmd,shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        with open(outputfile, "w") as f:
            f.write(res.stdout.decode())
            f.write(res.stderr.decode())

        if (res.returncode != 0):  #0 indicates success
            print(f"Error executing {cmd}.")
            if (quitOnError):
                sys.exit(1) 

def run_sra(sra_id, input_dir, output_dir):
    # Download SRA data and convert to FASTQ format.
    fastq_file = os.path.join(input_dir, f"{sra_id}.fastq")
    split_fastq_file = os.path.join(input_dir, f"{sra_id}_1.fastq")
    execute(f"fastq-dump {sra_id}", fastq_file)
    execute(f"fastq-dump {sra_id} --split-files", split_fastq_file)

def run_fastqc(sra_id, input_dir, output_dir):
    # Run FastQC on the downloaded FASTQ files to check quality of sequences.
    fastqc_out_dir = os.path.join(output_dir, "fastqc-out")
    os.makedirs(fastqc_out_dir, exist_ok=True)
    fastqc_input = [
        f"{sra_id}_1.fastq",
        f"{sra_id}_2.fastq"
    ]
    fastqc_exec = f"fastqc -o {fastqc_out_dir} -t4"
    fastqc_cmd = f"{fastqc_exec} {' '.join(fastqc_input)}"
    execute(fastqc_cmd, os.path.join(fastqc_out_dir, f"{sra_id}_1_fastqc.html"))

def run_star(gtf_file, genome_file, output_dir, sra_id, input_dir):
    #Indexing of reference genome using STAR.  (memory-intensive)
    run_mode = f"STAR --runMode genomeGenerate --genomeDir {output_dir} --genomeFastaFiles {genome_file} --sjdbGTFfile {gtf_file}"
    execute(run_mode) 
    # Aligning reads to the reference genome using STAR.
    align_mode = f"STAR --runThreadN 4 --genomeDir {output_dir} --readFilesIn {input_dir}/{sra_id}_1.fastq {input_dir}/{sra_id}_2.fastq --outSAMtype BAM SortedByCoordinate"
    execute(align_mode, os.path.join(output_dir, f"{sra_id}_aligned.bam"))

def run_samtools(input_bam, output_dir):
    #Use samtools to sort and index the BAM file.
    #sorted_bam = os.path.join(output_dir, f"{sra_id}_sorted.bam")
    #execute(f"samtools sort -o {sorted_bam} {input_bam}")
    pass

def run_htseq():
    #Run HTSeq to count reads per gene.
    pass  

#list of software preprocess functions stored in a dictionary
software_preprocess = {
    "sra" : run_sra,
    "fastqc": run_fastqc,
    "STAR indexing and aligner": run_star
    # "samtools": run_samtools,
    # "htseq": run_htseq
}

def select_operations():
    # Select operations to run based on user input.
    cmds = list(software_preprocess.keys()) 
    execs = [f"{i}. {cmd}" for i, cmd in enumerate(cmds, start=1)]
    print(f"Available operations: {execs}")
    response = input("Please select the operations you want to run by entering their numbers separated by commas (e.g., 1,2):")
    response = response.split(',')
    try:
        indices = [int(x.strip())for x in response]
    except ValueError:
        print("Invalid input use integers only")
        return []
    operate = []
    for x in indices:
        if 1 <= x <= len(cmds):
            operate.append(cmds[x-1])
        else:
            print(f"Invalid number: {x}")
    print(f" the follwing {operate} will run")
    return operate

def main():
    # Main function to handle the preprocessing pipeline.
    args = parse_args()
    input_dir = args.input_dir
    if not os.path.exists(input_dir):
        os.makedirs(input_dir, exist_ok=True)
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    gtf_file = args.GTF
    genome_file = args.genome
    input_bam = os.path.join(output_dir, f"{sra_id}_aligned.bam")
    sra_ids = []
    if args.sra:
        sra_ids.append(args.sra)
    elif args.list:
        with open(args.list, 'r') as f:
            sra_ids = [line.strip() for line in f if line.strip()]
    for sra_id in sra_ids:
        operate = select_operations()
        for operation in operate:
            if operation in software_preprocess:
                log(f"Running {operation} for {sra_id}")
                print(f"Running {operation} for {sra_id}")
                software_preprocess[operation](sra_id, input_dir, output_dir)
                print(f"Finished {operation} for {sra_id}")
            else:
                print(f"Operation {operation} is not defined in the software_preprocess dictionary.")
                log(f"Operation {operation} is not defined in the software_preprocess dictionary.")


if __name__ == "__main__":
    main()
    print("Preprocessing completed.")
    log("Preprocessing completed.")
        
