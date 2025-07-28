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
    return parser.parse_args()


def execute(cmd, outputfile = None, quitOnError = False):
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
    execute(f"fastq-dump {sra_id}", f"{sra_id}.fastq")
    execute(f"fastq-dump {sra_id} --split-files", f"{sra_id}_1.fastq")
    print(f"mv {sra_id}_1.fastq ")

def run_fastqc(sra_id, input_dir, output_dir):
    # Run FastQC on the downloaded FASTQ files to check quality of sequences.
    fastqc_out_dir = os.path.join(output_dir, "fastqc-out")
    fastqc_input = [
        f"{sra_id}_1.fastq",
        f"{sra_id}_2.fastq"
    ]
    fastqc_exec = f"fastqc -o {fastqc_out_dir} -t4"
    fastqc_cmd = f"{fastqc_exec} {' '.join(fastqc_input)}"
    execute(fastqc_cmd, f"{fastqc_out_dir}/{sra_id}_1_fastqc.html")

#list of software preprocess functions stored in a dictionary
software_preprocess = {
    "sra" : run_sra,
    "fastqc": run_fastqc
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
        


   

    
    
    
    
    