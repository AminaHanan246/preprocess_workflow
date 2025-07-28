import os
import sys
import subprocess
import argparse
from datetime import datetime

def log(message,logfile="preprocess_bulk.log"):
    timestamp = datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    with open(logfile, "a") as f:
        f.write(f"{timestamp}: {message}\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Bulk RNA-seq preprocessing pipeline.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sra", type=str, help="SRA ID to process (e.g. SRR1234567)")
    group.add_argument("--list", type=str, help="Path to file with one SRA ID per line")
    parser.add_argument("--input_dir", type=str, default="./data", help="Directory to store raw files")
    parser.add_argument("--output_dir", type=str, default="./output", help="Directory to store FastQC output")
    return parser.parse_args()


def execute(cmd, outputfile = None, quitOnError = False):
    """The preprocessing function 'x' ."""
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

if not os.path.exists(input_dir):
    os.makedirs(input_dir, exist_ok = True)
if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok = True)

def run_sra(sra_id, input_dir, output_dir):
    """Download SRA data and convert to FASTQ format."""
    if not os.path.exists(input_dir):
    os.makedirs(input_dir, exist_ok = True)
    if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok = True)
    execute(print(f"fastq-dump {sra_id}"),print(f"{sra_id}.fastq"))
    execute(print(f'fastq-dump {sra_id} --split-files'),print(f"{sra_id}_1.fastq"))
    print(f"mv {sra_id}_1.fastq ")

def run_fastqc(sra_id, input_dir, output_dir):
    fastqc_out_dir = os.path.join(output_dir, "fasqc-out")
    fastqc_input = [
        f"{sra_id}_1.fastq",
        f"{sra_id}_2.fastq"
    ]
    fastqc_exec = f"fastqc -o {fastqc_out_dir} -t4"
    fastqc_cmd = f"{fastqc_exec} {' '.join(fastqc_input)}"
    execute(fastqc_cmd, f"{fastqc_out_dir}/{input_base}_1_fastqc.html")


software_preprocess = {
    "sra" : run_sra,
    "fastqc": run_fastqc
}

def select_operations():
    cmds = ['sra','fastqc']
    print("Available operations: {cmds}")
    print("Please select the operations you want to run by entering their numbers separated by commas (e.g., 1,3,5):")
    
    operate = []
    for item in cmds:
        while True:
            response = input(f"Run {item}? (y/n): ").lower()
            if response == 'y':
                operate.append(item)
                break
            elif response == 'n':
                print(f"Skipping {item}.")
                break
            else:
                print("Invalid input, please enter 'y' or 'n'.")
                continue

    print(f" the follwing {operate} will run")

def main():
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    sra_ids = []
    if args.sra:
        sra_ids.append(args.sra)
    elif args.list:
        with open(args.list, 'r') as f:
            sra_ids = [line.strip() for line in f if line.strip()]
    for sra_id in sra_ids:
        select_operations()
        for operation in operate:
            if operation in software_preprocess:
                log(f"Running {operation} for {sra_id}")
                software_preprocess[operation](sra_id, input_dir, output_dir)
            else:
            print(f"Operation {operation} is not defined in the software_preprocess dictionary.")
            log(f"Operation {operation} is not defined in the software_preprocess dictionary.")
        