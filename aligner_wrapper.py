#!/usr/bin/env python3

import argparse, subprocess, os

def checkFASTQ(list1, list2):

    #Checking if lengths are the same, since they are paired end reads:
    if len(list1) == len(list2):
        #Since the lengths are the same, check if the files end with _1.fastq and _2.fastq
        for i in range(0, len(list1)):
            if list1[i].endswith("_1.fastq") and list2[i].endswith("_2.fastq"): #(Assuming file ends with .fastq and not .fq, add in .fq later)
                status = True
    else:
        status = False

    return status

def runKallisto(referenceindex, output_directory, input_directory, file1, file2):
    print("Running Kallisto for: "+file1.split("/")[-1].replace("_1.fastq", ""))
    subprocess.run("kallisto quant -i "+referenceindex+" -o "+output_directory+file1.split("/")[-1].replace("_1.fastq", "")+" "+input_directory+file1+" "+input_directory+file2+" -b 100", shell=True)

def runHISAT2(referenceindex, output_directory, input_directory, file1, file2):
    print("Running HISAT2 for: "+file1.split("/")[-1].replace("_1.fastq", ""))
    subprocess.run("hisat2 -x "+referenceindex+" -1 "+input_directory+file1+" -2 "+input_directory+file2+" -S "+output_directory+file1.split("/")[-1].replace("_1.fastq", "")+".sam", shell=True)

def main():

    #Argparse code
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--aligner-to-use", help="1: Uses Kallisto; 2: Uses HISAT2; Default: 1", required=True, type=int, default=1)
    parser.add_argument("-i", "--input-files-directory", help="Path to directory containing FASTQ files.", required=True)
    parser.add_argument("-r", "--reference-index", help="Path of indexed reference file.\n"
                                                         "NOTE: 1. If Kallisto is selected, then enter the path of indexed reference transcriptome (Ends with .idx);\n"
                                                         "2. If HISAT2 is selected, then enter the path of indexed reference genome with the prefix", required=True)
    parser.add_argument("-o", "--output-files-directory", help="Enter name of directory which will contain output of aligner.", required=True)
    args = vars(parser.parse_args())

    #Populating the variables
    aligner = args['aligner_to_use']
    input_path = args['input_files_directory']
    referenceindex_path = args['reference_index']
    output_path = args['output_files_directory']

    #Stripping "/" if present and adding "/", if not present will add "/" anyway: (Can do it in the subprocess commands itself, edit later)
    input_path = input_path.rstrip("/")
    input_path = input_path + "/"
    output_path = output_path.rstrip("/")
    output_path = output_path + "/"

    #Making output directory:
    if os.path.isdir(output_path) is False:
        os.mkdir(output_path)

    input_dir = sorted(os.listdir(input_path))
    mate1_list = [x for x in input_dir if "_1" in x]
    mate2_list = [x for x in input_dir if "_2" in x]

    #Checking FASTQ files:
    status = checkFASTQ(mate1_list, mate2_list)

    if status is True:
        if aligner == 1:
            #Running kallisto:
            for i,j in zip(mate1_list, mate2_list):
                runKallisto(referenceindex_path, output_path, input_path, i, j)
        if aligner == 2:
            #Running hisat2:
            for i,j in zip(mate1_list, mate2_list):
                runHISAT2(referenceindex_path, output_path, input_path, i, j)
    else:
        print("Error running aligners, check your files.")

if __name__ == "__main__":
    main()
