import argparse, os
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
import subprocess
import time

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=str, default="RV11",help='Name of BALIBASE reference set')
    parser.add_argument("--alignmode",type=str, default="local", help="Alignment mode, global or local")
    parser.add_argument("--repdetails", type=str,default="proteinbert_full")

    return parser.parse_args()


if __name__ == '__main__':
    start_time = time.time()

    args = parse_arguments()
    base_data_path = "data/sequences"
    balibase_ref_path = os.path.join(base_data_path,args.ref)
    msas = os.listdir(balibase_ref_path)
    representations_base_folder = "C:/Users/efeka/Documents/GitHub/549-term-project/protein_embeddings/output_representations/" 


    rep_peba_alignments_path = "data/alignments"+"/" + args.alignmode +'_' + args.repdetails
    if not os.path.exists(rep_peba_alignments_path):
        os.makedirs(rep_peba_alignments_path)
    rep_peba_alignment_ref_path = rep_peba_alignments_path + "/" + args.ref
    if not os.path.exists(rep_peba_alignment_ref_path):
        os.makedirs(rep_peba_alignment_ref_path)

    index = 0
    for msa in msas:
        print("ON MSA:", msa, index)
        index += 1
        if index == 4:
            break
        msa_path = os.path.join(balibase_ref_path,msa)
        sequence_list = os.listdir(msa_path)

        msa_alignments_directory = rep_peba_alignment_ref_path + "/" + msa # path to write alignments
        if not os.path.exists(msa_alignments_directory):
            os.makedirs(msa_alignments_directory)

        # Loop to detect pairs of sequences
        for i in range(len(sequence_list)):
            for j in range(i+1,len(sequence_list)):
                id1 = sequence_list[i].split('.')[0] 
                id2 = sequence_list[j].split('.')[0]
                seq1_path = os.path.join(msa_path,sequence_list[i])
                seq2_path = os.path.join(msa_path,sequence_list[j])
                
                rep1_path = representations_base_folder + args.ref + '/' + msa + '/' + id1 + ".txt"
                rep2_path = representations_base_folder + args.ref + '/' + msa + '/' + id2 + ".txt"

                command = ["python", "peba/peba.py", "-a", args.alignmode,"-f1", seq1_path,"-f2", seq2_path, "-e1", rep1_path, "-e2", rep2_path, "-s", msa_alignments_directory]
                # Call peba for each pair
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                # Wait for the process to finish and get the return code
                return_code = process.wait()

                # Check if the subprocess ran successfully
                if return_code == 0:
                    print("Script executed successfully")
                else:
                    print(f"Script failed with return code {return_code}")

                # Optionally, you can capture and print the output and errors
                output, errors = process.communicate()
                print("Output:", output.decode("utf-8"))
                print("Errors:", errors.decode("utf-8"))
        
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")