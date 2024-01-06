import argparse, os
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
import utils

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=str, default="RV11",help='Name of BALIBASE reference set')
    parser.add_argument("--alignmode",type=str, default="global", help="Alignment mode, global or local")
    parser.add_argument("--simscalefactor", type=float,default=10.0)
    return parser.parse_args()

def align_and_write_msf(seq1,seq2):
    ...

if __name__ == '__main__':
    args = parse_arguments()
    base_data_path = "data/sequences"
    balibase_ref_path = os.path.join(base_data_path,args.ref)
    msas = os.listdir(balibase_ref_path)

    # Create the directories
    blosum62_path = "data/alignments"+"/" + args.alignmode +"_blosum62" 
    if not os.path.exists(blosum62_path):
        os.makedirs(blosum62_path)
    blosum62_ref_path = blosum62_path + "/" + args.ref
    if not os.path.exists(blosum62_ref_path):
        os.makedirs(blosum62_ref_path)

    aaembedding_path = "data/alignments"+"/" + args.alignmode +"_aaembedding"
    if not os.path.exists(aaembedding_path):
        os.makedirs(aaembedding_path)
    aaembedding_ref_path = aaembedding_path + "/" + args.ref
    if not os.path.exists(aaembedding_ref_path):
        os.makedirs(aaembedding_ref_path)


    embedding_similarity_matrix = utils.get_scaled_similarity_matrix(factor=args.simscalefactor)
    blosum62_m = substitution_matrices.load("BLOSUM62")
    penalty_params = {"penalty_type": "affine", "gap_penalty": -11, "gap_extension_penalty": -1}

    for msa in msas:
        print("ON MSA:", msa)
        msa_path = os.path.join(balibase_ref_path,msa)
        sequence_list = os.listdir(msa_path)

        msa_alignment_directory_blosum = blosum62_ref_path + "/" + msa
        if not os.path.exists(msa_alignment_directory_blosum):
            os.makedirs(msa_alignment_directory_blosum)

        msa_alignment_directory_aaembedding = aaembedding_ref_path + "/" + msa
        if not os.path.exists(msa_alignment_directory_aaembedding):
            os.makedirs(msa_alignment_directory_aaembedding)

        # Loop to detect pairs of sequences
        for i in range(len(sequence_list)):
            for j in range(i+1,len(sequence_list)):
                id1 = sequence_list[i].split('.')[0] 
                id2 = sequence_list[j].split('.')[0]
                raw_seq1 = SeqIO.read(os.path.join(msa_path,sequence_list[i]), "fasta")
                raw_seq2 = SeqIO.read(os.path.join(msa_path,sequence_list[j]), "fasta")

                embedding_alignment = utils.align_sequences_biopython(raw_seq1, raw_seq2, penalty_params, embedding_similarity_matrix,mode=args.alignmode)
                blosum62_alignment = utils.align_sequences_biopython(raw_seq1, raw_seq2, penalty_params, blosum62_m,mode=args.alignmode)
                
                # EMBEDDING BASED ALIGNMENT
                seq1 = embedding_alignment[0][0].replace('-', '.')
                seq2 = embedding_alignment[0][1].replace('-', '.')
                if args.alignmode == "global":
                    # In case of global alignment, beginning indices are 0 and end is the len
                    beg, end = [0, 0], [len(seq1), len(seq2)]
                elif args.alignmode == "local":
                    # TODO: local alignments
                    ...
                utils.write_msf(seq1=seq1,seq2=seq2,id1=id1,id2=id2,method = "embeddingmatrix",gopen=-11.0, gext=-1.0,path= msa_alignment_directory_aaembedding, beg=beg, end=end )

                # BLOSUM ALIGNMENT
                seq1 = blosum62_alignment[0][0].replace('-', '.')
                seq2 = blosum62_alignment[0][1].replace('-', '.')
                if args.alignmode == "global":
                    # In case of global alignment, beginning indices are 0 and end is the len
                    beg, end = [0, 0], [len(seq1), len(seq2)]
                elif args.alignmode == "local":
                    # TODO: local alignments
                    ...
                utils.write_msf(seq1=seq1,seq2=seq2,id1=id1,id2=id2,method = "blosum62",gopen=-11.0, gext=-1.0,path= msa_alignment_directory_blosum, beg=beg, end=end )
    
