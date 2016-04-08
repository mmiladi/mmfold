import numpy as np

import re
import glob
import os, sys
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
import itertools

import pandas as pd

from subprocess import *

# settings !!!! IMPORTANT UPDATE IT to CisDataset folder !!!!!!
CIS_GENOME_PATH = '/home/milad/DataBase/CisReg/Cis_include_genome2/'
CIS_MRNA_PATH = '/home/milad/DataBase/CisReg/Cis_include_mRNA/'


VIENNA_BIN_PATH = '/home/milad/software/bin/'
QUAKE_PARAM_FILE = '/home/milad/workspace/mmfold/src/misc/rna_turner2004_ML_up_penalty.par '
RNAFOLD = 'RNAfold -p --noPS '
RNAPLFOLD = 'RNAplfold '


# -------------
def remove_gap_columns(malign):
    ''' Remove all-gap columns from a biopython multiple-alignment
    Returns pruned multiple-alignment'''
    regx_gaps = '[-.~_]'  # valid gap symbols all be converted to "-"

    if (len(malign) == 0):
        return
    for c in reversed(range(len(malign[0]))):
        if len(re.sub(regx_gaps, '', malign[:, c])) == 0:
            malign = malign[:, 0:c] + malign[:, c + 1:]  # concat left & right side of column c

    return malign

# End def remove_gap_columns


def sub_dotplot(dp, lcontext, rcontext):
    '''Returns an squared submatrix by removing pre and post section ....
    the aim is to extract dotplot of RNA from a dotplot of the context extened RNA'''
    assert (dp.shape[0] == dp.shape[1])
    assert (dp.shape[0] >= lcontext + rcontext)
    n = dp.shape[0]
    return dp[lcontext:n-rcontext, lcontext:n-rcontext]


def upper_part(dpX):
    ''' Returns the upper section of matrix in flattened form, expecting to have bp probabilities'''

    return dpX[np.triu_indices(dpX.shape[0], 1)]

# TODO: Add all repeatedly used functions into a library accros all tools


def dotbracket_to_dict(struct):
    '''Returns a dictionary where basepairs are keys with !ONE! based indices joined by ":" ,
    e.g. dict {'0:10': 1, '2:8': 1} '''
    assert len(struct.replace('.', '').replace('(', '').replace(')', '')) == 0
    stack = list()
    pairs = dict()
    for pos, ch in enumerate(list(struct)):
        #         print pos+1, ch
        if ch == '(':
            stack.append(pos)
        elif ch == ')':
            left = stack.pop()
            key = "{}:{}".format(left+1, pos+1)
            pairs[key] = 1

    assert len(stack) == 0
    return pairs


def compute_part_func(infile_fa, seq_names, outdir_path="./", use_plfold=False, quake_params=False,
                      use_cache=False):
    '''Runs Vienna RNAfold/RNAplfold with partition function for all sequences inside input fasta file
    If use_cache, it does nothing if If the ps file with same paramaters exists '''
    from subprocess import Popen, PIPE
    #     print "compute_part_func(", infile_fa, seq_names
    if use_plfold:
        out_dir = outdir_path + RNAPLFOLD.replace(' ', '')
    else:
        out_dir = outdir_path + RNAFOLD.replace(' ', '')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if not os.path.isfile(infile_fa):
        raise IOError("Fastafile not found: {}".format(infile_fa))

    all_in_cache = all([os.path.isfile(os.path.join(out_dir, sname+'_dp.ps')) for sname in seq_names])
    if all_in_cache and use_cache:
        raise NotImplementedError("Sequence names for caching are not correctly set")
        return out_dir

    with open(infile_fa) as in_rna:

        arg_param = ""
        if quake_params:
            arg_param = " -P %s " % QUAKE_PARAM_FILE

        if use_plfold:
            p = Popen(('cd %s;' %out_dir) + VIENNA_BIN_PATH + RNAPLFOLD + arg_param, stdin=in_rna, shell=True, stdout=PIPE, stderr=PIPE)
        else:
            p = Popen(('cd %s;' %out_dir) + VIENNA_BIN_PATH +  RNAFOLD + arg_param, stdin=in_rna, shell=True, stdout=PIPE, stderr=PIPE)

        out, err = p.communicate()
        if err:
            print "Error in calling RNAfold for ", infile_fa
            print out
            print err

            # With long sequences RNAfold prints scalign factor to stderr
            if (not use_plfold and not ("scaling factor" in err or "free energy" in err)):
                raise RuntimeError

    return out_dir


def parse_dp_ps(ps_file):
    '''Extracts base pair probabliies from vienna ps file
    returns: Dictinary of form dict[i:j]=p(i,j) '''

    # Extract sequence from ps file
    myseq = ""
    read_seq = False
    with open(ps_file) as in_ps:
        for line in in_ps:
            if "/sequence" in line:
                read_seq = True
            elif read_seq and ") } def" in line:
                read_seq = False
            elif read_seq:
                myseq += line.rstrip().rstrip("\\")
    #     print ps_file.rstrip("_dp.ps") , myseq

    ureg = re.compile(r'^(\d+)\s+(\d+)\s+(\d+\.\d+)\s+ubox\s*')
    bp_prob_dict = dict()
    bp_prob_mat = np.zeros((len(myseq), len(myseq)))

    with open(ps_file) as in_ps:
        for line in in_ps:
            if "ubox" in line:
                um = ureg.match(line)
                if um:
                    i, j, sqrp = um.groups()

                    #                     print i, j, sqrp

                    # keys are pair of indexes as smaller:larger
                    key = ":".join(sorted([i, j], reverse=True))
                    assert (key not in bp_prob_dict)
                    bpprob = float(sqrp)*float(sqrp)
                    bp_prob_dict[key] = bpprob

                    i, j = int(i), int(j)
                    bp_prob_mat[i-1, j-1] = bpprob
    return bp_prob_mat


def get_expected_accuracy(reference_struct, dp_matrix):
    '''dp_matrix is a numpy matrix where base indeices are ZERO based'''
    assert dp_matrix.shape[0] == dp_matrix.shape[1]
    assert dp_matrix.shape[0] == len(reference_struct)
    reference_struct_dict = dotbracket_to_dict(reference_struct)
    sum_TP_prob = 0.0
    for bp_key in reference_struct_dict:
        i, j = bp_key.split(":")
        i, j = int(i), int(j)
        sum_TP_prob += dp_matrix[i-1, j-1]
#         print i,j, dp_matrix[i-1,j-1]

#     print "    TP_score: %.2f" % (sum_TP_prob/len(reference_struct_dict))
    if len(reference_struct_dict) == 0:
        return 0
    return (sum_TP_prob/len(reference_struct_dict))


def get_left_right_context_lengths(motif_id, context_id):

    motif_splits = motif_id.replace('/', ' ').replace('_', ' ').replace('-', ' ').split()
    assert len(motif_splits) == 3

    context_splits = context_id.replace('/', ' ').replace('_', ' ').replace('-', ' ').split()
    assert len(context_splits) == 3

    motif_acc, motif_start, motif_end = motif_splits
    context_acc, context_start, context_end = context_splits

    if motif_acc != context_acc:
        raise RuntimeError("Mismtach motif and context accesions {} {}".format(motif_id, context_id))

    motif_start, motif_end = motif_splits[1:]
    context_start, context_end = context_splits[1:]
    motif_start, motif_end, context_start, context_end = [(int)(s) for s in
                                                          [motif_start, motif_end, context_start, context_end ]]
#     print motif_id, context_id, motif_start, motif_end, context_start, context_end
    context_len_left = motif_start - context_start
    context_len_right = context_end - motif_end
#     print  motif_id, context_id, "context_len_left {}, context_len_right {}".format(context_len_left, context_len_right)

    # Verify and adpat reverse strands
    # Sorry for the complication implemented below and imposed by the accesion encoding 
    on_reverse = False
    if motif_start > motif_end:
        if context_start > context_end:
            #       print "reverse strand"
            on_reverse = True
        else:
            raise RuntimeError("Mismatch1 of context right left positions for seq {}".format(motif_id))
    else:
        if context_start > context_end:
            raise RuntimeError("Mismatch2 of context right left positions for seq {}".format(motif_id))


#         print on_reverse
    if on_reverse:
        context_len_left *= -1
        context_len_right *= -1
    if not on_reverse:
        assert motif_start < motif_end
        assert context_start < context_end
        assert context_start <= motif_start
        assert motif_end <= context_end
    else:
        assert motif_start > motif_end
        assert context_start > context_end
        assert context_start >= motif_start
        assert motif_end >= context_end

    motif_len = abs(motif_end - motif_start) + 1
    return context_len_left, context_len_right, motif_len, context_start, context_end, motif_acc


def generate_assymetric_fasta(motif_seq_id, target_context_len, fasta_supercontext, output_path):

    if not os.path.isfile(fasta_supercontext):
        raise IOError("Not found fasta_supercontext:{} ".format(fasta_supercontext))
    # Open supercontext to get the id as well as the super-sequence
    with open(fasta_supercontext, "r") as in_fasta_handle:
        fa_recs = list(SeqIO.parse(in_fasta_handle, "fasta"))
    assert len(fa_recs) == 1
    super_seq = fa_recs[0]
    super_seq_str = str(super_seq.seq)

    super_left_len, super_right_len, motif_seq_len, super_start, super_end, accession = get_left_right_context_lengths(motif_seq_id, super_seq.id)

    if super_left_len + super_right_len < target_context_len:  # TODO: Check this case where context is missing
        #         raise RuntimeError
        print ("ERROR: Not enough super context available l:{} r:{} target:{}".format(super_left_len, super_right_len, target_context_len))
        len_to_remove = left_len_to_remove = right_len_to_remove = 0
        out_seq_str = super_seq_str  # Get whole availble
        # TODO: What to do with these cases? For compatibility with CisDatasetreturn all available, otherswise return
    else:
        len_to_remove = (super_left_len + super_right_len) - target_context_len

        avg_left_ratio = super_left_len / (float)(super_left_len + super_right_len)

        left_ratio = min(1.0, max(0, np.random.normal(avg_left_ratio, 0.1)))
    #     left_ratio =  avg_left_ratio
        left_len_to_remove = min(super_left_len, int(left_ratio * len_to_remove))
        right_len_to_remove = len_to_remove - left_len_to_remove
        if right_len_to_remove > super_right_len:
            left_len_to_remove = len_to_remove - super_right_len
            right_len_to_remove = super_right_len

        assert 0 <= right_len_to_remove and right_len_to_remove <= super_right_len
        assert 0 <= left_len_to_remove and left_len_to_remove <= super_left_len
#         print "right_len_to_remove + left_len_to_remove == len_to_remove", right_len_to_remove , left_len_to_remove , len_to_remove
        assert right_len_to_remove + left_len_to_remove == len_to_remove
        assert len(super_seq_str) == (len_to_remove + target_context_len + motif_seq_len)
#         print "Range: {}:{}".format( left_len_to_remove, (left_len_to_remove+target_context_len+motif_seq_len))
        out_seq_str = super_seq_str[left_len_to_remove:(left_len_to_remove+target_context_len+motif_seq_len)]
        assert len(out_seq_str) == target_context_len + motif_seq_len

    out_fasta_file = "{}/{}".format(output_path, os.path.basename(fasta_supercontext))
    if super_start < super_end:  # Direct strand
        out_seq_id = "{}_{}-{}".format(accession, super_start+left_len_to_remove, super_end-right_len_to_remove)
    else:  # reverse strand
        out_seq_id = "{}_{}-{}".format(accession, super_start-left_len_to_remove, super_end+right_len_to_remove )

    with open(out_fasta_file, "w") as out_fasta_handle:
        #         SeqIO.write([], out_fasta_file, "fasta")
        out_fasta_handle.write(">{}\n{}\n".format(out_seq_id, out_seq_str))


def get_extended_fasta_file(seq, famid, flank_len, dataset='genome', use_assymetric_context=False):
    # Get the fasta file of specific flanking range
    # TODO: Waht to do with unknown_nt i.e. sequnces with unknown nucleutides

    if dataset == 'genome':
        db_path = CIS_GENOME_PATH
    elif dataset == 'mrna':
        db_path = CIS_MRNA_PATH
    else:
        raise RuntimeError("Unknown dataset type: {}".format(dataset))
    if use_assymetric_context:
        fam_flank_path = '{}{}/Cis_flanks-Assym{}/'.format(db_path, famid, flank_len)
    else:
        fam_flank_path = '{}{}/Cis_flanks-{}/'.format(db_path, famid, flank_len)

    if not os.path.isdir(fam_flank_path):
        raise IOError("Family flanking dir does not exist: {}".format(fam_flank_path))

    fasta_flanked_seq = "{}/{}_known_nt.fasta".format(fam_flank_path, seq.id.replace("/", "_"))
    # Sometimes the starting ending positions are swapped, redfine fasta var
    if not os.path.isfile(fasta_flanked_seq):
        fasta_flanked_seq_unknown = fasta_flanked_seq.replace('_known_', '_unknown_')

        if os.path.isfile(fasta_flanked_seq_unknown):
            print "Warning: benchmarking flanked sequnce with unknown nucleotides {}".format(fasta_flanked_seq_unknown)
            fasta_flanked_seq = fasta_flanked_seq_unknown

        else:  # reverse case
            splits = seq.id.replace("/", " ").replace("-", " ").split()
            assert(len(splits) == 3)

            fasta_flanked_seq = "{}/{}_known_nt.fasta".format(fam_flank_path, "{}_{}-{}".format(splits[0], splits[2], splits[1]) )
            if not os.path.isfile(fasta_flanked_seq):
                fasta_flanked_seq_unknown = fasta_flanked_seq.replace('_known_', '_unknown_')
                if os.path.isfile(fasta_flanked_seq_unknown):
                    print "Warning: benchmarking flanked sequnce with unknown nucleotides {}".format(fasta_flanked_seq_unknown)
                    fasta_flanked_seq = fasta_flanked_seq_unknown
                else:  # none of : direct, reverse, known, unknown
                    raise IOError("Fasta file not found: {}".format(fasta_flanked_seq))
    #             print "reverse strand"

    return fasta_flanked_seq

