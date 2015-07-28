#!/usr/bin/env python

from __future__ import print_function

# The fasta parser.
import screed
# For taking arguments from the command line.
import sys
# Does what a dictionary does, essentially the same as
# dicti[word] = dicti.get(word,0) + 1
from collections import Counter, OrderedDict


"""
Ambiguous sites counter

need to have python-screed installed
outputs the file name, the sequence name,the length of the sequence -
including gaps and the proportion of ambigous sites

Run as follows:

    $ python ambig_counter.py file.fasta header_seq1 header_seq2 > output.txt
"""

__author__ = "lteasnail"


def count_ambig(filename):
    ambigs = OrderedDict()
    for seq in screed.open(filename):
        seq_name = seq.name
        seq = seq.sequence
        length = len(seq)

        counter = Counter()
        for base in seq:
            counter[base] += 1

        num_unambiguous_bp = 0
        unambiguous_bases = "ACGTacgtNnxX-~"
        for base_type in unambiguous_bases:
            num_unambiguous_bp += counter.get(base_type, 0)
        prop_ambig = (length - num_unambiguous_bp) / float(length)

        ambigs[seq_name] = (length, prop_ambig)
    return ambigs


def prop_ambig_threshold(files, threshold):
    for filename in files:
        ambigs = count_ambig(filename)
        for seq in screed.open(filename):
            seq_name = seq.name
            seq = seq.sequence
            length, prop_ambig = ambigs[seq_name]
            if prop_ambig >= threshold:
                seq = '~' * length
            print(">seq_{}\n{}\n".format(seq_name, seq))


def print_prop_ambig(files, output_file=sys.stdout):
    for filename in files:
        all_prop_ambig = count_ambig(filename)
        for seq_name, (length, prop_ambig) in all_prop_ambig.items():
            print(filename, seq_name, length, prop_ambig, sep='\t',
                  file=output_file)


# If I am being run as a script...
if __name__ == '__main__':
    threshold = sys.argv[1]
    files = sys.argv[2:]
    print_prop_ambig(files)
