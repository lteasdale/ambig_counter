#!/usr/bin/env python

# Ambiguous sites counter!

from __future__ import print_function
from os import path
# The fasta parser.
import screed
# Does what a dictionary does, essentially the same as
# dicti[word] = dicti.get(word,0) + 1
from collections import Counter, OrderedDict
import docopt
import sys


__author__ = "lteasnail"

CLI_ARGS = """
USAGE:
ambig_counter.py [options] <FASTAFILES> ...

OPTIONS:
 -t THRESHOLD   If the proportion of ambiguous sites for a sequence is above
                this threshold the sequence will be replaced with a dummy
                sequence (all ~). The threshold needs to be given as a
                proportion e.g. 0.03 for a cutoff of 3%. Ifyou only want to
                count the ambiguous sites don't specify this flag.


Ambiguous sites counter

This program goes through a fasta file and outputs a tab dimlimitated txt
file with the file name, the sequence header,the length of the sequence -
including gaps, and the proportion of ambigous sites. If you specify a
threshold value (i.e. a proportion) you can also replace all sequences where
the proportion of ambiguous sites is greater than or equal to the threshold
with a dummy sequence.

You need to have python-screed and python-docopt installed.
"""

# This function counts the number of ambiguous bases per sequence per file,
# then works out the proportion of ambiguous sites and put them all in a
# dictionary


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

# This function replaces a sequence with a dummy if the proportion of
# ambiguous sites is over a threshold


def prop_ambig_threshold(files, threshold):
    for filename in files:
        ambigs = count_ambig(filename)
        file_basename = path.splitext(filename)[0]
        new_filename = file_basename + '_processed.fasta'
        new_file = open(new_filename, 'w')
        for seq in screed.open(filename):
            seq_name = seq.name
            seq = seq.sequence
            length, prop_ambig = ambigs[seq_name]
            if prop_ambig >= threshold:
                seq = '~' * length
                seq_name = seq_name + '_dummy_high_prop_ambig'
            print(">{}\n{}\n".format(seq_name, seq), file=new_file)
        new_file.close()


# This function runs the count_ambig function for each file and outputs the
# stats to a tab demlim file


def print_prop_ambig(files, table_file=sys.stdout):
    print('Filename\tSeq_header\tLength_of_seq\tProportion_ambiguous_sites',
          file=table_file)
    for filename in files:
        all_prop_ambig = count_ambig(filename)
        for seq_name, (length, prop_ambig) in all_prop_ambig.items():
            print(filename, seq_name, length, prop_ambig, sep='\t',
                  file=table_file)


# If I am being run as a script...
if __name__ == '__main__':
    opts = docopt.docopt(CLI_ARGS)
    threshold = opts['-t']
    files = opts['<FASTAFILES>']
    print('counting...', file=sys.stderr)
    print_prop_ambig(files)
    print('Finished counting ambiguous sites', file=sys.stderr)
    if threshold is not None:
        threshold = float(threshold)
        pc = threshold * 100
        print("Replacing sequences with >= {}% ambiguous sites...".format(pc),
              file=sys.stderr)
        prop_ambig_threshold(files, threshold)
        print("Finished replacing sequences", file=sys.stderr)
