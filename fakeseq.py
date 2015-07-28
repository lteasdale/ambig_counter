from random import choice as randchoice
from ambig_counter import main as count_ambig


def make_magic_alignment(name):
    bp = [randchoice('ACGTN-') for x in range(100)]
    seq = ''.join(bp)
    return ">seq_{}\n{}\n".format(name, seq)


aln_files = []
for x in range(1000):
    aln = make_magic_alignment(x)
    aln_fname = "fakeseq/aln_{}.fasta".format(x)
    with open(aln_fname, 'w') as fh:
        fh.write(aln)
    aln_files.append(aln_fname)

stat_fname = "fakeseq/all.tab".format(x)
with open(stat_fname, 'w') as tabfh:
    count_ambig(aln_files, output_file=tabfh)
