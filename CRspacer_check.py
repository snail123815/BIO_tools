# For CRISPR spacer design: generate seqs for BLAST checking the specificity.
# Add blast feature later

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC

my_seq = Seq('ACAGGTCCTCCAGCATGAAGCGG', IUPAC.unambiguous_dna)
print(f'\nSpacer sequence:\n{my_seq[:-3]}\nOrigional PAM = {my_seq[-3:]}')
rc_sp = my_seq[:-3].reverse_complement()
my_seq = my_seq[8:-3]
print(f'Core sequence: {my_seq}+PAM\n')
print(f'Reverse complemented spacer:\n{rc_sp}\n')

A_seq = SeqRecord(my_seq + Seq('AGG', IUPAC.unambiguous_dna), 'A', description = '')
T_seq = SeqRecord(my_seq + Seq('TGG', IUPAC.unambiguous_dna), 'T', description = '')
C_seq = SeqRecord(my_seq + Seq('CGG', IUPAC.unambiguous_dna), 'C', description = '')
G_seq = SeqRecord(my_seq + Seq('GGG', IUPAC.unambiguous_dna), 'G', description = '')

all_sp = []
all_sp.append(A_seq)
all_sp.append(T_seq)
all_sp.append(C_seq)
all_sp.append(G_seq)

SeqIO.write(all_sp, 'seq_sp.fa', 'fasta')

#blastn -query seq_sp.fa -db '/media/e/WORKs/Resource_QL109/blastdb_QL109/nt_QL109' -out temp/sp_blast.out -word_size 7 -evalue 0.1

