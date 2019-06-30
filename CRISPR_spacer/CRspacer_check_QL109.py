#!/usr/local/bin/python3.6

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import os

# For CRISPR spacer design: generate seqs for BLAST checking the specificity.
# Add blast feature later

my_seq = Seq('GGTGAGCAGTCCACGAAGGGagg', IUPAC.unambiguous_dna)
nucldb_path = "/Users/durand.dc/Documents/works/Resources/Resource_QL109/Genome_QL109/blastdb_QL109/nt_QL109"
output_file = "/Users/durand.dc/Documents/works/Misc.files/temp/seq_sp_blast.out"
temp_file = "/Users/durand.dc/Documents/works/Misc.files/temp/seq_sp.fa"
blastn_program_path = "/Users/durand.dc/clitools/ncbi-blast-2.6.0+/bin/blastn"


def generate_4seq(seq):
	A_seq = SeqRecord(seq + Seq('AGG', IUPAC.unambiguous_dna), 'A', description = '')
	T_seq = SeqRecord(seq + Seq('TGG', IUPAC.unambiguous_dna), 'T', description = '')
	C_seq = SeqRecord(seq + Seq('CGG', IUPAC.unambiguous_dna), 'C', description = '')
	G_seq = SeqRecord(seq + Seq('GGG', IUPAC.unambiguous_dna), 'G', description = '')
	all_sp = []
	all_sp.append(A_seq)
	all_sp.append(T_seq)
	all_sp.append(C_seq)
	all_sp.append(G_seq)
	return all_sp

if len(my_seq) == 20:
	print(f'\nSpacer sequence:\n{my_seq}\nOrigional PAM unknown')
	rc_sp = my_seq.reverse_complement()
	my_seq_core = my_seq[8:]
	print(f'Core sequence: {my_seq_core}+PAM\n')
	print(f'Reverse complemented spacer:\n{rc_sp}\n')
	all_sp = generate_4seq(my_seq)
	spacer = my_seq
elif len(my_seq) == 23:
	print(f'\nSpacer sequence:\n{my_seq[:-3]}\nOrigional PAM = {my_seq[-3:]}')
	rc_sp = my_seq[:-3].reverse_complement()
	my_seq_core = my_seq[8:-3]
	print(f'Core sequence: {my_seq_core}+PAM\n')
	print(f'Reverse complemented spacer:\n{rc_sp}\n')
	all_sp = generate_4seq(my_seq[:-3])
	spacer = my_seq[:-3]
else:
	print('Wrong input!')
	raise

print('\nFor pCRISPomyces-2 one spacer design:')
print(f'Forward spacer: acgc{spacer}')
print(f'Reverse spacer: aaac{rc_sp}\n')

print('\nFor pCRISPR-Cas9 and pCRISPR-dCas9 one spacer design:')
print(f'CATGCCATGG{spacer}GTTTTAGAGCTAGAAATAGC')


SeqIO.write(all_sp, temp_file, 'fasta')

from Bio.Blast.Applications import NcbiblastnCommandline
num_threads = os.cpu_count()
run_blastn = NcbiblastnCommandline(blastn_program_path,
								   query = temp_file, 
								   out = output_file,
								   db = nucldb_path, 
								   remote = False,
								   evalue = 0.1,
								   num_threads = num_threads,
								   word_size = 7
								   )
stdout, stderr = run_blastn()
if stderr == '':
	print(f'\nSuccess!!\n{output_file}\n')
from subprocess import run

if os.name == 'posix':
	texteditor = ['open', '-a', '/Applications/TextEdit.app']
else:
	texteditor = ['notepad']
try:
	run([*texteditor, output_file])
except:
	pass