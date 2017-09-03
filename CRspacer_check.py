from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import os

# For CRISPR spacer design: generate seqs for BLAST checking the specificity.
# Add blast feature later

my_seq = Seq('CTCCGGGCAGATGAAAGGTCagg', IUPAC.unambiguous_dna)
nucldb_path = r"D:\WORKs\Resources\Streptomyces_genomes\blastdb\M145"
output_file = r"D:\WORKs\Resources\temp\seq_sp_blast.out"
temp_file = r"D:\WORKs\Resources\temp\seq_sp.fa"

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
elif len(my_seq) == 23:
	print(f'\nSpacer sequence:\n{my_seq[:-3]}\nOrigional PAM = {my_seq[-3:]}')
	rc_sp = my_seq[:-3].reverse_complement()
	my_seq_core = my_seq[8:-3]
	print(f'Core sequence: {my_seq_core}+PAM\n')
	print(f'Reverse complemented spacer:\n{rc_sp}\n')
	all_sp = generate_4seq(my_seq[:-3])
else:
	print('Wrong input!')
	raise

print('\nFor pCRISPomyces-2 one spacer design:')
print(f'Forward spacer: acgc{my_seq[:-3]}')
print(f'Reverse spacer: aaac{rc_sp}\n')

SeqIO.write(all_sp, temp_file, 'fasta')
blastn_program_path = r"D:\Program Files\NCBI\blast-2.6.0+\bin\blastn.exe"

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
try:
	run(['notepad', output_file])
except:
	pass