import os
from Bio.Seq import Seq

blastn_program_path = r"D:\Program Files\NCBI\blast-2.6.0+\bin\blastn.exe"
blastp_program_path = r"D:\Program Files\NCBI\blast-2.6.0+\bin\blastp.exe"
nucldb_path = r"D:\WORKs\Resources\Streptomyces_genomes\blastdb\all_strep_201705\blastdb_nucl\stre_all_nucldb"
protdb_path = r"D:\WORKs\Resources\Streptomyces_genomes\MGB_DB\All_strep_201705"
output_path = r"D:\WORKs\Resources\temp\blast"
num_threads = os.cpu_count()

if not os.path.exists(output_path):
	os.makedirs(output_path)

# determine seq type from user input
def seq_type():
	seq_type_in = input('Protein(p) or nucleic(n)?')
	print(seq_type_in)
	if seq_type_in != 'p' and seq_type_in != 'n':
		seq_type = ''
		print("Don't understand, I will determine by myself.\n")
	elif seq_type_in == 'p':
		seq_type = 'prot'
	elif seq_type_in == 'n':
		seq_type = 'nucl'
	return(seq_type)

def input_seqs(seq_type):
	num_seqs = 0
	seq_titles = []
	temp_input_seq_file = os.path.join(output_path, 'temp_input_seq.fa')
	
	# This block is to find whether to perform search with existing file
	import pathlib
	last_blast = 'n'
	if pathlib.Path(temp_input_seq_file).is_file() and os.stat(temp_input_seq_file).st_size != 0:
		last_blast = input('Found last search, do you want to do that again?(y/n)')
	if last_blast != 'y':
		print('Please enter/paste your gene(s)/protein(s) in <fasta> format.\nThen start a new line, press Ctrl-Z(or D) + <Enter> to save it.\nInput:')
		input_sequences = []
		while True:
			try:
				line = input("")
			except EOFError:
				break
			input_sequences.append(line)
	else:
		with open(temp_input_seq_file, 'r') as f:
			input_sequences = f.readlines()

	with open(temp_input_seq_file, 'w') as input_file_handle:
		for line in input_sequences:
			if line.startswith('>'):
				input_file_handle.write(line.strip() + '\n')
				num_seqs += 1
				seq_titles.append(line[1:])
			else:
				input_file_handle.write(line.strip() + '\n')
				while num_seqs == 1 and seq_type == '':
					from Bio.Data.CodonTable import TranslationError
					try:
						for i in range(len(line)-2):
							test_seq = line[i:]
							test_seq += 'n' * (3 - len(test_seq) % 3)
							Seq(test_seq).translate()
						seq_type = 'nucl'
					except TranslationError:
						seq_type = 'prot'
	return temp_input_seq_file, num_seqs, seq_titles, seq_type


def blastn(temp_input_seq_file, output_path, outfile):			
	from Bio.Blast.Applications import NcbiblastnCommandline
	out = os.path.join(output_path, outfile)
	run_blastn = NcbiblastnCommandline(blastn_program_path,
										query = temp_input_seq_file, 
										out = out,
										# outfmt = 5, 
										db = nucldb_path, 
										remote = False,
										# evalue = 1e-30,
										num_threads = num_threads
										)
	print(f'Running BLASTN command:\n{run_blastn}')
	stdout, stderr = run_blastn()
	if stderr == '':
		print(f'\nSuccess!!\n{out}\n')
	from subprocess import run
	try:
		run(['notepad', out])
	except:
		pass

def blastp(temp_input_seq_file, output_path, outfile):			
	from Bio.Blast.Applications import NcbiblastpCommandline
	out = os.path.join(output_path, outfile)
	run_blastp = NcbiblastpCommandline(blastp_program_path,
										query = temp_input_seq_file, 
										out = out,
										# outfmt = 5, 
										db = protdb_path, 
										remote = False,
										# evalue = 1e-30,
										num_threads = num_threads
										)
	print(f'Running BLASTP command:\n{run_blastp}')
	stdout, stderr = run_blastp()
	if stderr == '':
		print(f'\nSuccess!!\n{out}\n')
	from subprocess import run
	try:
		run(['notepad', out])
	except:
		pass

		
temp_input_seq_file, num_seqs, seq_titles, seq_type = input_seqs(seq_type())
import re
file_name = re.sub('[^\w_.)( -]', '.', seq_titles[0])[:30]
if seq_type == 'prot':
	blastp(temp_input_seq_file, output_path, outfile = f'allStrep201705_blastp_{file_name}_{num_seqs if num_seqs != 1 else ""}{"seqs_"if num_seqs != 1 else ""}.out')
elif seq_type == 'nucl':
	blastn(temp_input_seq_file, output_path, outfile = f'allStrep201705_blastn_{file_name}_{num_seqs if num_seqs != 1 else ""}{"seqs_"if num_seqs != 1 else ""}.out')
else:
	print('No BLAST done')