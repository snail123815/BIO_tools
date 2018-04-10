from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline, NcbiblastxCommandline, NcbitblastnCommandline, NcbitblastxCommandline
from os import path, stat, environ, getenv, makedirs
from time import time
from subprocess import Popen
from sys import exit
from platform import system
from csv import reader as csvreader
from os import cpu_count
import pathlib
import re

cwd = path.dirname(path.realpath(__file__))
output_path = path.join(cwd, "BLAST_Results")
temp_input_seq_file = path.join(output_path, 'temp_input_seq.fa')
input_seq_file = path.join(cwd, "BLAST_INPUT_fasta.txt")
print(f"\nPut your input <fasta format> sequence (multiple sequences allowed) in <{input_seq_file}>\nOR paste your sequence later.")

if system() in ['Darwin', 'Linux']:
	print('\nThis wrapper needs blast programs setted up in your PATH.')
	blastn_program_path  = 'blastn'
	blastp_program_path  = 'blastp'
	blastx_program_path  = 'blastx'
	tblastn_program_path = 'tblastn'
	tblastx_program_path = 'tblastx'
	if system() == 'Linux':
		fileReader = ['gedit']
	else:
		fileReader = ['open', '-a', 'TextWrangler']
else: # 'Windows'
	blastn_program_path  = path.join(cwd, "blastProgram/blastn.exe")
	blastp_program_path  = path.join(cwd, "blastProgram/blastp.exe")
	blastx_program_path  = path.join(cwd, "blastProgram/blastx.exe")
	tblastn_program_path = path.join(cwd, "blastProgram/tblastn.exe")
	tblastx_program_path = path.join(cwd, "blastProgram/tblastx.exe")
	fileReader = ['notepad']

def seq_type():
	seq_type_in = ''
	while seq_type_in not in ['n','p','x','tn','tx']:
		seq_type_in = input('''Which program to use:
\tblastn (n), blastp(p), blastx(x), tblastn(tn), tblastx(tx)
Your decision: ''')
		if seq_type_in not in ['n','p','x','tn','tx']:
			print("Don't understand, please re-enter.\n")
	seq_type = seq_type_in
	return(seq_type)

		
def input_seqs(seq_type):
	seq_titles = []
		# This block is to find whether to perform search with existing file
	last_blast = "n"
	input_seq_by_file = "n"

	if pathlib.Path(input_seq_file).is_file() and stat(input_seq_file).st_size != 0:
		input_seq_by_file = None
		while input_seq_by_file != 'y' and input_seq_by_file != 'n':
			input_seq_by_file = input(f'\nFound input file, use as input?(y/n)')
	if input_seq_by_file == 'n' and pathlib.Path(temp_input_seq_file).is_file() and stat(temp_input_seq_file).st_size != 0:
		last_blast = None
		while last_blast != 'y' and last_blast != 'n':
			last_blast = input(f'Found last search in, use as input?(y/n)')
	
	input_sequences = []
	if input_seq_by_file == 'n' and last_blast == 'n':
		print('\nPlease enter/paste your gene(s)/protein(s)\n\tFasta header is recommended.\n\tUse Fasta headers for multiple sequences.\nTo finish input: Start a new line, press Ctrl-Z (windows) OR Ctrl-D(OSX/Linux), press <Enter> to save it.\nInput:')
		while True:
			try:
				line = input("")
			except EOFError:
				break
			input_sequences.append(line)
	elif input_seq_by_file == 'n' and last_blast == 'y':
		with open(temp_input_seq_file, 'r') as f:
			input_sequences = f.readlines()
	elif input_seq_by_file == 'y':
		with open(input_seq_file, 'r') as f:
			input_sequences = f.readlines()
		print(f"\nCopying sequence(s) from <{input_seq_file}> to <{temp_input_seq_file}>\n")
		
	if not input_sequences[0].startswith('>'):
		input_sequences.insert(0, f">Input_Sequence{str(time()).split('.')[0]}")
	num_seqs = 0
	with open(temp_input_seq_file, 'w') as input_file_handle:
		for line in input_sequences:
			if line.startswith('>'):
				input_file_handle.write(line.strip() + '\n')
				num_seqs += 1
				seq_titles.append(line[1:])
			else:
				input_file_handle.write(line.strip() + '\n')
	return temp_input_seq_file, num_seqs, seq_titles, seq_type

def doBlast(seq_type, temp_input_seq_file, output_path, file_name, num_seqs):
	
	# Define parameters for blast command
	num_threads = cpu_count()
	
	## Read from file part of parameters: 'nucldb', 'protdb', 'outfmt', 'evalue', 'word_size'
	parameterDict = dict(nucldb='nucl', protdb='prot', outfmt=None, evalue=None, word_size=None)
	with open(path.join(cwd, 'Parameters.csv'), 'r') as parametersCSV:
		paras = csvreader(parametersCSV)
		for i, row in enumerate(paras):
			if row[0] in ['nucldb', 'protdb', 'outfmt', 'evalue', 'word_size'] and row[1].strip() != '':
				parameterDict[row[0]] = row[1]
	nucldb    = parameterDict['nucldb'   ]
	protdb    = parameterDict['protdb'   ]
	outfmt    = parameterDict['outfmt'   ]
	evalue    = parameterDict['evalue'   ]
	word_size = parameterDict['word_size']
	## Define parameters: blast_program, 'db', 'out'
	dynamicNamePart = f'{file_name}{("_"+str(num_seqs)+"seqs") if num_seqs != 1 else ""}'
	blastProgramsAndOutput = {'p' :[NcbiblastpCommandline,  blastp_program_path , protdb, f'blast_protDB_{protdb}_{dynamicNamePart}.txt'],
							  'n' :[NcbiblastnCommandline,  blastn_program_path , nucldb, f'blast_nuclDB_{nucldb}_{dynamicNamePart}.txt'],
							  'x' :[NcbiblastxCommandline,  blastx_program_path , protdb, f'blast_ProtDB_{protdb}_transNucl_{dynamicNamePart}.txt'],
							  'tn':[NcbitblastnCommandline, tblastn_program_path, nucldb, f'blast_NuclDB_{nucldb}_transProt_{dynamicNamePart}.txt'],
							  'tx':[NcbitblastxCommandline, tblastx_program_path, nucldb, f'blast_NuclDB_{nucldb}_transNucl_{dynamicNamePart}.txt']
							 }
	blast_program = blastProgramsAndOutput[seq_type][1]
	db            = blastProgramsAndOutput[seq_type][2]
	out = path.join(output_path, blastProgramsAndOutput[seq_type][3])
	
	# Put arguments together
	arguments = {'cmd'        : blast_program,
				 'query'      : temp_input_seq_file,
				 'out'        : out,
				 'db'         : db, 
				 'num_threads': num_threads
				}
	if outfmt != None:
		arguments['outfmt'] = outfmt
	if evalue != None:
		arguments['evalue'] = evalue
	if word_size != None:
		arguments['word_size'] = word_size

	# Run blast command line
	run_blast = blastProgramsAndOutput[seq_type][0](**arguments)
	print(f'\nRunning BLAST command:\n{run_blast}')
	stdout, stderr = run_blast()
	if stderr == '':
		print(f'\nSuccess!! Find your result here:\n{out}\n')
	try:
		Popen(([*fileReader, out]))
	except:
		pass

# Database path configuration
print('\nSetting up BLASTDB environment variable...')
if "BLASTDB" in environ:
	origionalBLASTDBenv = getenv('BLASTDB')
else:
	origionalBLASTDBenv = False

try:
	environ['BLASTDB'] = path.join(cwd, "GenomeProteomeDatabase")
	print('Success, environment variable BLASTDB will be restored in the end.\n')
except:
	print('BLASTDB environment variable set up failed.')
	input('<Enter> to exit')
	exit()

# Do blast, return path configuration at any time
try:
	if not path.exists(output_path):
		makedirs(output_path)
	temp_input_seq_file, num_seqs, seq_titles, seq_type = input_seqs(seq_type())
	file_name = re.sub('[^\w_.)( -]', '.', seq_titles[0])[:30]
	try:
		doBlast(seq_type, temp_input_seq_file, output_path, file_name, num_seqs)
	except Exception as e:
		print('No blast done!')
		print(e)
except Exception as e:
	print(e)
	input("<Enter> to exit")
finally:
	print('\nRestoring environment varialbe... ', end = '')
	if origionalBLASTDBenv:
		environ['BLASTDB'] = origionalBLASTDBenv
	else:
		environ.pop('BLASTDB')
	print('Restored!')
	input('<Enter> to exit')