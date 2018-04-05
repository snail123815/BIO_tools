import os, platform

cwd = os.path.dirname(os.path.realpath(__file__))
output_path = os.path.join(cwd, "BLAST_Results")

nucldb = "nt"
protdb = "prot"
outputFormat = 0
evalue = 1
num_threads = os.cpu_count()

if platform.system() in ['Darwin', 'Linux']:
	print('\nThis wrapper needs blast programs setted up in your PATH.')
	blastn_program_path  = 'blastn'
	blastp_program_path  = 'blastp'
	blastx_program_path  = 'blastx'
	tblastn_program_path = 'tblastn'
	tblastx_program_path = 'tblastx'
	if platform.system() == 'Linux':
		fileReader = ['gedit']
	else:
		fileReader = ['open', '-a', 'TextWrangler']
else: # 'Windows'
	blastn_program_path  = os.path.join(cwd, "blastProgram/blastn.exe")
	blastp_program_path  = os.path.join(cwd, "blastProgram/blastp.exe")
	blastx_program_path  = os.path.join(cwd, "blastProgram/blastx.exe")
	tblastn_program_path = os.path.join(cwd, "blastProgram/tblastn.exe")
	tblastx_program_path = os.path.join(cwd, "blastProgram/tblastx.exe")
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
	num_seqs = 0
	seq_titles = []
	temp_input_seq_file = os.path.join(output_path, 'temp_input_seq.fa')
	input_seq_file = os.path.join(cwd, "BLAST_INPUT_fasta.txt")
	print(f"\nPut your input <fasta> sequence in <{input_seq_file}>\nOR paste your sequence later (multiple sequences allowed).\n")
		# This block is to find whether to perform search with existing file
	import pathlib
	last_blast = "n"
	input_seq_by_file = "n"

	if pathlib.Path(input_seq_file).is_file() and os.stat(input_seq_file).st_size != 0:
		input_seq_by_file = None
		while input_seq_by_file != 'y' and input_seq_by_file != 'n':
			input_seq_by_file = input(f'Found input file <{input_seq_file}>\nUse as input?(y/n)')
	if input_seq_by_file == 'n' and pathlib.Path(temp_input_seq_file).is_file() and os.stat(temp_input_seq_file).st_size != 0:
		last_blast = None
		while last_blast != 'y' and last_blast != 'n':
			last_blast = input(f'Found last search in <{temp_input_seq_file}>\nUse as input?(y/n)')
	
	input_sequences = []

	if input_seq_by_file == 'n' and last_blast == 'n':
		print('\nPlease enter/paste your gene(s)/protein(s) in <fasta> format\nTo finish, start a new line, press Ctrl-Z (windows) OR Ctrl-D(OSX/Linux) then press <Enter> to save it.\nInput:')
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
		print(f"\nCopy sequence(s) from <{input_seq_file}> to <{temp_input_seq_file}>\n")
	
	
	with open(temp_input_seq_file, 'w') as input_file_handle:
		for line in input_sequences:
			if line.startswith('>'):
				input_file_handle.write(line.strip() + '\n')
				num_seqs += 1
				seq_titles.append(line[1:])
			else:
				input_file_handle.write(line.strip() + '\n')
	return temp_input_seq_file, num_seqs, seq_titles, seq_type

def doBlast(seq_type, temp_input_seq_file, output_path, file_name):	
	from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastpCommandline, NcbiblastxCommandline, NcbitblastnCommandline, NcbitblastxCommandline
	
	outputFile = f'{file_name}_{num_seqs if num_seqs != 1 else ""}{"seqs_"if num_seqs != 1 else ""}'
	blastProgramsAndOutput = {'p' :[NcbiblastpCommandline,  blastp_program_path , protdb, f'blastp_protein_{outputFile}.txt'],
							  'n' :[NcbiblastpCommandline,  blastn_program_path , nucldb, f'blastn_nucleotide_{outputFile}.txt'],
							  'x' :[NcbiblastxCommandline,  blastx_program_path , protdb, f'blastx_transNucl2Prot_{outputFile}.txt'],
							  'tn':[NcbitblastnCommandline, tblastn_program_path, nucldb, f'tblastn_transProt2Nucl_{outputFile}.txt'],
							  'tx':[NcbitblastxCommandline, tblastx_program_path, nucldb, f'tblastx_transNucl2transNucl_{outputFile}.txt']
							 }
	blast_program = blastProgramsAndOutput[seq_type][1]
	db            = blastProgramsAndOutput[seq_type][2]
	out = os.path.join(output_path, blastProgramsAndOutput[seq_type][3])
	arguments = {'cmd'        : blast_program,
				 'query'      : temp_input_seq_file,
				 'out'        : out,
				 'outfmt'     : outputFormat, 
				 'db'         : db, 
				 'remote'     : False,
				 'evalue'     : evalue,
				 'num_threads': num_threads
				}
	run_blast = blastProgramsAndOutput[seq_type][0](**arguments)
	print(f'\nRunning BLAST command:\n{run_blast}')
	stdout, stderr = run_blast()
	if stderr == '':
		print(f'\nSuccess!! Find your result here:\n{out}\n')
	from subprocess import Popen
	try:
		Popen(([*fileReader, out]))
	except:
		pass

# Database path configuration
print('\nSetting up BLASTDB environment variable...')
if "BLASTDB" in os.environ:
	origionalBLASTDBenv = os.getenv('BLASTDB')
else:
	origionalBLASTDBenv = False

try:
	os.environ['BLASTDB'] = os.path.join(cwd, "GenomeProteomeDatabase")
	print('Success, environment variable BLASTDB will be restored in the end.\n')
except:
	print('BLASTDB environment variable set up failed.')
	input('')
	from sys import exit
	exit()

# Do blast, return path configuration at any time
try:
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	temp_input_seq_file, num_seqs, seq_titles, seq_type = input_seqs(seq_type())
	import re
	file_name = re.sub('[^\w_.)( -]', '.', seq_titles[0])[:30]
	try:
		doBlast(seq_type, temp_input_seq_file, output_path, file_name)
	except Exception as e:
		print('No blast done!')
		print(e)
except Exception as e:
	print(e)
	input("")
finally:
	print('Restoring environment varialbe... ', end = '')
	if origionalBLASTDBenv:
		os.environ['BLASTDB'] = origionalBLASTDBenv
	else:
		os.environ.pop('BLASTDB')
	print('Restored!')
	input('')