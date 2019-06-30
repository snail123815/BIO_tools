from Bio import SeqIO
import os
import pathlib # for checking the existance of a file
import platform

# Define path
if platform.system() == 'Darwin':
	genomePath = '/Users/durand.dc/Documents/works/Resources/Resource_QL109/Genome_QL109/Genome_SeperateGenbankFiles'
	output_path = '/Users/durand.dc/Documents/works/Resources/Resource_QL109/GenomeSlices'
elif platform.system() == 'Linux':
	genomePath = '/mnt/d/WORKs/Resources/Resource_QL109/Genome_QL109/Genome_SeperateGenbankFiles'
	output_path = '/mnt/d/WORKs/Resources/Resource_QL109/GenomeSlices'	
else: # 'Windows'
	genomePath = r'D:\WORKs\Resources\Resource_QL109\Genome_QL109\Genome_SeperateGenbankFiles'
	output_path = r'D:\WORKs\Resources\Resource_QL109\GenomeSlices'
	
	
def user_input():
	def input_num(promote):
		try:
			num = int(input(promote))
			return num
		except ValueError:
			print('Please enter a number!\n')
			return input_num(promote)
	quest = input('Do you want 1 gene(a) or several genes in a row?(b):')
	while quest != "a" and quest != "b":
		quest = input('a or b please:')
	
	# Take input from user
	if quest == "a":
		sc = input_num('Scaffold number: ')
		target_start = input_num('ORF number: ')
		target_end = target_start
	else:
		sc = input_num('Scaffold number: ')
		target_start = input_num('ORF number of the start gene: ')
		target_end = input_num('ORF number of the end gene: ')
	flanking = input_num('Please enter a number for flanking:')
	return(sc, target_start, target_end, flanking)

# Clear the log file if exist
def clear_log(logfile):
	if pathlib.Path(logfile).is_file():
		log_handle = open(logfile, 'w')
		log_handle.close()
		write_log('Previous log cleaned\n')

	
def write_log(logstr, end = None):
	log_handle = open(logfile, 'a')
	print(logstr, end = end, file = log_handle)
	print(logstr, end = end)
	log_handle.close()

def slice(sc, target_start, target_end, flanking):
	# Convert numbers to strings needed in output
	sc_str = str(sc).zfill(2)
	target_start_str = str(target_start).zfill(4)
	target_end_str = str(target_end).zfill(4)
	flanking_str = f'f{flanking}' if flanking != 0 else ''

	# Creat a global vairable that can be accessed by write_log function
	global logfile
	logfile = os.path.join(output_path, f"{sc_str}_{target_start_str}-{target_end_str}{flanking_str}.log")
	clear_log(logfile)

	# Read the file
	sample = SeqIO.read(os.path.join(genomePath,f"sc{sc_str}.gb"),"genbank")
	write_log('The imported sequence has {0} features.\n'.format(len(sample.features)))
	
	# Check if the data is actually imported
	genes_start = sample.features[target_start].location.start
	genes_end = sample.features[target_end].location.end
	write_log(f'The gene(s) wanted from {genes_start} to {genes_end}.\n')
	
	# Set the desired start and end point on the imported sequence, deal with out of range exceptions, and creat a sub sample
	start = genes_start - flanking
	end = genes_end + flanking
	if start < 0:
		start = 0
	if end > len(sample):
		end= len(sample)
	sub_sample = sample[start:end]
	write_log(f'Target seqID:\tSC{sc_str}_ORF{target_start_str} to {target_end_str}\nFlanking:{flanking}\n')
	
	# Some extra genes might appear in the flanking regions that will be included in the final file, inform user
	genes_needed = target_end - target_start + 1
	genes_loaded = len(sub_sample.features)
	write_log(f'There are {genes_needed} genes needed, {genes_loaded - genes_needed} extra genes in the flanking reigon.\n')
	
	# There is no orf information for unknown protein in the file, so we need to determin the first protein's orf number in case of unknown protein in the sliced sequence.
	first_orfnum = target_start - len(sample[start:genes_start].features)
	
	write_log('################# Genes in the output file #################\n')
	for feat_num, feat in enumerate(sub_sample.features):
		write_log(f'No.{feat_num+1}:')
		strand = f'{feat.location.strand:+}'[0]
		write_log(f'\tStrand: {strand}\n')
		write_log(f'{feat.qualifiers["ID"][0]}')
		write_log(f"{feat.qualifiers['Name'][0]}")
		
	# Rename the sample, if the name is too long, shorten it
	LOCUS = f"{sc_str}_{target_start}+{genes_loaded}{flanking_str}"
	if len(LOCUS)<=16:
		sub_sample.name = LOCUS
	else:
		LOCUS = f"{sc_str}_{target_start}+{genes_loaded}{flanking_str[0]}"
		sub_sample.name = LOCUS
	
	output_file = os.path.join(output_path, f"{sc_str}_{target_start_str}-{target_end_str}{flanking_str}.gb")
	
	if pathlib.Path(output_file).is_file():
		caution = 'Caution file overwriten:\n'
	else:
		caution = ''
	
	SeqIO.write(sub_sample, output_file, "genbank")
	# Should return 1 if successful
	# Show ending information
	write_log(f"\n########## Output files ##########\n\n{caution}{output_file}\n{logfile}\n")

# Test output:
# slice(2,1234,1245, 100000)
slice(*user_input())

# pause the program

input('Press <ENTER> to continue.')

