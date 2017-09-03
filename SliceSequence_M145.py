from Bio import SeqIO
import os
import pathlib # for checking the existance of a file


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
		target_start = input_num('SCO number: ')
		target_end = target_start
	else:
		target_start = input_num('SCO number of the start gene: ')
		target_end = input_num('SCO number of the end gene: ')
	flanking = input_num('Please enter a number for flanking:')
	return(target_start, target_end, flanking)

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

def slice_M145(target_start, target_end, flanking, output_path = r'D:\WORKs\Resources\Streptomyces_genomes\Slices'):
	# Convert numbers to strings needed in output
	target_start_str = str(target_start).zfill(4)
	target_end_str = str(target_end).zfill(4)
	flanking_str = f'f{flanking}' if flanking != 0 else ''

	# Creat a global vairable that can be accessed by write_log function
	global logfile
	logfile = os.path.join(output_path, f"SCO{target_start_str}-{target_end_str}{flanking_str}.log")
	clear_log(logfile)

	# Read the file
	# Check if the data is actually imported
	sample = SeqIO.read(r'D:\WORKs\Resources\Streptomyces_genomes\M145.gb',"genbank")
	write_log(f'The imported sequence has {len(sample.features)} features.\n')
	
	for int, feat in enumerate(sample.features):
		if feat.type == 'gene' and 'locus_tag' in feat.qualifiers:
			if feat.qualifiers['locus_tag'][0] == f'SCO{target_start_str}':
				genes_start = feat.location.start
			if feat.qualifiers['locus_tag'][0] == f'SCO{target_end_str}':
				genes_end = feat.location.end
	if not 'genes_start' in locals() and not 'genes_end' in locals():
		write_log(f'Did not find SCO{target_start_str} to SCO{target_end_str}')
		return None
	else:
		write_log(f'The gene(s) wanted locate(s) from {genes_start} to {genes_end}.\n')

	# Set the desired start and end point on the imported sequence, deal with out of range exceptions, and creat a sub sample
	start = genes_start - flanking
	end = genes_end + flanking
	if start < 0:
		start = 0
	if end > len(sample):
		end= len(sample)
	sub_sample = sample[start:end]
	write_log(f'Target seqID:\tSCO{target_start_str} to {target_end_str}\nFlanking:{flanking}\n')
	
	# Some extra genes might appear in the flanking regions that will be included in the final file, inform user
	genes_needed = target_end - target_start + 1
	genes_loaded = 0
	for feat in sub_sample.features:
		genes_loaded += 1 if feat.type == 'gene' and 'locus_tag' in feat.qualifiers else 0
	write_log(f'There are {genes_needed} genes needed, {genes_loaded - genes_needed} extra genes in the flanking reigon.\n')
	
	write_log('################# Proteins in the output file #################\n')
	num_prots = 0
	for feat in sub_sample.features:
		if feat.type == 'CDS':
			num_prots += 1
			write_log(f"No.{num_prots}\n{feat.qualifiers['locus_tag'][0]}\tID: {feat.qualifiers['protein_id'][0]}\n{feat.qualifiers['product'][0]}")
			strand = f'{feat.location.strand:+}'[0]
			write_log(f'\tStrand: {strand}\n')
		
	# Rename the sample, if the name is too long, shorten it
	LOCUS = f"{sub_sample.name}_SCO{target_start_str}-{target_end_str}f{flanking_str}"
	if len(LOCUS)<=16:
		sub_sample.name = LOCUS
	else:
		LOCUS = LOCUS[:13] + '...'
		sub_sample.name = LOCUS
	
	output_file = os.path.join(output_path, f"SCO{target_start_str}-{target_end_str}{flanking_str}.gb")
	
	if pathlib.Path(output_file).is_file():
		caution = 'Caution file overwriten:\n'
	else:
		caution = ''
	
	SeqIO.write(sub_sample, output_file, "genbank")
	# Should return 1 if successful
	# Show ending information
	write_log(f"\n########## Output files ##########\n\n{caution}{output_file}\n{logfile}\n")

# Test output:
# slice(1839,1850, 100000)
slice_M145(*user_input())

# pause the program

input('Press <ENTER> to continue.')