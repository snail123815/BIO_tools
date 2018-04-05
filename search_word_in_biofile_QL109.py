from Bio import SeqIO
import re
import os
from time import strftime

import write_file

# Define path
import platform
if platform.system() == 'Darwin':
	workPath = '/Users/durand.dc/Documents/works'
	fileReader = ['open', '-a', 'TextWrangler']
elif platform.system() == 'Linux':
	workPath = '/mnt/d/WORKs'
else: # 'Windows'
	workPath = 'D:/WORKs'
	fileReader = ['notepad']



genome_dir = os.path.join(workPath, 'Resources/Resource_QL109/Genome_QL109')
search_word = ''
found_num = 0
found_proteins = []

while search_word == '':
	search_word = input('What do you want to search? (case insensitive, can be regular expression)\n')
	print("="*60)
	genome_name = genome_dir.split("/")[-1]
	out_file = os.path.join(workPath, f'Resources/temp/genome_txt_search/{genome_name}_{write_file.clean_str(search_word)}.txt')


time_stamp = f"{'='*6} {strftime('%X %d/%m/%Y %Z')} "
time_stamp = time_stamp +  '='*(60-len(time_stamp)) + '\n'
write_file.write(time_stamp, out_file)

for file in os.listdir(genome_dir):
	if file.endswith('.embl') and file.startswith('sc') and len(file) == 9:
		sc = file[2:4].zfill(3)
		file_path = os.path.join (genome_dir, file)
		my_seq = SeqIO.read(file_path, 'embl')
		find = re.compile(search_word, re.IGNORECASE)
		for feat in my_seq.features:
			found = False # Set found flag for each feature run, turns true if found and stop searching against other tags in the same feature
			for tag in feat.qualifiers:
				if found or tag in ['translation', 'transl_table', 'codon_start']:
					pass
				else:
					for txt in feat.qualifiers[tag]:
						if find.search(txt) != None:
							output_str = f'{feat}\nMatch:\t{txt}\nOn:\t{file}\n{"-"*60}'
							found = True
							txt_str = txt
					if found and 'translation' in feat.qualifiers:
						for note in feat.qualifiers['note']:
							if note.startswith('ORF number'):
								orf = note.split(' ')[2]
						protein_id = f'scaffold_{sc}_ORF_{orf}'
						found_proteins.append(protein_id)
						output_str = f'{feat}\nMatch:\t{txt_str}\nID:\t{protein_id}\n{"-"*60}'
					if found:
						write_file.write(output_str, out_file)
						found_num += 1
	else:
		pass

write_file.write(f'\nFound {found_num} features with the word "{search_word}", protein ids:', out_file)
for index, id in enumerate(found_proteins):
	output_str = f'{index+1}:\t{id}'
	write_file.write(output_str, out_file, end = '\n')
write_file.write(f'\n{"="*60}', out_file)
write_file.open_in_notepad(out_file)