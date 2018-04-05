from Bio import SeqIO
import re
import os
import platform
from time import strftime

import write_file

# Define path
if platform.system() == 'Darwin':
	genome = '/Users/durand.dc/Documents/works/Misc.files/Streptomyces_genomes/M145PlusStredbNotes.gb'
	output_path = '/Users/durand.dc/Documents/works/Resources/temp/genome_txt_search'
elif platform.system() == 'Linux':
	genome = '/mnt/d/WORKs/Misc.files/Streptomyces_genomes/M145PlusStredbNotes.gb'
	output_path = '/mnt/d/WORKs/Resources/temp/genome_txt_search'	
else: # 'Windows'
	genome = r'D:\WORKs\Misc.files\Streptomyces_genomes\M145PlusStredbNotes.gb'
	output_path = r'D:\WORKs\Resources\temp\genome_txt_search'



search_word = ''
found_num = 0

while search_word == '':
	search_word = input('What do you want to search? (case insensitive, can be regular expression)\n')
	print("="*60)
	out_file = os.path.join(output_path, f'M145_{write_file.clean_str(search_word)}.txt')
	out_tsv = os.path.join(output_path, f'M145_{write_file.clean_str(search_word)}.tsv')

time_stamp = f"{'='*6} {strftime('%X %d/%m/%Y %Z')} "
time_stamp = time_stamp +  '='*(60-len(time_stamp)) + '\n'
write_file.write(time_stamp, out_file)

my_seq = SeqIO.read(genome, 'genbank')
find = re.compile(search_word, re.IGNORECASE)
foundlist = {}
for feat in my_seq.features:
	found = False # Set found flag for each feature run, turns true if found and stop searching against other tags in the same feature
	if feat.type in ['misc_feature', 'regulatory'] :
		continue
	for tag in feat.qualifiers:
		if found or tag in ['translation', 'transl_table', 'codon_start']:
			pass
		else:
			for txt in feat.qualifiers[tag]:
				if find.search(txt) != None:
					output_str = f'{feat}\nMatch:\t{txt}\n{"-"*60}'
					found = True
					try:
						SCO_num = feat.qualifiers['locus_tag'][0]
					except Exception as e:
						print(e)
						print(feat)
						break
					foundlist[SCO_num] = feat.type
			if found:
				write_file.write(output_str, out_file)
				found_num += 1

write_file.write(f'\nFound {found_num} features with the word "{search_word}"', out_file)
write_file.write(f'\nFound features:\n{foundlist}', out_file)

with open(out_tsv, 'w') as csv_file:
	for key, value in foundlist.items():
		csv_file.write(f'{key}\t{value}\n')
		
print(f'\ntsv file: {out_tsv}\n')
write_file.write(f'\n{"="*60}', out_file)


write_file.open_in_notepad(out_file)
