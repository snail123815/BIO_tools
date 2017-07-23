from Bio import SeqIO
import re
import os
from time import strftime

import write_file

genome = r"D:\WORKs\Resources\Streptomyces_genomes\M145.gb"
search_word = ''
found_num = 0

while search_word == '':
	search_word = input('What do you want to search? (case insensitive, can be regular expression)\n')
	print("="*60)
	out_file = f'D:/WORKs/Resources/temp/genome_txt_search/M145_{write_file.clean_str(search_word)}.txt'
	out_tsv = f'D:/WORKs/Resources/temp/genome_txt_search/M145_{write_file.clean_str(search_word)}.tsv'

time_stamp = f"{'='*6} {strftime('%X %d/%m/%Y %Z')} "
time_stamp = time_stamp +  '='*(60-len(time_stamp)) + '\n'
write_file.write(time_stamp, out_file)

my_seq = SeqIO.read(genome, 'genbank')
find = re.compile(search_word, re.IGNORECASE)
foundlist = {}
for feat in my_seq.features:
	found = False # Set found flag for each feature run, turns true if found and stop searching against other tags in the same feature
	if feat.type == 'misc_feature':
		continue
	for tag in feat.qualifiers:
		if found or tag in ['translation', 'transl_table', 'codon_start']:
			pass
		else:
			for txt in feat.qualifiers[tag]:
				if find.search(txt) != None:
					output_str = f'{feat}\nMatch:\t{txt}\n{"-"*60}'
					found = True
					SCO_num = feat.qualifiers['locus_tag'][0]
					foundlist[SCO_num] = feat.type
			if found:
				write_file.write(output_str, out_file)
				found_num += 1

write_file.write(f'\nFound {found_num} features with the word "{search_word}"', out_file)
write_file.write(f'\nFound features:\n{foundlist}', out_file)

with open(out_tsv, 'w') as csv_file:
	for key, value in foundlist.items():
		csv_file.write(f'{key}\t{value}\n')
		
write_file.write(f'\ntsv file: {out_tsv}\n')
write_file.write(f'\n{"="*60}', out_file)


write_file.open_in_notepad(out_file)
