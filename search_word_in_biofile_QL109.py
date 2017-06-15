from Bio import SeqIO
import re
import os

import write_file

genome_dir = r'D:\WORKs\Resources\Resource_QL109\Genome_QL109'
search_word = None
found_num = 0

for file in os.listdir(genome_dir):
	if file.endswith('.embl') and file.startswith('sc') and len(file) == 9:
		file_path = os.path.join (genome_dir, file)
		my_seq = SeqIO.read(file_path, 'embl')
		if search_word == None:
			search_word = input('What do you want to search? (case insensitive, can be regular expression)\n')
			print("="*60)
			genome_name = genome_dir.split("\\")[-1]
			out_file = f'D:/WORKs/BIOBIO/temp/{genome_name}_{write_file.clean_str(search_word)}.txt'
		else:
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
								found_num += 1
								write_file.write(output_str, out_file)
	else:
		pass

write_file.write(f'\nFound {found_num} features with the word "{search_word}".\n{"="*60}', out_file)
write_file.open_in_notepad(out_file)