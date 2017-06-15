from Bio import SeqIO
import re
import os

genome_dir = r'D:\WORKs\Resources\Resource_QL109\Genome_QL109'
search_word = r''
found_num = 0
for file in os.listdir(genome_dir):
	if file.endswith('.embl') and file.startswith('sc') and len(file) == 9:
		file_path = os.path.join (genome_dir, file)
		my_seq = SeqIO.read(file_path, 'embl')
		if search_word == r'':
			search_word = input('What do you want to search? (case insensitive, can be regular expression)\n')
			print("="*60)
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
								print(feat,f'\nMatch:\t{txt}\nOn:\t{file}\n{"-"*60}')
								found = True
								found_num += 1
	else:
		pass
print(f'\nFound {found_num} features with the word "{search_word}".\n{"="*60}')