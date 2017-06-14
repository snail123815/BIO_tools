from Bio import SeqIO
import re

my_seq = SeqIO.read('D:/WORKs/BIOBIO/temp/CP000255.1.gb', 'genbank')
search_word = input('What do you want to search? (case insensitive, can be regular expression)\n')
print("="*60)
found_num = 0
find = re.compile(search_word, re.IGNORECASE)

for feat in my_seq.features:
	found = False # Set found flag for each feature run, turns true if found and stop searching against other tags in the same feature
	for tag in feat.qualifiers:
		if found or tag in ['translation', 'transl_table', 'codon_start']:
			pass
		else:
			for txt in feat.qualifiers[tag]:
				if find.search(txt) != None:
					print(feat,f'\nMatch: {txt}\n{"-"*60}')
					found = True
					found_num += 1
print(f'\nFound {found_num} features with the word "{search_word}".\n{"="*60}')