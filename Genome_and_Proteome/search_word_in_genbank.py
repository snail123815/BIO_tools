description = """
Search a word or regular expression pattern in a genbank file,
output a table of match features and a text file for detail of the matches.

By DU Chao
c.du@biology.leidenuniv.nl
durand.dc@hotmail.com

"""

from argparse import ArgumentParser
import re
import sys
import subprocess
import platform
import argparse
from time import strftime, sleep
from tempfile import NamedTemporaryFile

from Bio import SeqIO

if platform.system() == 'Darwin':
	fileReader = ['open', '-a', 'TextEdit']
elif platform.system() == 'Linux':
	fileReader = ['cat']
else: # 'Windows'
	fileReader = ['notepad']

ArgParser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
ArgParser.add_argument('genome', type=str, help="genbank file to search")
ArgParser.add_argument('-nameTag', type=str, help="Tag of feature name, default 'locus_tag'", default='locus_tag')
args = ArgParser.parse_args()
genome = args.genome.strip()
nametag = args.nameTag.strip()

# read genome
my_seq = SeqIO.read(genome, 'genbank')

while True:
	search_word = ''
	found_num = 0
	while search_word == '':
		search_word = input('What do you want to search? (case insensitive, can be regular expression)\n')
		print("="*60)
	out_text = NamedTemporaryFile(mode='w')
	out_tsv = NamedTemporaryFile(mode='w')

	time_stamp = f"{'='*6} {strftime('%X %d/%m/%Y %Z')}"
	time_stamp += f"{'='*(60-len(time_stamp))}\n"
	out_text.write(time_stamp)

	find = re.compile(search_word, re.IGNORECASE)

	foundlist = set()

	for feat in my_seq.features:
		output_str = False
		#if feat.type == 'misc_feature':
		#	continue
		for tag in feat.qualifiers:
			if output_str or tag in ['translation', 'transl_table', 'codon_start']:
				continue
			else:
				for txt in feat.qualifiers[tag]:
					if find.search(txt) != None:
						output_str = f'{feat}\nMatch:\t{txt}\n{"-"*60}'
						try:
							featTag = feat.qualifiers[nametag][0]
						except KeyError as err:
							print(err, feat.qualifiers[tag])
							sys.exit()
						foundlist.add(featTag)
						found = True
						break
		if output_str:
			out_text.write(output_str)
			found_num += 1

	out_text.write(f'\nFound {found_num} features with the word "{search_word}"')
	out_text.write(f'\nFound features:\n{foundlist}')

	with open(out_tsv.name, 'w') as csv_file:
		for value in foundlist:
			csv_file.write(f'{value}\n')
			
	print(f'\ntsv file: {out_tsv.name}\n')
	out_text.write(f'\n{"="*60}')

	out_text.seek(0)
	out_tsv.seek(0)

	subprocess.run([*fileReader, out_tsv.name])
	sleep(2) # leave some time for app to open
	subprocess.run([*fileReader, out_text.name])
	input(f"copy at your will:\ntsv: {out_tsv.name}\ntext: {out_text.name}\n")
	out_text.close()
	out_tsv.close()
