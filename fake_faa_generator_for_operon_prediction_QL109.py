from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import pathlib # for checking the existance of a file


def generate_faa(sc):

	scnum = str(sc).zfill(2)
	
	gb_file_dir = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\gb_file"
	gb_file = os.path.join(gb_file_dir, f'sc{scnum}.gb')
	sequence = SeqIO.read(gb_file, 'genbank')

	faa_dir = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\operon_prediction_files"
	faa_file = os.path.join(faa_dir, f'sc{scnum}.faa')
	output_handle = open(faa_file, 'w')
	
	for num, feat in enumerate(sequence.features):
		# if num == 3:
			# break
		if feat.type == 'CDS':
			feat_header = f'gi|{scnum}{feat.qualifiers["locus_tag"][0][6:]}|ref|LNBE01{scnum}{feat.qualifiers["locus_tag"][0][8:]}.1| {feat.qualifiers["product"][0]} [MBT76]'
			feat_seq = ''
			for i in range(0,len(feat.qualifiers['translation'][0]),70):
				feat_seq += (feat.qualifiers['translation'][0][i:i+70])
				feat_seq += ('\n')
			feat_str = f'>{feat_header}\n{feat_seq}'
			print(feat_str, file = output_handle)
		else:
			continue

	output_handle.close()
	
# generate_faa(2)
for i in range(1,14):
	generate_faa(i)