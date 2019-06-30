from Bio import SeqIO
import os

def generate_fake_fna(sc):
	scnum = str(sc).zfill(2)
	gb_file_dir = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\gb_file"
	gb_file = os.path.join(gb_file_dir, f'sc{scnum}.gb')
	fna_dir = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\operon_prediction_files"
	fna_file = os.path.join(fna_dir, f'sc{scnum}.fna')

	output_handle = open(fna_file, 'w')
	sequence = SeqIO.read(gb_file, 'genbank')
	
	cds_num = 0
	for feat in sequence.features:
		cds_num += 1 if feat.type == 'CDS' else 0
	sequence.id = f'gi|{scnum}{str(cds_num).zfill(6)}|ref|LNBE010000{scnum}.1| {sequence.id}, complete sequence'
	# gi number is fake, ref should link to the correct molecule on RefSeq, yet the sequence is actually from ourselves, and may be different slightly from that uploaded to genbank (some gaps still in our sequence)
	SeqIO.write(sequence,output_handle, 'fasta')

	output_handle.close()
	
# generate_fake_fna(2)
for i in range(1,14):
	generate_fake_fna(i)