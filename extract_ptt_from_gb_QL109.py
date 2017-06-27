from Bio import SeqIO
import os
import pathlib # for checking the existance of a file


def generate_ptt(sc):
	scnum = str(sc).zfill(2)
	gb_file_dir = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\gb_file"
	gb_file = os.path.join(gb_file_dir, f'sc{scnum}.gb')
	ptt_dir = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\operon_prediction_files"
	ptt_file = os.path.join(ptt_dir, f'sc{scnum}.ptt')

	output_handle = open(ptt_file, 'w')
	sequence = SeqIO.read(gb_file, 'genbank')

	cds_num = 0
	for feat in sequence.features:
		cds_num += 1 if feat.type == 'CDS' else 0

	header_str = f'{sequence.id}, complete sequence - 1..{len(sequence)}\n{cds_num} proteins\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct'
	print(header_str,file = output_handle)

	for num, feat in enumerate(sequence.features):
		# if num == 3:
			# break
		if feat.type == 'CDS':
			Location = f'{feat.location.start}..{feat.location.end}'
			Strand = f'{feat.strand:+}'[0]
			Length = f'{len(feat)}'
			PID = f'{feat.qualifiers["protein_id"][0]}'
			Gene = '-'
			Synonym = f'{feat.qualifiers["locus_tag"][0]}'
			Code = '-'
			COG = '-'
			Product = f'{feat.qualifiers["product"][0]}'
			feat_str = f'{Location}\t{Strand}\t{Length}\t{PID}\t{Gene}\t{Synonym}\t{Code}\t{COG}\t{Product}'
			print(feat_str, file = output_handle)
		else:
			continue

	output_handle.close()
	
	
for i in range(1,14):
	generate_ptt(i)