from Bio.Blast.Applications import NcbiblastpCommandline
import os

protdb_path = r"D:\WORKsSoftware\COGdata2014\COG_2003-2014_pDB"

blastp_program_path = r"D:\Program Files\NCBI\blast-2.6.0+\bin\blastp.exe"
num_threads = os.cpu_count()
protein_seq_path = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\operon_prediction_files"
blastp_output_path = r"D:\WORKs\Resources\Resource_QL109\Genome_QL109\operon_prediction_files\blast_COG"

def run_blastp(sc):
	scnum = str(sc).zfill(2)
	seq_file = os.path.join(protein_seq_path, f'sc{scnum}.faa')
	blastp_output = os.path.join(blastp_output_path, f'COG_sc{scnum}.blast_xml')

	run_blastp_run = NcbiblastpCommandline(blastp_program_path,
										   query = seq_file, 
										   out = blastp_output,
										   outfmt = 5, 
										   db = protdb_path,
										   evalue = 1e-5,
										   num_threads = num_threads,
										   max_target_seqs = 1,
										   max_hsps = 1
										   )
	print(f'Running BLASTP command:\n{run_blastp_run}')
	stdout, stderr = run_blastp_run()
	if stderr == '':
		print(f'\nSuccess!!\n{blastp_output}\n')

def parse_blastp(sc):
	scnum = str(sc).zfill(2)
	output_tab_file = os.path.join(protein_seq_path, 'temp', f'COG_sc{scnum}.table_out')
	blastp_output = os.path.join(blastp_output_path, f'COG_sc{scnum}.blast_xml')
	from Bio.Blast import NCBIXML
	with open(blastp_output, 'r') as blastp_file_handle, open(output_tab_file,'w') as output_handle:
		print('QL109_ID\tHit_ID\tE-value', file = output_handle)
		blastp_records = NCBIXML.parse(blastp_file_handle)
		for idx, blastp_record in enumerate(blastp_records):
			ql109_id = blastp_record.query.split("|")[1]
			hit_id = '-' if blastp_record.alignments == [] else blastp_record.alignments[0].hit_def.split("|")[1]
			e_value = '-' if blastp_record.alignments == [] else blastp_record.alignments[0].hsps[0].expect
			write_str = f'{ql109_id}\t{hit_id}\t{e_value}'
			print(write_str, file = output_handle)

# for sc in range(4,14):
	# run_blastp(sc)
	# parse_blastp(sc)
# run_blastp(sc)
# parse_blastp(3)
