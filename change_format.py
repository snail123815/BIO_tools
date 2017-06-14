from Bio import SeqIO

in_path = r"D:\WORKs\Project_Sporulation\Genomes\dcw cluster.embl"
out_path = r"D:\WORKs\Project_Sporulation\Genomes\dcw cluster.gb"

with open(out_path, 'w') as out_file:
	seq_handle = SeqIO.read(in_path, 'embl')
	SeqIO.write(seq_handle, out_file, 'genbank')