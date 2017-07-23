from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import FeatureLocation
import os

protein_seqs = []



MBT27 = SeqIO.parse(r'D:\WORKs\Resources\Streptomyces_genomes\MBT27\MBT27.gbk',"genbank")
output = r"D:\WORKs\Resources\Streptomyces_genomes\MBT27\MBT27_Protein.fa"



for contig_num, contig in enumerate(MBT27):

	for idx, feat in enumerate(contig.features):
		if feat.type == 'CDS':
			protein_id = feat.qualifiers['locus_tag'][0]
			product = feat.qualifiers['product'][0]
			translation = feat.qualifiers['translation'][0].replace('\n','')
			protein_seq = Seq(translation, generic_protein)
			protein_record = SeqRecord(protein_seq, id = protein_id, description = product)
			
			protein_seqs.append(protein_record)
		
		
		
with open(output, 'w') as file_handle:
	SeqIO.write(protein_seqs, file_handle, 'fasta')