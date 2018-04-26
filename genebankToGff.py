## convert genbank to GFF formart

from Bio import SeqIO
from BCBio import GFF

with open('/mnt/d/WORKs/Resources/Resource_QL109/Genome_QL109/QL109_all_in_1.gbk','r') as ql109gb, open('/mnt/d/WORKs/Resources/Resource_QL109/Genome_QL109/QL109_all_in_1.gff','w') as ql109gff:
		GFF.write(SeqIO.parse(ql109gb,'genbank'), ql109gff)