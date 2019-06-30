## convert genbank to GFF formart

from Bio import SeqIO
from BCBio import GFF

with open('/Users/durand.dc/Documents/works/Resources/Resource_M145/M145.gb','r') as gbFile, open('/Users/durand.dc/Documents/works/Resources/Resource_M145/M145.gff','w') as gffFile:
		GFF.write(SeqIO.parse(gbFile,'genbank'), gffFile)
