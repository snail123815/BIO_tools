from Bio import SeqIO
import os
from Bio.SeqFeature import FeatureLocation

M145 = SeqIO.read(r'D:\WORKs\Resources\Streptomyces_genomes\M145.gb',"genbank")
output = r"D:\WORKs\Resources\Streptomyces_genomes\M145_gene_table.tsv"

# Define global variable as theses information can be used in different cycles of loop
def clear():
	global locus_tag
	global protein_id
	global GI
	global GeneID
	global product
	global translation
	global location
	locus_tag = ''
	protein_id = ''
	GI = ''
	GeneID = ''
	product = ''
	translation = ''
	location = FeatureLocation(0,0,strand = 1)

clear()
# After each print to file statement, a cleariance of the clobal variables will be done.

with open(output, 'w') as file_handle:
	print('locus_tag\tprotein_id\tGI\tEntrez_GeneID\tproduct\ttranslation', file = file_handle)
	for idx, feat in enumerate(M145.features):
		if feat.type == 'gene': # Unlike "CDS", "gene" type should contain the full list of SCO numbers, including those pseudo genes.
			location = feat.location # location fixed for identify 'CDS' or 'misc_feature' is for the same gene.
			locus_tag = feat.qualifiers['locus_tag'][0]
			if 'db_xref' in feat.qualifiers:
				for db_xref in feat.qualifiers['db_xref']:
					if db_xref.startswith('GI:'):
						GI = db_xref.split(':')[1]
					elif db_xref.startswith('GeneID:'):
						GeneID = db_xref.split(':')[1]
		elif feat.type == 'CDS' and feat.location.start == location.start and feat.location.end == location.end and feat.location.strand == location.strand: # location == location doesn't seem to work.
			protein_id = feat.qualifiers['protein_id'][0]
			if 'db_xref' in feat.qualifiers:
				for db_xref in feat.qualifiers['db_xref']:
					if db_xref.startswith('GI:'):
						GI = db_xref.split(':')[1]
					elif db_xref.startswith('GeneID:'):
						GeneID = db_xref.split(':')[1]
			product = feat.qualifiers['product'][0]
			translation = feat.qualifiers['translation'][0].replace('\n','')
			print(f"{locus_tag}\t{protein_id}\t{GI}\t{GeneID}\t{product}\t{translation}", file = file_handle)
			clear()
		elif feat.type == 'misc_feature' and feat.location.start == location.start and feat.location.end == location.end and feat.location.strand == location.strand:
			if 'db_xref' in feat.qualifiers:
				for db_xref in feat.qualifiers['db_xref']:
					if db_xref.startswith('GI:'):
						GI = db_xref.split(':')[1]
					elif db_xref.startswith('GeneID:'):
						GeneID = db_xref.split(':')[1]
			product = feat.qualifiers['note'][0].split(';')[-1] # Grab the last item in note, which is normally a brief description of predicted function
			translation = ''
			print(f"{locus_tag}\t{protein_id}\t{GI}\t{GeneID}\t{product}\t{translation}", file = file_handle)
			clear()