from Bio import SeqIO
import os
from Bio.SeqFeature import SeqFeature, FeatureLocation
from math import ceil
import sys

try:
	genomePath = ' '.join(sys.argv[1:sys.argv.index('from')])

	# location transform to programming location
	start = int(sys.argv[sys.argv.index('from') + 1]) - 1
	end = int(sys.argv[sys.argv.index('to') + 1])
	if end <= start:
		raise
	if 'format' in sys.argv:
		format = sys.argv[sys.argv.index('format') + 1]
	else:
		format = 'genbank'
	if 'includePartial' in sys.argv:
		includePartialFeature = int(sys.argv[sys.argv.index('includePartial') + 1])
		if includePartialFeature not in [0,1]:
			raise
	else:
		includePartialFeature = 1
except:
	print("\nInput error.\n\
Usage:\npython3 SliceSequence_general.py [file path] from [start] to [end] {format [format] includePartial [1 or 0]}\n\
\nIf you want to split files, makeup 'from' and 'to' numbers, they will be ignored.\n")
	exit()

def doSlicing(genomePath, start, end, format, includePartialFeature):
	if format == 'genbank':
		ext = '.gb'
	else:
		ext = f'.{format}'
	sourceFile = os.path.basename(genomePath)
	sourceExt = os.path.splitext(genomePath)[1]
	if sourceExt in ['.gb', '.genbank','.gbk']:
		sourceFormat = 'genbank'
	elif sourceExt in ['.fa', '.fasta', '.fna']:
		sourceFormat = 'fasta'
		format = 'fasta'
	else:
		print('Unknown source file or format, check source:\n')
		print(f"path and extension = \n{os.path.splitext(genomePath)}")
		print()
		exit()
	outputPath = os.path.join(os.path.dirname(genomePath), f'{sourceFile}_slice_{start+1}_to_{end}{ext}')
	print('\nReading File...')
	source = list(SeqIO.parse(genomePath, sourceFormat))
	if len(source) > 1:
		seperate = None
		while seperate not in ('y', 'n'):
			seperate = input('More than one record found in source file.\n\
	Do you want to seperate records into seperate files? (y/n)')
		if seperate == 'n':
			exit()
		else:
			outputPath = os.path.join(os.path.dirname(genomePath), f'{sourceFile} parts')
			os.mkdir(outputPath)
			baseName = os.path.splitext(os.path.basename(genomePath))[0]
			print(f'Writing seperate files to\n{outputPath}/...')
			for idx, singleSeq in enumerate(source):
				singleSeqFilePath = os.path.join(outputPath, f'{baseName}_{idx+1}{sourceExt}')
				SeqIO.write(singleSeq, singleSeqFilePath, sourceFormat)
			print('\nDone!\n')
			exit()
	else:
		sourceSeq = source[0]
	targetSeq = sourceSeq[start:end]
	if includePartialFeature == 1 and sourceFormat == 'genbank':
		targetSeq.features.extend(findSpan(sourceSeq, start, end))
	print(f'\nWriting file:\n{outputPath}')
	SeqIO.write(targetSeq, outputPath, format)
	print('\nDone!\n')
# doSlicing

def findSpan(sourceSeq, start, end):
	print('\nLooking for partial features...')

	spanFeats = []
	for feat in sourceSeq.features:
		spanStart = start in feat
		spanEnd   = end   in feat
		spanFeat  = SeqFeature(type = feat.type)

		if spanStart and spanEnd:
			# Target position is inside feature
			if feat.type == 'CDS':
				newStart = (3 - abs(start - feat.location.start)) % 3
				newEnd   = end - start - abs(end - feat.location.start) % 3
			else:
				newStart = 0
				newEnd = end - start
			spanFeat.location = FeatureLocation(newStart,
			 									newEnd,
												strand = feat.location.strand)
			for key in feat.qualifiers:
				if key in ['gene_synonym', 'locus_tag', 'product', 'protein_id', 'db_xref', 'mol_type']:
					spanFeat.qualifiers[key] = [f'{keyStr} (slice)' for keyStr in feat.qualifiers[key]]
				elif key == 'translation':
					spanFeat.qualifiers[key] = list(feat.qualifiers[key])
					if feat.location.strand == 1:
						cutPointA = ceil((start - feat.location.start)/3)
						cutPointB = (end - feat.location.start) // 3
					else:
						len(spanFeat.qualifiers[key][0])
						cutPointA = len(spanFeat.qualifiers[key][0]) + 1 - (end - feat.location.start) // 3
						cutPointB = len(spanFeat.qualifiers[key][0]) + 1 - ceil((start - feat.location.start)/3)
					spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][cutPointA:cutPointB]
				else:
					spanFeat.qualifiers[key] = feat.qualifiers[key]
			spanFeats.append(spanFeat)

		elif spanStart:
			# Start position inside feature, feature ends in this range
			if feat.type == 'CDS':
				newStart = (3 - abs(start - feat.location.start)) % 3
			else:
				newStart = 0
			newEnd = feat.location.end - start
			spanFeat.location = FeatureLocation(newStart,
			 									newEnd,
												strand = feat.location.strand)
			for key in feat.qualifiers:
				if key in ['gene_synonym', 'locus_tag', 'product', 'protein_id', 'db_xref', 'mol_type']:
					spanFeat.qualifiers[key] = [f'{keyStr} (right part)' for keyStr in feat.qualifiers[key]]
				elif key == 'translation':
					spanFeat.qualifiers[key] = [keyStr for keyStr in feat.qualifiers[key]]
					if feat.location.strand == 1:
						cutPoint = len(spanFeat.qualifiers[key][0]) - ceil(len(spanFeat) / 3) + 1
						spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][cutPoint:]
					else:
						cutPoint = len(spanFeat) // 3
						spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][:cutPoint]
				else:
					spanFeat.qualifiers[key] = feat.qualifiers[key]
			spanFeats.append(spanFeat)

		elif spanEnd:
			# End position inside feature, feature ends in this range
			if feat.type == 'CDS':
				newEnd   = end - start - abs(end - feat.location.start) % 3
			else:
				newEnd = end - start
			newStart = feat.location.start - start
			spanFeat.location = FeatureLocation(newStart,
			 									newEnd,
												strand = feat.location.strand)
			for key in feat.qualifiers:
				if key in ['gene_synonym', 'locus_tag', 'product', 'protein_id', 'db_xref', 'mol_type']:
					spanFeat.qualifiers[key] = [f'{keyStr} (left part)' for keyStr in feat.qualifiers[key]]
				elif key == 'translation':
					spanFeat.qualifiers[key] = list(feat.qualifiers[key])
					if feat.location.strand == 1:
						cutPoint = len(spanFeat) // 3
						spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][:cutPoint]
					else:
						cutPoint = ceil((len(feat) - len(spanFeat)) / 3)
						spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][cutPoint:]
				else:
					spanFeat.qualifiers[key] = feat.qualifiers[key]
			spanFeats.append(spanFeat)

		else:
			# Not in range, ignore
			continue

	return spanFeats
# findSpan

doSlicing(genomePath, start, end, format, includePartialFeature)
