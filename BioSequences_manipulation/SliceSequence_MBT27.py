from Bio import SeqIO

MBT27 = SeqIO.parse('./MBT27.gbk',"genbank")

protein_seqs = []

def user_input():
	def input_num(promote):
		try:
			num = int(input(promote))
			return num
		except ValueError:
			print('Please enter a number!\n')
			return input_num(promote)
	quest = input('Do you want 1 gene(a) or several genes in a row?(b):')
	while quest != "a" and quest != "b":
		quest = input('a or b please:')
	
	# Take input from user
	if quest == "a":
		target_start = input_num('Prokka number: ')
		target_end = target_start
	else:
		target_start = input_num('Prokka number of the start gene: ')
		target_end = input_num('Prokka number of the end gene: ')
	flanking = input_num('Please enter a number for flanking:')
	return(str(target_start).zfill(5), str(target_end).zfill(5), flanking)

prokkaNumFrom, prokkaNumTo, flanking = user_input()

print(f'From gene Prokka_{prokkaNumFrom} to {prokkaNumTo} with {flanking}bp flanking:')

for contig, seq in enumerate(MBT27):
	start = None
	end = None
	for feat in seq.features:
		if feat.type == 'CDS':
			if feat.qualifiers['locus_tag'][0] == f'Prokka_{prokkaNumFrom}':
				start = feat.location.start
			if feat.qualifiers['locus_tag'][0] == f'Prokka_{prokkaNumTo}':
				end = feat.location.end
	if start == None and end == None: # means nothing is found
		continue
	elif start == None or end == None: # means one gene (start) is found in this contig while the other (end) not.
		print('Please make sure the required genes are in the same contig!\nNothing created, run this program again.')
		break
	else:
		print(f'Found in contig No. {contig+1}')
		slicedSeq = seq[start-flanking:end+flanking]
		output = f'./SlicedSeqs/MBT27contig{contig}-{prokkaNumFrom}-{prokkaNumTo}f{flanking}.gbk'
		SeqIO.write(slicedSeq, output, 'genbank')
		print(f'File\n{output}\ncreated!')

input('Finished! Press <ENTER> to exit.')