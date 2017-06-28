from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def print_translation(nucl, line_width,three_lettered = False, analysis = False):
	line_num = len(nucl)//line_width + 1
	prot = nucl.translate()
	trilett = seq3(prot,custom_map={"*": "***"})
	print(prot)
	print(trilett) if three_lettered == True else None
	print()
	
	for i in range(line_num):
		nucl_loci = list(map(lambda x: x*line_width, [i, i+1]))
		nuclsli = nucl[nucl_loci[0]:nucl_loci[1]]
		print(nucl_loci[0]+1,'\t',
			  str(nuclsli).ljust(line_width,' '),
			  '\t',nucl_loci[0]+len(nuclsli))
		
		prot_loci = [a * int(line_width/3) for a in [i, i+1]]
		protsli = prot[prot_loci[0]:prot_loci[1]]
		print(prot_loci[0]+1,'\t',
			  str(protsli).replace('','  ').rstrip().ljust(line_width,' '),
			  '\t',prot_loci[0]+len(protsli))
		if three_lettered:
			trilettsli = trilett[nucl_loci[0]:nucl_loci[1]]
			print(prot_loci[0]+1,'\t',
				  trilettsli.ljust(line_width,' '),
				  '\t',prot_loci[0]+len(protsli))
		print()
	if analysis:
		print('Protein properties:\n')
		prot_str = str(prot).split('*')[0]
		prot_ana = ProteinAnalysis(prot_str)
		print(f'Molecular weight:\t{prot_ana.molecular_weight():,.2f}')
		print(f'Isoelectric point:\t{prot_ana.isoelectric_point():,.2f}')
		print(f'Instability index:\t{prot_ana.instability_index():,.2f}')
		print(f'Aromaticity:\t\t{prot_ana.aromaticity():,.4f}')
		

nucl = Seq('''ATGGACTACAAGGACCACGACGGCGACTACAAGGACCACGACATCGACTACAAGGACGACGACGACAAGTAGACCAGAGCCCGCTCACCCGGCCCCAGATTGCGGTTGAAGTCC'''.replace('\n', '').replace('\r', ''))

print(nucl)
print(f'Nucleotide length: {len(nucl)}')

print_translation(nucl, 60, analysis = True, three_lettered = True)		  

