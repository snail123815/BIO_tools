

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


readFile = SeqIO.read('Streptomyces_genomes/M145.gb', 'genbank')


proteinSeq = []
for feature in readFile.features:
    if feature.type == 'CDS':
        protein = SeqRecord(Seq(feature.qualifiers['translation'][0], IUPAC.protein),
                            id=feature.qualifiers['locus_tag'][0],
                            description=feature.qualifiers['product'][0])
        proteinSeq.append(protein)


with open('Streptomyces_genomes/M145Proteome.fasta', 'w') as output:
    SeqIO.write(proteinSeq, output, "fasta")
