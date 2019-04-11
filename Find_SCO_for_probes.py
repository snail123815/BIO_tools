import pandas as pd
import lzma
import pickle
from Bio.Seq import Seq

probeFile = '/Users/durand.dc/Downloads/GPL9417_sequences.txt'
probeSeqCol = 4

probeDf = pd.read_csv(probeFile, sep='\t', index_col=None, header=0)

# load M145 genome from pickle file
with lzma.open('/Users/durand.dc/Documents/works/Resources/Resource_M145/M145.gb.bioseq.pickle.lzma', 'rb') as handle:
    genome = pickle.load(handle)


def findSeq(seq, site):
    '''
    Return [start, end, strand] of the query sequence (site)
    '''
    s = seq.find(site)
    strand = 1
    if s == -1:
        s = seq.find(site.reverse_complement())
        strand = -1
    if s == -1:
        return None, None, None
    else:
        pass
    e = s + len(site)
    return [s, e, strand]
# findSeq


probeStart = pd.Series(index=probeDf.index, dtype='int32', name='GenomeStart')
probeEnd = pd.Series(index=probeDf.index, dtype='int32', name='GenomeEnd')
probeStrand = pd.Series(index=probeDf.index, dtype='object', name='GenomeStrand')


for idx, seq in probeDf.iloc[:, probeSeqCol].iteritems():
    if idx % 1000 == 0:
        print(idx, end='|', flush=True)
    probe = Seq(seq)
    start, end, strand = findSeq(genome.seq, probe)
    if start == None:
        continue
    probeStart[idx] = start
    probeEnd[idx] = end
    probeStrand[idx] = ('+' if strand == 1 else '-')

newDf = pd.concat((probeDf, probeStart, probeEnd, probeStrand), axis=1)
newDf.to_csv('/Users/durand.dc/Downloads/GPL9417_sequences.genomePos.txt', sep='\t')
