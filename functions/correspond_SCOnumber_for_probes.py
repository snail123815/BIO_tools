import pandas as pd
import lzma
import pickle
from Bio.Seq import Seq
from multiprocessing import Pool
from queue import Queue
import numpy as np

probeFile = '/Users/durand.dc/Downloads/GPL9417_sequences.txt'
probeSeqCol = 4

probeDf = pd.read_csv(probeFile, sep='\t', index_col=None, header=0)

# load M145 genome from pickle file
with lzma.open('/Users/durand.dc/Documents/works/Resources/Resource_M145/M145.gb.bioseq.pickle.lzma', 'rb') as handle:
    genome = pickle.load(handle)


def findSeq(seq, site, idx):
    '''
    Return [start, end, strand] of the query sequence (site)
    '''
    if idx % 1000 == 0:
        print(idx, end='|', flush=True)

    s = seq.find(site)
    strand = 1
    if s == -1:
        s = seq.find(site.reverse_complement())
        strand = -1
    if s == -1:
        return (idx, (np.nan, np.nan, ''))
    else:
        pass
    e = s + len(site)
    return (idx, (s, e, strand))
# findSeq


def main():
    probeStart = pd.Series(index=probeDf.index, dtype='int32', name='GenomeStart')
    probeEnd = pd.Series(index=probeDf.index, dtype='int32', name='GenomeEnd')
    probeStrand = pd.Series(index=probeDf.index, dtype='object', name='GenomeStrand')

    pool = Pool()
    print('Put everything in...')
    results = []
    for idx, seq in probeDf.iloc[:, probeSeqCol].iteritems():
        probe = Seq(seq)
        results.append(pool.apply_async(findSeq, (genome.seq, probe, idx)))
    pool.close()
    pool.join()

    for res in results:
        idx, (start, end, strand) = res.get()
        probeStart[idx] = start
        probeEnd[idx] = end
        probeStrand[idx] = strand

    newDf = pd.concat((probeDf, probeStart, probeEnd, probeStrand), axis=1)
    newDf.to_csv('/Users/durand.dc/Downloads/GPL9417_sequences.genomePos.txt', sep='\t')
    print('\nDone\n')


if __name__ == "__main__":
    main()
