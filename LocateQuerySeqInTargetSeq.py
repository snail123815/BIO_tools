import lzma
import pickle
from Bio.Seq import Seq
from Bio import SeqIO

# load M145 genome from pickle file
with lzma.open('/Users/durand.dc/Documents/works/Resources/Resource_M145/M145.gb.bioseq.pickle.lzma', 'rb') as handle:
    genome = pickle.load(handle)

# the sequence to query
seqFull = Seq('''
GACTGGTACCCTCCGTTATGACTGCTGAAGATG
'''.strip())

# restriction sites to remove from query sequence
allResSites = {
    'KpnI': Seq('GGTACC'),
    'HindIII': Seq('AAGCTT')
}


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
        return None
    else:
        pass
    e = s + len(site)
    return [s, e, strand]
# findSeq


def findResLoc(allResSites):
    '''Input a dictionary of restriction sites {'nameA': SeqA, 'nameB': SeqB}
    Return a dictionary of found sites {'nameA': [s,e,strand], 'nameB': [s,e,strand]}'''
    locResSites = {}
    for resSite in allResSites:
        seqLoc = findSeq(seqFull, allResSites[resSite])
        if seqLoc != None:
            locResSites[resSite] = seqLoc
            # search for a second hit:
            secondHit = findSeq(seqFull[seqLoc[1]:], allResSites[resSite])
            if secondHit != None:
                secondHit[:-1] = [loc + seqLoc[1] for loc in secondHit[:-1]]
                print(f'There is multiple hits of {allResSites[resSite]} in target seq:')
                print(seqFull[:100])
                print('Second hit location:')
                print(f'{secondHit[0]} to {secondHit[1]-1} on strand {secondHit[2]}')
    return locResSites
# findResLoc


locResSites = findResLoc(allResSites)


if len(locResSites) == 1:
    singleResLoc = locResSites[list(locResSites.keys())[0]]
    seqClean = seqFull[singleResLoc[1]:]
elif len(locResSites) > 1:
    print(f'{len(resSites)} sites')
    exit()
else:
    print('No restriction sites found')
    seqClean = seqFull


targetLoc = findSeq(genome.seq, seqClean)
if targetLoc == None:
    print('Did not find on target')
elif targetLoc[2] == 1:
    print(targetLoc[0] + 1)
    print(targetLoc[0] + 1, targetLoc[1], targetLoc[2])
else:
    print(f'{targetLoc[1]} C')
    print(targetLoc[0] + 1, targetLoc[1], targetLoc[2])
