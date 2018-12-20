from Bio import ExPASy
import pandas as pd
from Bio import SeqIO
from Bio import SwissProt
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from io import StringIO
from math import ceil
from threading import BoundedSemaphore, Thread
from queue import Queue
from time import time
from os import cpu_count, path, environ


def writeProtUniprot(ids, output):
    seqList = []
    notFound = []
    maxConnections = 5
    startTime = time()

    def getSeq(semaphore, id, seqQueue, notFoundIdQueue):
        with semaphore:
            print(f'Getting {id}...')
            try:
                with ExPASy.get_sprot_raw(id) as handle:
                    seqRec = SeqIO.read(handle, 'swiss')
                    seqQueue.put(seqRec)
            except:
                print(f'{id} not found')
                notFoundIdQueue.put(id)
    semaphore = BoundedSemaphore(value=maxConnections)
    seqQueue = Queue()
    notFoundIdQueue = Queue()
    threads = []
    for id in ids:
        get = Thread(target=getSeq, args=(semaphore, id, seqQueue, notFoundIdQueue))
        get.start()
        threads.append(get)
    for get in threads:
        get.join()
    while not seqQueue.empty():
        seqList.append(seqQueue.get())
    while not notFoundIdQueue.empty():
        notFound.append(notFoundIdQueue.get())
    endTime = time()
    print(endTime - startTime)
    SeqIO.write(seqList, output, 'fasta')
# writeProtUniprot


def setBlastDbEnv(dbDir='', origionalBLASTDBenv='', restore=False):
    if restore:
        print('\nRestoring BLASTDB environment variable...')
        if origionalBLASTDBenv:
            environ['BLASTDB'] = origionalBLASTDBenv
        else:
            environ.pop('BLASTDB')
        print('Restored!')

    else:
        print('\nSetting up BLASTDB environment variable...')
        if "BLASTDB" in environ:
            origionalBLASTDBenv = getenv('BLASTDB')
        else:
            origionalBLASTDBenv = False

        try:
            environ['BLASTDB'] = dbDir
            print('Success, environment variable BLASTDB will be restored in the end.\n')
        except:
            print('BLASTDB environment variable set up failed.')
            input('<Enter> to exit')
            exit()
        return origionalBLASTDBenv
# setBlastDbEnv


def findSco(queryFile, logFile):
    dictFound = {}
    blastpCline = NcbiblastpCommandline(query=queryFile,
                                        db='/Users/durand.dc/Documents/works/Resources/Resource_M145/BlastM145/GenomeProteomeDatabase/S.coelicolor',
                                        evalue=1e-10,
                                        outfmt=5,
                                        num_threads=cpu_count()
                                        )
    stdout, stderr = blastpCline()
    results = NCBIXML.parse(StringIO(stdout))
    log = open(logFile, 'w')
    for blastRec in results:
        scoId = ''
        uniprotId = blastRec.query.split()[0]
        if len(blastRec.alignments) == 0:
            dictFound[uniprotId] = ''
            continue
        for alignment in blastRec.alignments:
            # iterate to find scoId (low chance if first do not fit)
            if scoId == '':
                for hsp in alignment.hsps:
                    if 0.98 * blastRec.query_length < hsp.identities < 1.02 * blastRec.query_length:
                        pass
                    else:
                        # will lookfor the next hsp that meets the standard (probably none)
                        continue
                    scoId = '_'.join([text for text in alignment.title.split()
                                      if text.startswith('SCO') or text.startswith('SCP')])
                    if hsp.align_length / blastRec.query_length != 1:
                        print('*' * 100, file=log)
                        print(uniprotId, file=log)
                        print(scoId, file=log)
                        print(f'coverage: {hsp.align_length/blastRec.query_length:0.1%}',
                              file=log)
                        print(f'Query total {blastRec.query_letters}', file=log)
                        print(f'Hit: {hsp.query_start} - {hsp.query_end}',
                              file=log)
                        print(f'Target total {alignment.length}', file=log)
                        print(f'Hit: {hsp.sbjct_start} - {hsp.sbjct_end}',
                              file=log)
                        for x in range(ceil(hsp.align_length / 100)):
                            print(hsp.query[x * 100:(x + 1) * 100], file=log)
                            print(hsp.match[x * 100:(x + 1) * 100], file=log)
                            print(hsp.sbjct[x * 100:(x + 1) * 100], file=log)
                            print(file=log)
                    break  # only first hsp that meets the standard
        dictFound[uniprotId] = scoId
    log.close()
    return dictFound
# findSco


if __name__ == '__main__':
    output_unmappedUniprot = '/Users/durand.dc/Documents/works/Resources/Resource_M145/S.coelicolor_uniprot_unmapped_prot.fa'
    uniprot_coelicolor = '/Users/durand.dc/Documents/works/Resources/Resource_M145/S.coelicolor_uniprot.tsv'
    output_unmappedUniprot = '/Users/durand.dc/Documents/works/Resources/Resource_M145/S.coelicolor_uniprot_unmapped_prot.fa'
    dfUniprot = pd.read_csv(uniprot_coelicolor, delimiter='\t', index_col='Entry')

    haveScoOrScp = dfUniprot.loc[:, 'Gene names'].str.contains(
        'SCO') | dfUniprot.loc[:, 'Gene names'].str.contains('SCP')

    unmappedUniprotIds = dfUniprot[-haveScoOrScp].index

    # writeProtUniprot(unmappedUniprotIds, output_unmappedUniprot)
    logFindSco = '/Users/durand.dc/Documents/works/Resources/Resource_M145/S.coelicolor_uniprot_unmapped_prot_mapped.log'
    dictFound = findSco(output_unmappedUniprot, logFindSco)
    for uniprotId in unmappedUniprotIds:
        if dictFound[uniprotId] == '':
            continue
        newName = [dictFound[uniprotId], ]
        geneNameOld = dfUniprot.loc[uniprotId, 'Gene names']
        if type(geneNameOld) == str:
            newName.append(geneNameOld)
        newName = ' '.join(newName)
        dfUniprot.loc[uniprotId, 'Gene names'] = newName
    newUniprotTsv = '/Users/durand.dc/Documents/works/Resources/Resource_M145/S.coelicolor_uniprot_unmapped_prot_mapped.tsv'
    dfUniprot.to_csv(newUniprotTsv, sep='\t')
