from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from io import StringIO


def parse_GFF(f):
    """Parse Gff annotation into list of seqRecords, only accept gff-version 3

    Args:
        f (str): gff file path

    Returns:
        list: list of SeqRecords
    """
    recs = OrderedDict()
    fasta = False
    # switch on if fasta header is found, then we stop parsing the line as gff annotation
    # Prokka output gffs usually include the sequence
    fastring = ''
    with open(f) as fh:
        firstLine = fh.readline()
        assert firstLine.lower().startswith('##gff-version 3'), 'Only accept gff version 3'
        for l in fh.readlines():
            # Process header
            if l.startswith('##'):
                if l.lower().startswith('##fasta'):
                    fasta = True
                else:
                    continue
            if fasta:
                fastring += l
                continue
            row = l.split('\t')

            # Create SeqRecord
            recid = row[0]
            if recid in recs:
                pass
            else:
                recs[recid] = SeqRecord('', id=recid, features=[], annotations={"molecule_type": ""})
                # molecule type is needed for writing, since we don't know the type, keep it as empty string
            
            # Create feature
            feat = SeqFeature(
                FeatureLocation(int(row[3]), int(row[4])),
                type=row[2],
                strand=(1 if row[6] == "+" else -1 if row[6] == "-" else 0)
            )
            # Parse qualifiers
            qualifiers = row[-1].split(';')
            for q in qualifiers:
                if '=' not in q:
                    continue
                qname, qvalue = q.split('=')
                qvalue = qvalue.strip() # in case qvalue is in the end of the line
                if qname not in feat.qualifiers:
                    feat.qualifiers[qname] = []
                feat.qualifiers[qname].append(qvalue)
            recs[recid].features.append(feat)
    
    # Write fasta if have
    if fasta:
        fastaFile = StringIO(fastring)
        for s in SeqIO.parse(fastaFile, 'fasta'):
            recs[s.id].seq = s.seq
    else:
        for id in recs:
            recs[id].seq = Seq('')
            
    return list(recs.values())
