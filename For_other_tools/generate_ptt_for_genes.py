from Bio import SeqIO
import os


def generatePtt(seqRec, outputHandle, onlyProteins=False):
    geneTable = {}
    code = '-'
    cog = '-'
    geneNum = 0
    for feat in seqRec.features:
        if feat.type == "gene":
            location = f'{feat.location.start}..{feat.location.end}'
            strand = f'{feat.strand:+}'[0]
            length = f'{len(feat)}'
            pid = '-'
            geneId = feat.qualifiers["locus_tag"][0]
            try:
                synonym = feat.qualifiers["gene_synonym"][0]
            except:
                synonym = '_'
            try:
                product = feat.qualifiers['gene'][0]
            except:
                product = '-'
            if geneId in geneTable:
                count = len([g for g in geneTable if g.startswith(geneId)])
                geneId = f'{geneId}_{count+1}'
            geneTable[geneId] = [location, strand, length,
                                 pid, geneId, synonym, code, cog, product]

        if feat.type == 'CDS':
            geneId = feat.qualifiers["locus_tag"][0]
            if geneId in geneTable:
                count = len([g for g in geneTable if g.startswith(geneId)])
                if count > 1:
                    geneId = f'{geneId}_{count}'
            else:
                print(f"Genome file error, no corresponding gene for protein {geneId}")
                continue
            pid = f'{feat.qualifiers["protein_id"][0]}'
            product = f'{feat.qualifiers["product"][0]}'
            geneTable[geneId][3] = pid
            geneTable[geneId][8] = product
        if feat.type in ['tRNA', 'rRNA', 'misc_RNA']:
            geneId = feat.qualifiers["locus_tag"][0]
            if geneId in geneTable:
                count = len([g for g in geneTable if g.startswith(geneId)])
                if count > 1:
                    geneId = f'{geneId}_{count}'
            else:
                print(f"Genome file error, no corresponding gene for tRNA {geneId}")
                continue
            if geneTable[geneId][8] == '-':
                try:
                    product = feat.qualifiers['product'][0]
                except:
                    product = '-'
                note = feat.qualifiers['note'][0]  # for coelicolor genome, contains anticodon info
                if 'anticodon' in note:
                    geneTable[geneId][8] = note  # use complete info
                elif product != '-':
                    geneTable[geneId][8] = product
                else:
                    geneTable[geneId][8] = note.split(';')[0].split(',')[0]

    num = 0
    for feat in seqRec.features:
        num += 1 if feat.type == 'gene' else 0
    header = f'{seqRec.id}, complete sequence - 1..{len(seqRec)}\n{num} {"proteins" if onlyProteins else "genes"}\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct'
    print(header, file=outputHandle)

    for geneId in geneTable:
        if not onlyProteins:
            pass
        else:
            if geneTable[geneId][3] == '-':
                continue
        print("\t".join(geneTable[geneId]), file=outputHandle)


if __name__ == "__main__":
    seqRec = SeqIO.read(
        '/Users/durand.dc/Documents/works/Resources/Resource_M145/M145.gb', 'genbank')
    with open('/Users/durand.dc/Documents/works/Resources/Resource_M145/M145_genes.ptt', 'w') as handle:
        generatePtt(seqRec, handle, onlyProteins=False)
    with open('/Users/durand.dc/Documents/works/Resources/Resource_M145/M145_proteins.ptt', 'w') as handle:
        generatePtt(seqRec, handle, onlyProteins=True)
