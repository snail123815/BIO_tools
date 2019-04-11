# calculate DNA extinction coefficient
# calculate corrected concentration
# https://doi.org/10.1016/j.bpc.2007.12.004

from Bio.Seq import Seq

nearestNeighborCoefficientA260 = {'AA': 27400, 'AT': 22800, 'AC': 21200, 'AG': 25000,
                                  'TA': 23400, 'TT': 16800, 'TC': 16200, 'TG': 19000,
                                  'CA': 21200, 'CT': 15200, 'CC': 14600, 'CG': 18000,
                                  'GA': 25200, 'GT': 20000, 'GC': 17600, 'GG': 21600, }
singleBaseCoefficientA260 = {'A': 15400, 'T': 8700, 'C': 7400, 'G': 11500}

modificationExCoef = {'FAM': 20960, 'TET': 16255, 'HEX': 31580}
modificationMw = {'FAM': 537.46, 'TET': 675.24, 'HEX': 744.13}


def ssExtinctionA260(seq):
    exCoef = 0
    noneEndingSeq = seq[1:-1]
    for i in range(len(seq) - 1):
        exCoef += nearestNeighborCoefficientA260[seq[i:i + 2]]
        if i != len(seq) - 2:
            exCoef -= singleBaseCoefficientA260[noneEndingSeq[i]]
    return exCoef
# ssExtinctionA260


def extinctionA260(seq, nType='ds', printE=True, modifications=[None, ]):
    """nType can be 'ss'/'ds'"""
    seq = seq.replace('\n', '').replace('\t', '').replace(' ', '')
    seq = seq.upper()
    exCoef = ssExtinctionA260(seq)
    if nType == 'ds':
        rcSeq = str(Seq(seq).reverse_complement())
        exCoef += ssExtinctionA260(rcSeq)
        fGC = (seq.count('G') + seq.count('C')) / len(seq)
        fAT = (seq.count('A') + seq.count('T')) / len(seq)
        hypochromicity = (0.059 * fGC) + (0.287 * fAT)
        exCoef = exCoef * (1 - hypochromicity)
    for mod in modifications:
        if mod == None:
            break
        exCoef += modificationExCoef[mod]
    if printE:
        print(f'\nA260 extinction coefficient:\n{exCoef}')
    return exCoef
# extinctionA260


def ssMolecularWeight(seq):
    """g/mole
    default assuming synthesized DNA with no 5' phosphate"""
    seq = seq.upper()
    mw = (seq.count('A') * 313.21) + (seq.count('T') * 304.2) + \
        (seq.count('C') * 289.18) + (seq.count('G') * 329.21) - 61.96
    return mw
# ssMolecularWeight


def molecularWeight(seq, nType='ds', phosphate=None, modifications=[None, ]):
    """g/mole
    phosphate can be None/'single'/'double'
    modifications = ['FAM', 'TET','HEX']
        One keyword means one modification
    """
    if nType == 'ss' and phosphate == 'double':
        print("single strand don't have double phosphate")
        return
    seq = seq.replace('\n', '').replace('\t', '').replace(' ', '')
    seq = seq.upper()
    mw = ssMolecularWeight(seq)
    if nType == 'ds':
        rcSeq = str(Seq(seq).reverse_complement())
        mw += ssMolecularWeight(rcSeq)
    if phosphate == 'single':
        mw += 79
    elif phosphate == 'double':
        mw += 159

    for mod in modifications:
        if mod == None:
            break
        mw += modificationMw[mod]

    print(f"\nMolecular weight:\n{mw:.2f}")
    return mw
# molecularWeight


def molarConcentration(seq, a260, nType='ds', l=1, modifications=[None, ], printE=True):
    """c = a260/(e260*l) assuming l = 1 cm
    molar = mole / L"""
    exCoef = extinctionA260(seq, nType=nType,
                            modifications=modifications,  printE=printE)
    molarC = a260 / (exCoef * l)
    return molarC
# molarConcentration


def massConcentration(seq, a260, nType='ds', phosphate=None, modifications=[None, ],):
    """return concentration in g/L or mg/mL or μg/μL
    phosphate can be None/'single'/'double'
    """
    molarC = molarConcentration(seq, a260, nType=nType,
                                modifications=modifications, printE=False)
    mw = molecularWeight(seq, nType=nType, phosphate=phosphate,
                         modifications=modifications)
    massC = molarC * mw
    return massC
# massConcentration


if __name__ == '__main__':
    # seq = input('Please type DNA sequence:')
    # a260 = float(input('A260:'))
    # nType = input('Double strand (ds) or Single strand (ss):')

    seq = '''
TCCttAAttAGAATGATactTTttcaCTtt
'''
    a260 = 1
    nType = 'ds'
    modifications = [None, ]
    modifications = ['FAM', 'FAM']

    print(f"\n{nType}DNA sequence:")
    print(seq)

    molarC = molarConcentration(seq, a260, nType, modifications=modifications)
    massC = massConcentration(seq, a260, nType, modifications=modifications)
    print(f"\nConcentration A260 = {a260:.2f}:\n{molarC:.3} M")
    print(f"{massC:.2f} g/L")
    print()
    print(f"{molarC*1000000:.3f} μM")
    print(f"{massC*1000:.2f} ng/μL")
    print()
