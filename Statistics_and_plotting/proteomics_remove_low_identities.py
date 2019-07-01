'''This is for the pre processing of results from ISOQuant,
identThresh = one protein is identified in % of all samples,
intenSumThresh = total sum of the expression values'''


def removeLowIdentify(df, identThresh=0.7, intenSumThresh=1e5):
    numSamples = len(df.columns)
    print(f'Remove protein identified less than {identThresh:0.0%} samples:')
    print(f'Total identification: {len(df.index)}')
    identThresh = identThresh * numSamples
    totalIdentifyFilter = (df.count(axis=1) >= identThresh)
    print(f'{totalIdentifyFilter.sum()} with more than {int(np.ceil(identThresh))}/{len(df.columns)} identification ')
    totalIntensitySaver = (~totalIdentifyFilter & (df.sum(axis=1) >= intenSumThresh))
    print(f'{totalIntensitySaver.sum()} rescued with fewer identification but high quantification (sum >= {intenSumThresh:.0})\n')
    cleanDf = df.loc[(totalIdentifyFilter | totalIntensitySaver), :]
    return cleanDf
# removeLowIdentify
