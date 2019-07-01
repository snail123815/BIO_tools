import pandas as pd


def getMeanErr(df, colIDs, ids=None):
    '''
    getMeanErr(df, ids=None)
    return
    dictStatDf[strain] = {'mean': strainMean.copy(),
                          'err': strainErr.copy()}
    mean and standard error will be calculated
    one value = nan, two values and more = calculate

    colIDS = {'strain1': [column names related to strain1],
              'strain2': [column names related to strain2],}
    number of ids equal to df.index
    '''
    dictStatDf = {}
    if type(ids) == type(None):  # infer ids from index if not provided
        ids = df.index

    # assuming you have different strains,
    # which need to be calculated separately
    for strain in colIDs:
        cols = colIDs[strain]
        strainMean = df.loc[ids, cols].mean(axis=1)
        strainSampleNum = df.loc[ids, cols].count(axis=1)
        strainErr = pd.Series(stats.sem(df.loc[ids, cols], axis=1, nan_policy='omit'),
                              index=ids)
        strainMean.name = f'{strain}_mean'
        strainErr.name = f'{strain}_err'
        # one value = nan, two values = calculate
        dictStatDf[strain] = {'mean': strainMean.copy(),
                              'err': strainErr.copy()}
    return dictStatDf
# getMeanErr(df)
