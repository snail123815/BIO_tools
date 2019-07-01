# now only works for SCO numbers, adjust before use
import matplotlib.pyplot as plt
import pandas as pd
import getMeanErr


def getIds(scoNumbers, index):
    # convert 4 digit numbers to SCO+4 numbers, then search index if it exists
    # return found [ids] and not found [notDetected]
    ids = []
    found = {}
    notDetected = []
    for num in scoNumbers:
        found[num] = False
    for ind in index:
        for num in scoNumbers:
            sco = f'SCO{str(num).zfill(4)}'
            if ind.startswith(sco):
                ids.append(ind)
                found[num] = True
    for num in found:
        if not found[num]:
            notDetected.append(f'SCO{str(num).zfill(4)}')
    ids.sort()
    notDetected.sort()
    return ids, notDetected
# getIds


def queryProteins(scoNumbers,
                  querySamples=colIDs.keys(),
                  export=False,
                  exportFile='queryProtein(s).xlsx',
                  logScale=False,
                  xtickLabelRotation=None,
                  plotTitle=None,
                  returnValue=False,
                  figsize=None,
                  doPlotting=True,
                  ):
    '''
    needed:
    colIDs = {
    'WT': [f'z{str(i).zfill(2)} 1' for i in range(1, 4)],
    '2A': [f'z{str(i).zfill(2)} 1' for i in range(4, 7)],
    '2I': [f'z{str(i).zfill(2)} 1' for i in range(7, 10)],
    '8B': [f'z{str(i).zfill(2)} 1' for i in range(10, 13)],
    '9A': [f'z{str(i).zfill(2)} 1' for i in range(13, 16)],
    '9B': [f'z{str(i).zfill(2)} 1' for i in range(16, 19)]
    }
    '''
    ids, notDetected = getIds(scoNumbers, dfData.index)
    if len(ids) == 0:
        print('None detected.')
        if returnValue:
            return numQuery, 0, pd.DataFrame()
        else:
            return

    numFound = len(ids)

    dictStatDf = getMeanErr(dfData, colIDs, ids=ids)
    fig, ax = plt.subplots(1, 1)
    gotDataAverage = pd.DataFrame()
    gotDataErrors = pd.DataFrame()

    for sample in querySamples:
        mean = dictStatDf[sample]['mean']
        mean.name = sample
        err = dictStatDf[sample]['err']
        err.name = sample  # error bar needs to have the same name as mean

        gotDataAverage = pd.concat((gotDataAverage, mean), axis=1, sort=False)
        gotDataErrors = pd.concat((gotDataErrors, err), axis=1, sort=False)

    plotQuery = gotDataAverage.plot.bar(logy=False,
                                        sort_columns=True,
                                        yerr=gotDataErrors,
                                        width=0.8,
                                        title=plotTitle,
                                        antialiased=True,
                                        figsize=figsize
                                        )

    if xtickLabelRotation == None:
        if figsize:
            xlabelNonRotationLim = figsize[0] / 1.6
        else:
            xlabelNonRotationLim = 5
        if numFound <= xlabelNonRotationLim:
            plotQuery.set_xticklabels(plotQuery.get_xticklabels(), rotation=0)
        else:
            plotQuery.set_xticklabels(plotQuery.get_xticklabels(),
                                      rotation=45, horizontalalignment='right')
    else:
        plotQuery.set_xticklabels(plotQuery.get_xticklabels(
        ), rotation=xtickLabelRotation, horizontalalignment='right')

    plt.show()

    # if returnValue:
    #     return numQuery, numFound, dataResStat.reindex(proteinIDs).copy()
