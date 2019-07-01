import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import getMeanErr


# When you need to assign different colors to some data points:
'''
colorSpaceAb = {
    'Actinorhodin': [range(5071, 5093), 'blue', 'cornflowerblue'],
    # 'TW95a': [range(5314, 5321), 'black', 'dimgray'],
    'Prodiginines': [range(5877, 5899), 'red', 'lightcoral'],
    'CDA': [range(3210, 3250), 'magenta', 'violet'],
    # 'Coelichelin': [range(489, 500), 'dodgerblue', 'skyblue'],
    # 'Unknown NRPS': [range(6429, 6439), 'turquoise', 'paleturquoise'],
    'Coelimycin P1': [[6265] + list(range(6273, 6289)), 'orange', 'navajowhite'],
    # 'Geosmin': [range(6073, 6074), 'olive', 'darkkhaki'],
    # 'Hopanoids': [range(6759, 6772), 'blueviolet', 'thistle'],
    # 'Unknown Deoxysugar': [range(381, 402),  'orange', 'peachpuff'],
}


def isAbAsignColor(scoNum, isSig):
    found = False
    for ab in colorSpaceAb:
        if scoNum in colorSpaceAb[ab][0]:
            found = True
            if isSig:
                return colorSpaceAb[ab][1]
            else:
                return colorSpaceAb[ab][2]
    if not found:
        return None

'''


def plotVolcano(sample,
                plotLegend=False,
                xlims=None,
                ylims=None,
                plotShow=True,
                outputTable=False
                ):
    colorDict = {}
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    dataMeanErr = getMeanErr(dfData)
    wtMeanErr = dataMeanErr['WT']
    sampleMeanErr = dataMeanErr[sample]

    ttest2wt = stats.ttest_ind(dfData.loc[:, colIDs[sample]], dfData.loc[:, colIDs['WT']],
                               axis=1,
                               equal_var=False
                               )
    ttest2wt = pd.DataFrame({f'{sample}2wtTstat': ttest2wt[0], f'{sample}2wtTpval': ttest2wt[1]},
                            index=dfData.index
                            )

    xAxis = pd.Series(np.log2(sampleMeanErr['mean'] / wtMeanErr['mean']),
                      index=dfData.index)
    yAxis = pd.Series(-np.log10(ttest2wt[f'{sample}2wtTpval']),
                      index=dfData.index)
    data2plot = pd.concat((xAxis, yAxis), axis=1).replace(
        [np.inf, -np.inf], np.nan).dropna(axis=0, how='any')
    data2plot.columns = [0, 1]

    if outputTable:
        fileName = f'volcano_{sample}.xlsx'
        dataOut = pd.concat((data2plot,
                             dfData.loc[:, colIDs[sample]],
                             sampleMeanErr['mean'],
                             dfData.loc[:, colIDs['WT']],
                             wtMeanErr['mean']
                             ),
                            axis=1,
                            join='inner'
                            )
        dataOut.columns = ['log2_Fold', '-log10_p-Value',
                           f'{sample}_1', f'{sample}_2', f'{sample}_3', f'{sample}_mean',
                           'WT_1',   'WT_2',   'WT_3', 'WT_mean']
        dataOut.to_excel(fileName, sheet_name=f'{sample}_WT')

    colorDict[ax] = pd.Series(index=data2plot.index, dtype=np.dtype('U'))
    for index in colorDict[ax].index:
        scoNum = int(re.findall(r'[0-9]+', index)[0])
        c = None  # for adding function of specify color of data points
        if data2plot.loc[index, 1] < -np.log10(thresholdPvalue) or \
                -np.log2(thresholdFold) < data2plot.loc[index, 0] < np.log2(thresholdFold):
            # c = isAbAsignColor(scoNum, isSig=False)
            # colorDict[ax].at[index] = c
            if c == None:
                colorDict[ax].at[index] = 'lightgrey'
        else:
            # c = isAbAsignColor(scoNum, isSig=True)
            # colorDict[ax].at[index] = c
            if c == None:
                colorDict[ax].at[index] = 'grey'

    # Dot size calculation dependent on the actuall value
    pointSizeData = sampleMeanErr['mean'].loc[colorDict[ax].index]

    def calculateSize(x): return x / 4 * np.pi / 1000 + 1
    pointSize = calculateSize(pointSizeData)

    # PLOTTING
    ax.scatter(data2plot[0], data2plot[1], s=pointSize, marker='o',
               alpha=0.8, c=colorDict[ax], edgecolors='none')

    if xlims != None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims != None:
        ax.set_ylim(ylims[0], ylims[1])

    ax.set_title(f'{sample}')
    ax.set_xlabel('$\log_2$ fold change')
    ax.set_ylabel('-$\log_{10}$ p-value')
    ax.axhline(y=-np.log10(thresholdPvalue),
               linestyle='--',
               linewidth=0.3,
               color='gray'
               )
    ax.text(ax.get_xlim()[0], -np.log10(thresholdPvalue)
            + 0.1, '$\\mathit{p} = $' + str(thresholdPvalue))
    legendMaxMin = [200000, 10000]

    def fmt(x):
        a, b = f'{x:.1e}'.split('e')
        b = int(b)
        return r'${a} \times 10^{{{b}}}$'.format(a=a, b=b)

    legendElements = [Line2D([0], [0], marker='o', color='w',
                             markersize=np.sqrt(calculateSize(legendMaxMin[0])),
                             markeredgecolor='k',
                             markerfacecolor='w', alpha=1, linewidth=0,
                             label=fmt(legendMaxMin[0])
                             ),
                      Line2D([0], [0], marker='o', color='w',
                             markersize=np.sqrt(calculateSize(legendMaxMin[1])),
                             markeredgecolor='k',
                             markerfacecolor='w', alpha=1, linewidth=0,
                             label=fmt(legendMaxMin[1])
                             )
                      ]

    abColorLegends = []
    for ab in colorSpaceAb:
        newLeg = Line2D([0], [0], marker='o', color='w',
                        markersize=10,
                        markerfacecolor=colorSpaceAb[ab][1],
                        alpha=1, linewidth=0,
                        label=ab)
        abColorLegends.append(newLeg)
    if plotLegend:
        sizeLegend = plt.legend(handles=legendElements, title='Protein level',
                                labelspacing=1.3, loc='lower left')
        ax.add_artist(sizeLegend)
        colorLegend = ax.legend(handles=abColorLegends, title='Color map', loc='upper left')

    # set aspect ratio of the axes
    ratio = 1
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    # the abs method is used to make sure that all numbers are positive
    # because x and y axis of an axes maybe inversed.
    ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)

    fig.tight_layout()

    HoverClick = interactionHoverClick(colorDict, fig, ax)

    fig.canvas.mpl_connect("motion_notify_event", HoverClick.hover)
    fig.canvas.mpl_connect("button_press_event", HoverClick.click)

    if plotShow:
        plt.show()
    else:
        plt.savefig(f'volcano_{sample}.svg')
