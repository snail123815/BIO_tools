import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.axes_divider import AxesDivider
import bz2
import pickle


def readBedgraph(bedgraph):
    compressedDataPath = f'{os.path.splitext(bedgraph)[0]}.pickle.bz2'
    if os.path.isfile(compressedDataPath):
        with bz2.open(compressedDataPath, 'rb') as data:
            dataArray = pickle.load(data)
    else:
        print(f'Reading {bedgraph}')
        dataArray = np.zeros(8667507, dtype=[('loc', 'i4'), ('value', 'f4')])
        # This will generate a structured array
        dataArray['loc'] = np.arange(1, 8667508)
        data = np.genfromtxt(
            bedgraph, dtype=['i4', 'i4', 'f4'],
            delimiter='\t', usecols=[1, 2, 3], skip_header=1,
        )
        for start, end, num in data:
            dataArray[start - 1:end]['value'] = num

        print(f'Writing to file {compressedDataPath}')
        with bz2.open(compressedDataPath, 'wb') as data:
            pickle.dump(dataArray, data)

    return dataArray


def getAnnotation(annotationFile):
    dfAnnotation = pd.read_csv(annotationFile, delimiter='\t', skiprows=[0, 1])
    dfLoc = pd.DataFrame.from_dict(dict(zip(
        dfAnnotation.index, dfAnnotation.Location.str.split('.').values))).transpose().drop(1, axis=1)
    dfLoc.columns = ['start', 'end']
    dfLoc['start'] = dfLoc['start'].str.extract('(\d+)', expand=False)
    dfLoc['end'] = dfLoc['end'].str.extract('(\d+)', expand=False)
    dfLoc = dfLoc.apply(pd.to_numeric)
    return dfLoc, dfAnnotation


def getLocation(queryWord):
    positionsStr = dfAnnotation[dfAnnotation.Gene.str.contains(
        queryWord) | dfAnnotation.Synonym.str.contains(queryWord)]['Location']
    if len(positionsStr) == 0:
        print(f'"{queryWord}" does not exist in annotation file')
        return 'Not exist.', None
    elif len(positionsStr) > 1:
        print(f'Threre are {len(positionsStr)} instances for {queryWord}:')
        print(dfAnnotation.loc[positionsStr.index, ['Synonym', 'Location']].values)
        print('Please specify.')
        return 'Multiple hit', None
    else:
        posStr = positionsStr.values.tolist()[0].split('..')
        position = [int(i) for i in posStr]
        # now the index serves only as a test for hit
        index = positionsStr.index
        return position, index
# getLocation


def queryAndPlot(df, dfLoc, queryWord=None, dataLabel=''):
    index = None
    maxy = None
    isRangeQuery = False
    while index == None:
        if queryWord == None:
            print('\nPlease type in your query ("exit()" to exit)')
            print("query [maxy]")
            print('or')
            print("'range' start end [maxy]")
            queryWord = input('')
        if queryWord == 'exit()':
            exit()

        queryWordList = queryWord.split(' ')

        if len(queryWordList) == 2:  # maxy specified
            queryWord, maxy = queryWordList
            maxy = int(maxy)

        elif len(queryWordList) > 2 and queryWordList[0] == 'range':
            # range query
            isRangeQuery = True
            queryStart = int(queryWordList[1])
            queryEnd = int(queryWordList[2])
            # maxy can also be specified
            maxy = (int(queryWordList[3]) if len(queryWordList) == 4 else None)

        if not isRangeQuery:
            # use id query get position
            position, index = getLocation(queryWord)
            try:
                flanking = int(input('Flanking?(1000)\n'))
            except:
                flanking = 1000
            if index == None:
                queryWord = None  # clear previous input
                continue
            expandPosition = [position[0] - flanking, position[1] + flanking]
        else:
            try:
                flanking = int(input('Flanking?(0)\n'))
            except:
                flanking = 0
            index = 'range'
            expandPosition = [queryStart - flanking, queryEnd + flanking]

    plotSize = input('Plotting size? (inch) ("7,3")\n')
    try:
        figx, figy = plotSize.split(',')
        figx = int(figx)
        figy = int(figy)
    except:
        figx, figy = (7, 3)

    fig, ax = plt.subplots(1, 1, figsize=(figx * 9 / 7, figy * 4 / 3))

    # plot data
    ax.plot(range(expandPosition[0], expandPosition[1]),
            df[expandPosition[0]:expandPosition[1]]['value'], label=dataLabel)

    ax.set_xlim(*expandPosition)
    xlims = expandPosition
    if maxy:
        ax.set_ylim(ax.get_ylim()[0], maxy * 1.05)
    ylims = ax.get_ylim()
    length = xlims[1] - xlims[0]

    # plot genes as arrows
    # find the limits of gene index from dfLoc
    firstGeneIndex = dfLoc[dfLoc.end > expandPosition[0]].index[0]
    lastGeneIndex = dfLoc[dfLoc.start < expandPosition[1]].index[-1]
    for i in range(firstGeneIndex, lastGeneIndex + 1):
        label = dfAnnotation.loc[i, 'Gene']
        start, end = dfLoc.loc[i, ['start', 'end']]
        # determine the arrow start and end position
        if dfAnnotation.loc[i, 'Strand'] == '+':
            arrowStart, isPart1 = ((start, False) if start >=
                                   expandPosition[0] else (expandPosition[0], True))
            arrowEnd, isPart2 = ((end, False) if end <=
                                 expandPosition[1] else (expandPosition[1], True))
            arrowLength = arrowEnd - arrowStart + 1
            shape = 'right'
        elif dfAnnotation.loc[i, 'Strand'] == '-':
            arrowStart, isPart1 = ((end, False) if end <=
                                   expandPosition[1] else (expandPosition[1], True))
            arrowEnd, isPart2 = ((start, False) if start >=
                                 expandPosition[0] else (expandPosition[0], True))
            arrowLength = arrowEnd - arrowStart + 1
            shape = 'left'
        label = (label if not any((isPart1, isPart2)) else f'{label}(part)')
        height = (ylims[1] - ylims[0]) * 0.85
        height = (-height if sum(ylims) < 0 else height)
        width = (ylims[1] - ylims[0]) * 0.08
        xrange = xlims[1] - xlims[0]

        # arrow head width based on arrow length
        # xrange * 0.02 <= arrowLength * 0.8 --> xrange <= arrowLength * 40
        if xrange <= abs(arrowLength * 40):
            arrowHeadLength = xrange * 0.02
        else:
            arrowHeadLength = abs(arrowLength * 0.8)

        # plot arrows for this gene
        ax.arrow(arrowStart, height, arrowLength, 0, length_includes_head=True,
                 width=width, head_width=width * 1.5, head_length=arrowHeadLength,
                 shape=shape, linewidth=0)
        ax.text((arrowStart + arrowLength + arrowStart) / 2, height * 1.06,
                label, horizontalalignment='center', verticalalignment='bottom',
                bbox=dict(boxstyle="round",
                          fc=(1., 1, 1, 0.7),
                          ec=(1., 1, 1, 0.7)
                          ))

        # design a bar for none range query
        if isRangeQuery:
            pass  # the bar will plot later
        else:  # restriction with 1000 bp extended positions
            if abs(arrowStart - expandPosition[0]) == 1000 or abs(arrowStart - expandPosition[1]) == 1000:
                barStart = arrowStart
                if 0.8 * (xlims[1] - xlims[0]) > 300:
                    barLength = 300
                else:
                    barLength = (xlims[1] - xlims[0]) * 0.8 // 100 * 100
                barEnd = barStart - arrowLength / abs(arrowLength) * barLength
                # barEnd can be negative indecating the direction from bar start

    # design a bar for range query
    if isRangeQuery:
        barStart = xlims[0] + 0.8 * (xlims[1] - xlims[0])
        # not anymore start point of query
        if (xlims[1] - xlims[0]) * 0.8 > 300:
            barLength = 300
        else:
            barLength = (xlims[1] - xlims[0]) * 0.8 // 100 * 100
        barEnd = barStart - barLength
        # barEnd can be negative indecating the direction from bar start

    # plot bar
    barHeight = (ylims[1] - ylims[0]) * 0.8
    barHeight = (-barHeight * 0.95 if sum(ylims) < 0 else barHeight)
    barTextHeight = (barHeight * 0.9 if sum(ylims) < 0 else barHeight * 0.95)
    ax.plot([barStart, barEnd], [barHeight,
                                 barHeight], color='k', linewidth=4)
    ax.text(barStart + (barEnd - barStart) / 2, barTextHeight,
            f'{int(barLength)} bp', horizontalalignment='center', verticalalignment='top')
    ax.legend(loc=(0.05, 0.5))
    # plt.tight_layout()

    # adjust exact figure size
    figPos = AxesDivider(ax).get_position()
    h = [Size.Fixed(figx)]  # has to be fixed
    v = [Size.Fixed(figy)]

    divider = Divider(fig, figPos, h, v, aspect=False)

    ax.set_axes_locator(divider.new_locator(nx=0, ny=0))

    plt.show()
    saveFile = input(f"Save figure?\n{outputDir}\n(y/[n]/svg/png/pdf)")
    format = 'svg'
    if saveFile in ['y', 'svg', 'png', 'pdf']:
        if saveFile != 'y':
            format = saveFile
        fileName = f"{queryWord}{dataLabel}.{format}"
        outputFilePath = os.path.join(outputDir, fileName)
        if os.path.isfile(outputFilePath):
            overWrite = input(f"File {fileName} exists, over write? (y/n)")
            if overWrite != 'y':
                return
        fig.savefig(outputFilePath, format=format)
        print(f'Saved to:\n{outputFilePath}')
    else:
        print('Image discarded.')

    saveData = input(f"\nSave data?\n{outputDir}\n(y/[n])")
    if saveData == 'y':
        dataOut = pd.DataFrame(
            df[expandPosition[0]:expandPosition[1]],
        )
        dataOut.set_index('loc').to_excel(
            os.path.join(outputDir, f'{queryWord}{dataLabel}.xlsx')
        )
        print(f'Saved to:\n{os.path.join(outputDir, f"{queryWord}.xlsx")}')
    else:
        print('Data discarded.')
    return
# queryAndPlot


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot bedgraph with SCO genome')
    parser.add_argument('path')
    args = parser.parse_args()
    path = args.path
    outputDir = os.path.dirname(path)

    bedData = readBedgraph(path)

    annotationFile = "/Users/durand.dc/Documents/works/Resources/Resource_M145/M145_genes.ptt"
    dfLoc, dfAnnotation = getAnnotation(annotationFile)

    while True:
        queryAndPlot(bedData, dfLoc, dataLabel=os.path.split(path)[1])
