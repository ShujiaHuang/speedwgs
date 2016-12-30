#!/usr/bin/env python

import sys
import math
import matplotlib.pyplot as plt
import numpy as np

sys.stdin.readline() # skip header

reducerCount = 0
# reducerID inputCount totalTime rmdupTime indexTime realignTargetCreateTime indelRealignTime baseRecTime printReadsTime UniGenotyTime
reducerStats = {}
for line in sys.stdin:
    reducerCount += 1
    vals = map(int, line.split(' '))
    reducerStats[vals[0]] = vals[1:]

# inputCountArray, totalTimeArray, rmdupArray, indexArray, realignArray, indelArray, baseArray, printArray, uniArray
infoArrays = map(np.array, zip(*[value for (key, value) in sorted(reducerStats.items())]))

# Create chr range
idxData = [
"1:249250621",
"2:243199373",
"3:198022430",
"4:191154276",
"5:180915260",
"6:171115067",
"7:159138663",
"8:146364022",
"9:141213431",
"10:135534747",
"11:135006516",
"12:133851895",
"13:115169878",
"14:107349540",
"15:102531392",
"16:90354753",
"17:81195210",
"18:78077248",
"19:59128983",
"20:63025520",
"21:48129895",
"22:51304566",
"X:155270560",
"Y:59373566",
"MT:16569"
]

chrSeq = []
chrLength = {}
chrRange = {}
sumLength = 0
for idxLine in idxData:
    tmp = idxLine.split(':') # chr:length
    chrSeq.append(tmp[0])
    chrLength[tmp[0]] = long(tmp[1])
    chrRange[tmp[0]] = [sumLength + 1, sumLength + long(tmp[1])]
    sumLength += long(tmp[1])

eachPart = int(math.ceil(1.0 * sumLength / reducerCount))
for chr in chrSeq:
    chrRange[chr][0] = chrRange[chr][0] / eachPart
    chrRange[chr][1] = chrRange[chr][1] / eachPart

## Plot
xArray = np.arange(0, reducerCount)

fig = plt.figure(figsize=(20,8))

ax2 = fig.add_subplot(111)
accArray = np.array([0] * reducerCount)
# inputCountArray, totalTimeArray, rmdupArray, indexArray, realignArray, indelArray, baseArray, printArray, uniArray
width = 1
colors = ['b','g','r','c','m','y','k']
colorIdx = 0
bars = []
for i in range(2, len(infoArrays)):
    bars.append(ax2.bar(xArray, infoArrays[i], width, color=colors[colorIdx], bottom=accArray, linewidth=0)[0])
    colorIdx = (colorIdx + 1) % len(colors)
    accArray += infoArrays[i]

ax2.plot(xArray, infoArrays[1], 'g+')

labels = ['rmdupTime','indexTime','realignTargetCreateTime','indelRealignTime','baseRecTime','printReadsTime','UniGenotyTime']
leg = ax2.legend(bars, labels, loc='best', fancybox=True)
leg.get_frame().set_alpha(0.5)

upBound = max(infoArrays[1])
for r in chrRange.values():
    ax2.plot([r[0], r[0]], [0, upBound], 'b--')
    ax2.plot([r[1], r[1]], [0, upBound], 'b--')
ax2.set_ylabel('time (s)');
ax2.set_ylim([0, 3000])

ax1 = ax2.twinx()
ax1.plot(xArray, infoArrays[0], '+')
ylim = ax1.get_ylim()
ax1.set_ylim([-ylim[1], ylim[1]])
ax1.set_ylabel('reads count')

plt.show()
plt.close()
