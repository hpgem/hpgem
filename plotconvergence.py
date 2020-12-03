#!/bin/python

import math

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import csv
import time

def plotConverged():
    iters = []
    converged=[]
    with open('test-convergence-number-progress.csv',newline='') as csvfile:
        convReader = csv.reader(csvfile,delimiter=',')
        for row in convReader:
           (i,c) = row
           iters.append(float(i))
           converged.append(float(c))

    fig, ax = plt.subplots()

    ax.plot(iters, converged)  # Matplotlib plot.
    ax.set_xlabel('Iteration')
    ax.set_ylabel('#Converged eigenvalues')
    ax.set_xlim(min(iters)-1, max(iters)+1)
    ax.set_ylim(min(converged)-1, max(converged)+1)

def plotEigenvalues():
    # Implicitly assumes that each eigenvalue is present at each iteration
    stoppedLines = []
    currentLines = []
    convergedLines = []
    convergencePoints = []
    with open('test-eigenvalue-progress.csv', newline='') as evFile, \
         open('test-residual-progress.csv', newline='') as resFile, \
         open('test-convergence-number-progress.csv', newline='') as convFile:
        evReader = csv.reader(evFile, delimiter=',')
        resReader = csv.reader(resFile, delimiter=',')
        convReader = csv.reader(convFile, delimiter=',')
        lastConverged = 0
        for (evRow,resRow,convRow) in zip(evReader,resReader,convReader):
            # Handle the current iteration count
            curIter = float(evRow.pop(0))
            resRow.pop(0)
            # Number of converged lines
            converged = int(convRow[1])

            # Current eigenvalue points
            evs = []
            for i,evres in enumerate(zip(evRow,resRow)):
                # Append the tuple (iter, eigenValue, order, residual, unconverged order)
                evs.append((curIter, float(evres[0]), i, float(evres[1]), max(-1,i-converged)))

            # Handle that we have converged lines
            for i in range(0,lastConverged):
                ev = evs.pop(0)
                convergedLines[i].append(ev)
            # Newly converged lines
            if converged is not lastConverged:
                for i in range(lastConverged,converged):
                    # So far it seems only the EV with lowest order can converge
                    line = min(currentLines, key=lambda l: l[-1][2])
                    currentLines.remove(line)
                    ev = evs.pop(0)
                    line.append(ev)
                    convergedLines.append(line)
                    convergencePoints.append(ev)
                    print("Converged %f (%d) at iter %d" % (ev[1],i,curIter))
                lastConverged = converged
            # Sort by eigenvalues to allow matching
            evs = sorted(evs, key=lambda e: e[1])

            # Match to previous
            while len(currentLines) > 1 and \
                abs(currentLines[0][-1][1] - evs[0][1]) > \
                    abs(currentLines[1][-1][1] - evs[0][1]):
                # Second one matches better
                line = currentLines.pop(0)
                stoppedLines.append(line)
                print("Dropped line at %d" % curIter)
            # Append it all to the previous lines
            i = 0
            for ev in evs:
                if i < len(currentLines):
                    currentLines[i].append(ev)
                else:
                    currentLines.append([ev])
                    print("New line %d at %d" % (i,curIter))
                i = i+1
    # All lines have been processed, create one set of them
    lines = stoppedLines
    lines.extend(currentLines)
    lines.extend(convergedLines)

    def computeSegments(xind,yind,cind, **kwargs):
        segments = []
        segmentColors = []

        pred = lambda p0,p1: True

        for key,value in kwargs.items():
            if key == "pred":
                pred = value
            else:
                raise Exception("Unknown key %" % key)

        for line in lines:
            # Create a segment between subsequent iterations
            for i in range(0, len(line)-1):
                p0 = line[i]
                p1 = line[i+1]
                if not pred(p0,p1):
                    continue
                segment = [ (p0[xind], p0[yind]), (p1[xind], p1[yind]) ]
                segments.append(segment)
                segmentColors.append(p1[cind])

        cps = {
                "x" : list(map(lambda p: p[xind], convergencePoints)),
                "y" : list(map(lambda p: p[yind], convergencePoints))
            }
        return  \
            { "segments" : segments
            , "segmentColors" : segmentColors
            , "cps" : cps
            }

    def plot1():

        plotData = computeSegments(0, 1, 3)
        plotData['segmentColors'] \
                    = list(map(lambda c: math.log(c,10), plotData['segmentColors']))
        norm = plt.Normalize(-7,1)
        lc = LineCollection(plotData['segments'], cmap='viridis', norm=norm)
        lc.set_array(np.array(plotData['segmentColors']))

        fig,ax = plt.subplots()
        

        line = ax.add_collection(lc)
        cbar = fig.colorbar(line, ax=ax)
        cbar.ax.set_ylabel('Log(residual)')

        ax.plot(plotData['cps']["x"], plotData['cps']["y"], 'ko')

        ax.set_xlabel('Iteration')
        ax.set_ylabel('Eigenvalue')
        ax.set_xlim(0, 1500)
        ax.set_ylim(0.1, 100.0)
        ax.set_yscale('log')
    plot1()
   
    # Plot with residuals colored according to order
    def plot2():
        plotData = computeSegments(0, 3, 4, pred=lambda p0,p1: min(p0[1],p1[1]) <= 100)
        
        norm = plt.Normalize(-1.5,8.5)
        lc = LineCollection(plotData['segments'], cmap='tab10', norm=norm)
        #lc = LineCollection(segments, cmap='viridis', norm=norm)
        lc.set_array(np.array(plotData["segmentColors"]))

        fig,ax = plt.subplots()
        
        line = ax.add_collection(lc)
        cbar = fig.colorbar(line, ax=ax, ticks=[-1,0,2,4,6,8])
        cbar.ax.set_yticklabels(['Conv.', 0, 2, 4, 6, 8])
        cbar.ax.set_ylabel('Eigenvalue ordering')

        ax.plot(plotData["cps"]["x"],plotData["cps"]["y"],'ko')

        ax.set_xlabel('Iteration')
        ax.set_ylabel('Residual')
        ax.set_xlim(0, 1500)
        ax.set_ylim(1e-8, 10)
        ax.set_yscale('log')
    plot2()


plotConverged()
plotEigenvalues()
plt.show()


