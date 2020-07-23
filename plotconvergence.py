#!/bin/python

import math

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import csv
import time

# Set some default readable font sizes
SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rc('lines', linewidth=2) # Increase the linewidth
plt.rc('text', usetex=True)

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
    ax.set_ylabel('\#Converged eigenvalues')
    ax.set_xlim(min(iters)-1, max(iters)+1)
    ax.set_ylim(min(converged)-1, max(converged)+1)

def plotEigenvalues():
    # Implicitly assumes that each eigenvalue is present at each iteration
    stoppedLines = []
    currentLines = []
    convergedLines = []
    convergencePoints = []
    maxIter = 0
    with open('test-eigenvalue-progress.csv', newline='') as evFile, \
         open('test-residual-progress.csv', newline='') as resFile, \
         open('test-orthogonality-progress.csv', newline='') as orthoFile, \
         open('test-convergence-number-progress.csv', newline='') as convFile:
        evReader = csv.reader(evFile, delimiter=',')
        resReader = csv.reader(resFile, delimiter=',')
        orthoReader = csv.reader(orthoFile, delimiter=',')
        convReader = csv.reader(convFile, delimiter=',')
        lastConverged = 0
        for (evRow,resRow,orthoRow,convRow) in zip(evReader,resReader,orthoReader,convReader):
            # Handle the current iteration count
            curIter = float(evRow.pop(0))
            resRow.pop(0)
            orthoRow.pop(0)
            # Number of converged lines
            converged = int(convRow[1])

            maxIter = max(maxIter,curIter)

            # Current eigenvalue points
            evs = []
            for i,evres in enumerate(zip(evRow,resRow,orthoRow)):
                # Append the tuple (iter, eigenValue, order, residual, unconverged order)
                evs.append((curIter, float(evres[0]), i, float(evres[1]), \
                        max(-1,i-converged),float(evres[2])))

            # Handle that we have converged lines
            for i in range(0,lastConverged):
                ev = evs.pop(0)
            # Newly converged lines
            if converged is not lastConverged:
                for i in range(lastConverged,converged):
                    ev = evs.pop(0)
                    # Find the line that is the closest match in eigenvalue
                    line = min(currentLines, key=lambda l: abs(l[-1][1] - ev[1]))
                    currentLines.remove(line)
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

        pred = lambda p1: True

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
                if not pred(p1):
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

    def outputSegments(filename, xind, yind, cind, **kwargs):

        pred = lambda p0: True
        color = lambda c: c
        ymax = math.inf

        for key,value in kwargs.items():
            if key == "pred":
                pred = value
            elif key == "color":
                color = value
            elif key == "ymax":
                ymax = value
            else:
                raise Exception("Unknown key %" % key)

        print(ymax)
        with  open(filename, 'w') as f:
            f.write('x\ty\tcolor\n')
            skipping = False
            for line in lines:
                for i in range(0, len(line)):
                    p0 = line[i]
                    skip = not pred(p0)
                    # Skip middle point if all three above max
                    # Note this could be improved to handle the first/last points and points
                    # near jumps, but this is probably enough
                    skip = skip or (i > 0 and i < len(line)-1 and \
                            p0[yind] > ymax and line[i-1][yind] > ymax and line[i+1][yind] > ymax)
                    if skip:
                        if not skipping:
                            f.write('\n') # Create discontinuity line
                        skipping = True
                        continue
                    else:
                        skipping = False
                    f.write(f'{p0[xind]}\t{p0[yind]}\t{color(p0[cind])}\n')
                f.write('\n') # Separator line
    def outputConvergencePoints(filename, xind, yind):
        with open(filename,'w') as f:
            f.write('x\ty\n')
            for p in convergencePoints:
                f.write(f'{p[xind]}\t{p[yind]}\n')
    def ouptutConvergenceLines(filename, xind, yind, xtarget, cind, **kwargs):
        pred = lambda p0,p1: True
        color = lambda c: c

        for key,value in kwargs.items():
            if key == "pred":
                pred = value
            elif key == "color":
                color = value
            else:
                raise Exception("Unknown key %" % key)

        with  open(filename, 'w') as f:
            f.write('x\ty\tcolor\n')
            for cp in convergencePoints:
                f.write(f'{cp[xind]}\t{cp[yind]}\t{color(cp[cind])}\n')
                f.write(f'{xtarget }\t{cp[yind]}\t{color(cp[cind])}\n')
                f.write('\n') # Separator line


    # Compute the straight lines after convergece
    def computeConvergedLines(xind,yind,xtarget,cind):
        lines = []
        colors = []
        for cp in convergencePoints:
            lines.append([[cp[xind],cp[yind]], [xtarget,cp[yind]]])
            colors.append(cp[cind])
        return \
            { "lines" : lines
            , "colors" : colors
            }

    def plot1():

        plotData = computeSegments(0, 1, 4)
        norm = plt.Normalize(-9,1)
        norm = plt.Normalize(-1.5,8.5)
        lc = LineCollection(plotData['segments'], cmap='tab10', norm=norm)
        lc.set_array(np.array(plotData['segmentColors']))

        convergedData = computeConvergedLines(0,1,maxIter,4)
        lcc = LineCollection(convergedData['lines'], cmap='tab10', linestyles='dashed')
        lcc.set_array(np.array(convergedData['colors']))

        fig,ax = plt.subplots()
        

        ax.add_collection(lcc)
        line = ax.add_collection(lc)
        cbar = fig.colorbar(line, ax=ax, ticks=[-1,0,2,4,6,8])
        cbar.ax.set_yticklabels(['Conv.', 0, 2, 4, 6, 8])
        cbar.ax.set_ylabel('Eigenvalue ordering')

        ax.plot(plotData['cps']["x"], plotData['cps']["y"], 'ko')

        ax.set_xlabel('Iteration')
        ax.set_ylabel('Eigenvalue')
        ax.set_xlim(0, maxIter)
        ax.set_ylim(0.1, 1000.0)
        ax.set_yscale('log')
    plot1()

    def plot3():

        colorKey = 5
        plotData = computeSegments(0, 1, colorKey)
        plotData['segmentColors'] \
                    = list(map(lambda c: math.log(c,10), plotData['segmentColors']))
        norm = plt.Normalize(-9,1)
        lc = LineCollection(plotData['segments'], cmap='jet', norm=norm)
        lc.set_array(np.array(plotData['segmentColors']))
        
        convergedData = computeConvergedLines(0,1,maxIter,colorKey)
        convergedData['colors'] = list(map(lambda c:math.log(c,10), convergedData['colors']))
        lcc = LineCollection(convergedData['lines'], cmap='jet', \
                linestyles='dotted', norm=norm)
        lcc.set_array(np.array(convergedData['colors']))

        fig,ax = plt.subplots()
        

        ax.add_collection(lcc)
        line = ax.add_collection(lc)
        cbar = fig.colorbar(line, ax=ax)
        cbar.ax.set_ylabel('$\log(||Bx||_2)$')

        ax.plot(plotData['cps']["x"], plotData['cps']["y"], 'ko')

        ax.set_xlabel('Iteration')
        ax.set_ylabel('Eigenvalue')
        ax.set_xlim(0, maxIter)
        ax.set_ylim(0.1, 1000.0)
        ax.set_yscale('log')
        # Outputing it to file
        outputSegments('ev-res-trace.dat', 0, 1, 3, color=lambda c: round(math.log(c,10),2),ymax=10000)
        outputConvergencePoints('ev-res-trace-cp.dat', 0, 1)
        ouptutConvergenceLines('ev-res-trace-cl.dat', 0, 1, maxIter, 3, \
                color = lambda c:math.log(c,10))
    plot3()
   
    # Plot with residuals colored according to order
    def plot2():
        plotData = computeSegments(0, 3, 4, pred=lambda p1: p1[4] <= 10)
        
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
        ax.set_xlim(0, maxIter)
        ax.set_ylim(1e-8, 10)
        ax.set_yscale('log')
    plot2()


plotConverged()
plotEigenvalues()
plt.show()


