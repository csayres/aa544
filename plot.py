# -*- coding: utf-8 -*-
import numpy
import itertools
import time
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import glob
import os
import bisect

# from matplotlib import rc
# rc('text', usetex=True)
# rc('font', family='serif')

class DataMuncher(object):
    def __init__(self, uvpFile, lagFile, xFile, yFile, residFile):
        """Object for interacting with fortran code output, various plotting methods, etc.

        @param[in] uvpFile: string, path to [uvpFile].dat, output from fortran routine
        @param[in] lagFile: string, path to file containing lagrangian point poisions
        @param[in] forceFile: string, path to [uvpFile].dat, output from fortran routine
        """
        self.uvpFile2Mat(uvpFile)
        self.lagPoints = self.lagFile2Mat(lagFile)
        self.xPts = numpy.loadtxt(xFile)
        self.yPts = numpy.loadtxt(yFile)
        self.residMat = self.resid2Mat(residFile)
        assert len(self.xPts)==self.gridsize[0]
        assert len(self.yPts)==self.gridsize[1]
        #boxLength = uvpFile.split("_")[-1].split(".")[0]
        self.figSuffix = "grid (%ix%i) box aspect (%s) reynolds (%.2f)" % (self.gridsize[0], self.gridsize[1], self.bluffSize[1]/self.bluffSize[0], self.reynolds)
        self.uProfileIndex = bisect.bisect(self.xPts, numpy.max(self.lagPoints[:,1]))
        self.vProfileIndex = len(self.yPts)/2
        # get first true index

    def getTimeFromLine(self, line):
        """@param[in] line from uvpFile containing a new time point
            @return float, the time specified in the line

        line looks like this: ZONE T="t =     0.1596000000000005E+00" F=POINT, I=  200 J=   80
        """
        # keep the strin`g between the ""
        return float(line.split('"')[1].split()[-1])

    def parseHeader(self, line):
        """From a line determine the correct grid size
        """
        splitted = line.split()
        self.gridsize = [int(splitted[2]), int(splitted[3])]
        self.reynolds = float(splitted[6])
        self.bluffSize = [float(splitted[8]), float(splitted[10])]

    def resid2Mat(self, residualFile):
        """Convert an residual file (output from fortran code) into a 2D numpy matrix

        @param[in] residualFile: string, path to [residualFile].dat, output from fortran routine
        @return a 2D numpy array of shape n x [n, i, j, resid]
        """
        outArray = []
        with open(residualFile, 'r') as f:
            lines = f.readlines()
        for ind, line in enumerate(lines):
            # skip first two lines which only contain header info
            # split on whitespace, cast to floats
            resid = float(line.split()[-1])
            # append to outArray
            outArray.append([ind, resid])
        # return a numpy matrix
        return numpy.asarray(outArray, dtype=float)

    def plotResidSemiLog(self, figName):
        """Plot the residual vs number of steps
        """
        plt.figure()
        plt.plot(self.residMat[:,0], numpy.log(self.residMat[:,-1]), '.k')
        plt.xlabel("Step Number")
        plt.ylabel("Log Residual")
        plt.savefig(figName + '.eps', format='eps')
        plt.close()

    def lagFile2Mat(self, lagFile):
        """Convert a lagFile into a 2D array
        """
        return numpy.loadtxt(lagFile)

    def uvpFile2Mat(self, uvpFile):
        """Convert an uvp file (output from fortran code) into a 2D numpy matrix

        @param[in] uvpFile: string, path to [uvpFile].dat, output from fortran routine
        @return outArray, tVector
            out array: a 3D numpy array of shape timeSteps x nPoints x [x, y, u, v, p]
            tVector: time vector of length timeSteps
        """
        tVector = [] # will be 1D
        uOut = [] # will be 2D
        vOut = []
        pOut = []
        xArray = []
        yArray = []
        outArray = [] # will be 3D (will hold an array of gridArrays)
        with open(uvpFile, 'r') as f:
            line1 = f.readline()
            self.parseHeader(line1)
            line2 = f.readline() # ignored
            line3 = f.readline()
            tVector.append(self.getTimeFromLine(line3))
            uArray, vArray, pArray = [], [], []
            firstArray = True
            while True:
                line = f.readline()
                if not line:
                    # end of file
                    if len(uArray) == self.gridsize[0] * self.gridsize[1]:
                        print 'got all data!'
                        uOut.append(numpy.asarray(uArray, dtype=float))
                        vOut.append(numpy.asarray(vArray, dtype=float))
                        pOut.append(numpy.asarray(pArray, dtype=float))
                    else:
                        # remove the last time point, we dont have
                        # all the data
                        tVector.pop(-1)
                    break
                if "VARIABLES" in line:
                    continue
                if "ZONE" in line:
                    # new time step encountered, parse the time, start a  gridArray
                    uOut.append(numpy.asarray(uArray, dtype=float))
                    vOut.append(numpy.asarray(vArray, dtype=float))
                    pOut.append(numpy.asarray(pArray, dtype=float))
                    uArray, vArray, pArray = [], [], []
                    tVector.append(self.getTimeFromLine(line))
                    firstArray = False
                else:
                    # split on whitespace, cast to floats
                    lineArray = [float(x) for x in line.split()]
                    if firstArray:
                        xArray.append(lineArray[0])
                        yArray.append(lineArray[1])
                    uArray.append(lineArray[2])
                    vArray.append(lineArray[3])
                    pArray.append(lineArray[4])
        # return a numpy matrix
        self.tVector = tVector
        self.u = uOut
        self.v = vOut
        self.p = pOut
        self.x = xArray
        self.y = yArray

    def reshapeZ(self, Z):
        Z = numpy.reshape(Z, (self.gridsize[0],self.gridsize[1]), order="F")
        return Z

    def getUVorP(self, u_v_or_p, timeStep):
        assert u_v_or_p in ["u", "v", "p"]
        # fig = plt.figure()
        if u_v_or_p == "u":
            z = self.u[timeStep]
        elif u_v_or_p == "v":
            z = self.v[timeStep]
        else:
            z = self.p[timeStep]
        return z

    def getGridData(self, timeStep, u_v_or_p):
        z = self.getUVorP(u_v_or_p, timeStep)
        X,Y = numpy.meshgrid(self.xPts,self.yPts)
        Z = self.reshapeZ(z)
        return X,Y,Z

    def plotVelocityProfiles(self, timeStep):
        fig = plt.figure()
        x,y,zu = self.getGridData(timeStep, "u")
        x,y,zv = self.getGridData(timeStep, "v")
        ax = fig.add_subplot(1,2,1)
        plt.plot(self.xPts, zu[:,self.uProfileIndex])
        plt.title("u")
        ax = fig.add_subplot(1,2,2)
        plt.plot(self.yPts, zv[self.vProfileIndex,:])
        plt.title("v")
        plt.savefig('velocityprofiles.pdf', format='pdf')
        plt.close()



    def makeColorMaps(self, fast=True, tVector=None, saveDir=""):
        if not tVector:
            tVector = self.tVector
        for i in range(len(tVector)):
            for j in ["u"]:#, "v", "p"]:
                plotStr = j + " velocity" + " frame(%i) "%i + self.figSuffix
                fig = plt.figure(figsize=(10, 5))
                ax = fig.add_subplot(111,aspect="equal")
                self.plotColorMap(i, j, figTitle=plotStr, fast=fast)
                plt.savefig(saveDir+plotStr + '.png', format='png')
                plt.close()

    def plotColorMap(self, timeStep, u_v_or_p, figTitle="", fast=True):
        """Plot a 2D color contour

        @param[in] int, timeStep to use (-1 is last time step)
        @param[in] u_v_or_p, one of 'u', 'v', 'p'
        @param[in] figName: name of the figure
        @param[in] figTitle: title for the figure
        """

        # grab the 2D matrix to plot
        plt.hold(True)

        # reduce the array sizes by a factor of 10 to contain
        # x = x[::100]
        # y = y[::100]
        # z = z[::100]
        # xi = numpy.linspace(min(x), max(x), 1000)
        # yi = numpy.linspace(min(y), max(y), 2000)

        cmap = plt.get_cmap('winter')

        # X, Y = numpy.meshgrid(xi, yi)
        # X, Y = numpy.meshgrid(x,y)
        # Z = ml.griddata(x, y, z, xi, yi)
        if fast:
            z = self.getUVorP(u_v_or_p, timeStep)
            Z = self.reshapeZ(z)

            plt.imshow(Z.T, vmin=-0.5, vmax=1.5)
        else:
            X,Y,Z = self.getGridData(timeStep, u_v_or_p)
            # Z = scipy.interpolate.griddata(x,y,z,(xi,yi))

            # X, Y = numpy.meshgrid(x, y)
            # z = numpy.sin(X)
            # note
            # plt.contourf(X, Y, Z, cmap=cmap, vmin=-0.5, vmax=1.5)#, norm=norm)
            # plt.imshow(Z.T, cmap=cmap, vmin=-.5, vmax=1.5)
            plt.pcolormesh(X, Y, Z.T, vmin=-0.5, vmax=1.5)#, cmap=cmap)#, norm=norm)
            # img = plot.imshow()
            self.plotLagPoints()

        plt.colorbar()
        plt.title(figTitle)
        plt.xlabel("x location")
        plt.ylabel("y location")
        # im = plt.imshow(value)

    def plotQuiver(self, timeStep, figName=None):
        plt.quiver(self.x,self.y,self.u[timeStep],self.v[timeStep])
        if figName:
            plt.savefig(figName + '.png', format='png')
            plt.close()
        else:
            plt.show()

    def plotLagPoints(self):
        plt.plot(self.lagPoints[:,0], self.lagPoints[:,1], ".k-", alpha=0.5)


    def plotAll(self, timestep):
        fig = plt.figure();
        ax = fig.add_subplot(111, aspect='equal');
        plt.hold(True);
        x.plotQuiver(timestep);
        x.plotPressure(timestep);
        x.plotLagPoints();
        x.plotQuiverForce(timestep)
        plt.show()

def createDataMuncher(gridsize, boxLength):
    if gridsize < 100:
        gs = ' %i'%gridsize
    else:
        gs = '%i'%gridsize
    return DataMuncher(
        uvpFile="_output/UVP_%s_%i.dat"%(gs, boxLength),
        lagFile="_output/lagrangian_points%s_%i.dat"%(gs, boxLength),
        forceFile="_output/force_grid%s_%i.dat"%(gs, boxLength),
        xFile="_output/x_points%s_%i.dat"%(gs, boxLength),
        yFile="_output/y_points%s_%i.dat"%(gs, boxLength),
        residFile="_output/residual_%s_%i.dat"%(gs, boxLength),
        )

def makeDataMunchDict(gridsizes, boxLengths):
    outDict = {}
    for gridsize in strgs:
        outDict[str(gridsize)] = {}
        for boxLength in boxLengths:
            outDict[str(gridsize)][str(boxLength)] = createDataMuncher(gridsize, boxLength)
    return outDict

class elJefe(object):
    """Object for managing / plotting all runs!
    """
    def __init__(self, fileDir):
        self.fileDir = fileDir
        allFiles = glob.glob(fileDir + "/*")
        uvpFiles = []
        residFiles = []
        xFiles = []
        yFiles = []
        lagFiles = []
        for f in allFiles:
            if 'lagrangian' in f:
                lagFiles.append(f)
            elif 'UVP' in f:
                uvpFiles.append(f)
            elif 'residual' in f:
                residFiles.append(f)
            elif 'x_points' in f:
                xFiles.append(f)
            elif 'y_points' in f:
                yFiles.append(f)
        # sortem
        uvpFiles.sort()
        residFiles.sort()
        xFiles.sort()
        yFiles.sort()
        lagFiles.sort()
        jefeList = []
        ii = 0
        for uvp, lag, resid, xf, yf in itertools.izip(uvpFiles,lagFiles,residFiles,xFiles,yFiles):
            dm = DataMuncher(
                uvpFile=uvp,
                lagFile = lag,
                xFile = xf,
                yFile = yf,
                residFile = resid
                )
            dirName = "kknumber %i "%ii + dm.figSuffix
            os.mkdir(dirName)
            dm.makeColorMaps(saveDir=dirName + "/")
            ii += 1
        self.jefeList = jefeList

    def plotUVP_resArray(self):
        inds = [22,4,7]
        fig = plt.figure()
        plotnum = 1
        for ii,ind in enumerate(inds):
            dm = self.jefeList[ind]
            for jj, j in enumerate(["u", "v", "p"]):
                plotStr = j
                ax = fig.add_subplot(3, 3, plotnum, aspect="equal")
                dm.plotColorMap(-1, j, figTitle=plotStr, fast=False)
                plotnum += 1
        plt.savefig('resarray.pdf', format='pdf')
        plt.close()

    def plotProfiles(self):
        dm = self.jefeList[0]
        dm.plotVelocityProfiles(-1)

    def dump2dirs(self, base):
        for int, dm in enumerate(self.jefeList):
            dirName = base + str(int)
            os.mkdir(dirName)
            dm.makeColorMaps(saveDir=dirName+"/")





def cleanUpDir(d):
    allFiles = glob.glob(d+"/*")
    for f in allFiles:
        fnew = f[:]
        fnew = fnew.replace(" ", "")
        fnew = fnew.replace("..", ".")
        print f, fnew
        os.rename(f, fnew)


if __name__ == "__main__":
    # makeFigures()
    #x = createDataMuncher(256, 5)
    # x = DataMuncher(
    #     uvpFile="_output/UVP.dat",
    #     lagFile="_output/lagrangian_points.dat",
    #     xFile="_output/x_points.dat",
    #     yFile="_output/y_points.dat",
    #     residFile="_output/residual.dat",
    #     )
    # x.plotAll(-1)
    # x.makeColorMaps()
    # x.plotResidSemiLog("resid")

    elJefe = elJefe("_output")
    #elJefe.jefeList[0].makeColorMaps(fast=False)
    #elJefe.plotUVP_resArray()
    #elJefe.plotProfiles()
    #elJefe.dump2dirs("flow")

