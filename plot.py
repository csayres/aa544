# -*- coding: utf-8 -*-
import numpy
import itertools
import time
import scipy.interpolate
import numpy.fft
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import glob
import os
import bisect
import cPickle as pickle

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
        self.strouhalBox = (270,280,100,200)
        self._tVector = None
        self._u = None
        self._v = None
        self._p = None
        self._x = None
        self._y = None
        self._st = None # strohl number effectively frequency
        self.uvpFile = uvpFile
        self.parseHeader(uvpFile)
        self.figSuffix = "grid (%ix%i) box aspect (%s) reynolds (%.2f)" % (self.gridsize[0], self.gridsize[1], self.bluffDim[1]/self.bluffDim[0], self.re)
        # self.uvpFile2Mat(uvpFile)
        self.lagPoints = self.lagFile2Mat(lagFile) + 1.5
        self.xPts = numpy.loadtxt(xFile) + 1.5
        self.yPts = numpy.loadtxt(yFile) + 1.5
        self.residMat = self.resid2Mat(residFile)
        assert len(self.xPts)==self.gridsize[0]
        assert len(self.yPts)==self.gridsize[1]
        #boxLength = uvpFile.split("_")[-1].split(".")[0]
        self.uProfileIndex = bisect.bisect(self.xPts, numpy.max(self.lagPoints[:,1]))
        self.vProfileIndex = len(self.yPts)/2
        # get first true index

    def saveToPickle(self):
        pickleDict = {
            "_tVector": self._tVector,
            "_u": self._u,
            "_v": self._v,
            "_p": self._p,
            "_x": self._x,
            "_y": self._y,
        }
        output = open(self.figSuffix + ".pk", "wb")
        pickle.dump(pickleDict, output)
        output.close()

    def loadFromPickle(self):
            pkFile = open("./"+self.figSuffix + ".pk", "rb")
            pickleDict = pickle.load(pkFile)
            pkFile.close()
            self._tVector = pickleDict["_tVector"]
            self._u = pickleDict["_u"]
            self._v = pickleDict["_v"]
            self._p = pickleDict["_p"]
            self._x = pickleDict["_x"]
            self._y = pickleDict["_y"]
            print 'loaded from pickle'

    def tryOrLoad(self):
        try:
            self.loadFromPickle()
        except:
            print "picklefile not found"
            self.uvpFile2Mat(self.uvpFile)
            self.saveToPickle()

    @property
    def tVector(self):
        if not self._tVector:
            self.tryOrLoad()
        return self._tVector

    @property
    def u(self):
        if not self._u:
            self.tryOrLoad()
        return self._u

    @property
    def v(self):
        if not self._v:
            self.tryOrLoad()
        return self._v

    @property
    def p(self):
        if not self._p:
            self.tryOrLoad()
        return self._p

    @property
    def x(self):
        if not self._x:
            self.tryOrLoad()
        return self._x

    @property
    def y(self):
        if not self._y:
            self.tryOrLoad()
        return self._y

    def getFFTSumAndFreq(self):
        xMin = self.strouhalBox[0]
        xMax = self.strouhalBox[1]
        yMin = self.strouhalBox[2]
        yMax = self.strouhalBox[3]
        xIndMin = numpy.argmin(numpy.abs(xMin-self.xPts))
        xIndMax = numpy.argmin(numpy.abs(xMax-self.xPts))
        xInds = range(xIndMin, xIndMax)
        yIndMax = numpy.argmin(numpy.abs(yMax-self.yPts))
        yIndMin = numpy.argmin(numpy.abs(yMin-self.yPts))
        yInds = range(yIndMin,yIndMax)

        # FFT stuff
        N = len(self.tVector)
        dt = self.tVector[1]
        #freq = numpy.linspace(0.0, N/2, N/2)
        freq = numpy.linspace(0.0, 1.0/(2.0*dt), N/2)
        # ft = numpy.zeros(freq.shape, dtype="complex")
        ft = numpy.zeros(freq.shape, dtype=float)
        mid_ii = len(xInds)*len(yInds)//2
        ii = 0
        for i in xInds:
            for j in yInds:

                flatInd = j*len(self.xPts)+i
                timeSeries = numpy.asarray([x[flatInd] for x in self.u])
                if ii == mid_ii:
                    ts = timeSeries
                # _ft = numpy.fft.fft(timeSeries)[0:N/2]
                # ft = ft + _ft*numpy.conj(_ft)
                ft = ft + numpy.abs(numpy.fft.fft(timeSeries)[0:N/2])
                ii += 1
        # set strohl number
        trimmedFt = ft[2:]
        trimmedFreq = freq[2:]
        self._st = trimmedFreq[numpy.argmax(trimmedFt)]*self.bluffDim[0]
        return freq, ft, ts

    @property
    def st(self):
        if not self._st:
            self.getFFTSumAndFreq()
        return self._st

    # def getTimeseries(self, u_v_or_p, xLoc, yLoc):
    #     # find closest index to xLoc, yLoc
    #     xInd = numpy.argmin(numpy.abs(xLoc-self.xPts))
    #     yInd = numpy.argmin(numpy.abs(yLoc-self.yPts))
    #     # find index in _p array corresponding to xInd, yInd
    #     flatInd = yInd*len(self.xPts)+xInd
    #     return numpy.asarray([x[flatInd] for x in getattr(self, u_v_or_p)])

    # def getFFTAndFreq(self, dt, timeSeries):
    #     N = len(timeSeries)
    #     ft = (numpy.fft.fft(timeSeries)[0:N/2])**2
    #     freq = numpy.linspace(0.0, 1.0/(2.0*dt), N/2)
    #     # freq = numpy.fft.fftfreq(len(timeSeries), d=dt)
    #     return freq, ft

    def plotFFTSum(self):
        freq, ft, ts = self.getFFTSumAndFreq()
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(211)
        tv = numpy.asarray(self.tVector)
        ax.set_xlabel("time")
        ax.set_ylabel("pressure")
        plt.plot(tv, ts, "k")

        ax = fig.add_subplot(212)
        plt.plot(freq[1:], ft[1:], "k")
        puthere=numpy.argmax(ft[1:])
        ax.text(freq[puthere] + .002, 0.95*numpy.max(ft[1:]), 'St = %.2f'%self.st)
        ax.set_xlabel("frequency")
        ax.set_ylabel("power")

        plt.savefig(self.figSuffix + ' fftsum.png', format='png')
        plt.close()

    # def plotTimeseries(self, u_v_or_p, xLoc, yLoc, figNum):
    #     fig = plt.figure()
    #     ax = fig.add_subplot(211)
    #     x = self.getTimeseries(u_v_or_p, xLoc, yLoc)
    #     plt.plot(self._tVector, x, ".k")

    #     ax = fig.add_subplot(212)
    #     # print self._tVector[0], self._tVector[1]
    #     dt = self.tVector[1]
    #     freq, ft = self.getFFTAndFreq(dt, x)
    #     plt.plot(freq, ft, ".k")

    #     plt.savefig('timeseries_%i.png'%figNum, format='png')
    #     plt.close()

    def sumFFT(self, xLoc, yRange):
        pass


    def getTimeFromLine(self, line):
        """@param[in] line from uvpFile containing a new time point
            @return float, the time specified in the line

        line looks like this: ZONE T="t =     0.1596000000000005E+00" F=POINT, I=  200 J=   80
        """
        # keep the strin`g between the ""
        return float(line.split('"')[1].split()[-1])

    def parseHeader(self, uvpFile):
        """From a line determine the correct grid size
        """
        with open(uvpFile, 'r') as f:
            line = f.readline()
        splitted = line.split()
        self.gridsize = [int(splitted[2]), int(splitted[3])]
        self.re = float(splitted[6])
        self.bluffDim = [float(splitted[8]), float(splitted[10])]
        self.aspectRatio = self.bluffDim[1] / self.bluffDim[0]

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
        self._tVector = tVector
        self._u = uOut
        self._v = vOut
        self._p = pOut
        self._x = xArray
        self._y = yArray
        print 'loaded from uvp file'

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

    # def makeColorMaps(self, fast=True, tVector=None, saveDir=""):
    #     if not tVector:
    #         tVector = self.tVector
    #     for i in range(len(tVector)):
    #         for j in ["u"]:#, "v", "p"]:
    #             plotStr = j + " velocity" + " frame(%i) "%i + self.figSuffix
    #             fig = plt.figure(figsize=(10, 5))
    #             ax = fig.add_subplot(111,aspect="equal")
    #             lims = (-.1, .1)
    #             self.plotColorMap(i, j, figTitle=plotStr, fast=fast, lims=lims)
    #             plt.savefig(saveDir+plotStr + '.png', format='png')
    #             plt.close()

    def crop(self, data):
        ds = data.shape
        xMin = 150.
        yMin = 100.
        yMax = 200.
        xIndMin = numpy.argmin(numpy.abs(xMin-self.xPts))
        yIndMax = numpy.argmin(numpy.abs(yMax-self.yPts))
        yIndMin = numpy.argmin(numpy.abs(yMin-self.yPts))


        return data[yIndMin:yIndMax, xIndMin:]

    def plotColorMap(self, timeStep, u_v_or_p, figTitle="", fast=True, lims=(None, None), zoom=False):
        """Plot a 2D color contour

        @param[in] int, timeStep to use (-1 is last time step)
        @param[in] u_v_or_p, one of 'u', 'v', 'p'
        @param[in] figName: name of the figure
        @param[in] figTitle: title for the figure
        """
        vmin=lims[0]
        vmax=lims[1]

        # grab the 2D matrix to plot
        plt.hold(True)

        # reduce the array sizes by a factor of 10 to contain
        # x = x[::100]
        # y = y[::100]
        # z = z[::100]
        # xi = numpy.linspace(min(x), max(x), 1000)
        # yi = numpy.linspace(min(y), max(y), 2000)

        cmap = plt.get_cmap('winter')
        #cmap = plt.get_cmap('hot')

        # X, Y = numpy.meshgrid(xi, yi)
        # X, Y = numpy.meshgrid(x,y)
        # Z = ml.griddata(x, y, z, xi, yi)
        if fast:
            z = self.getUVorP(u_v_or_p, timeStep)
            Z = self.reshapeZ(z).T
            if zoom:
                Z = self.crop(Z)
            plt.imshow(Z, vmin=vmin, vmax=vmax)
        else:
            X,Y,Z = self.getGridData(timeStep, u_v_or_p)
            # Z = scipy.interpolate.griddata(x,y,z,(xi,yi))
            Z = Z.T
            if zoom:
                X = self.crop(X)
                Y = self.crop(Y)
                Z = self.crop(Z)
            # X, Y = numpy.meshgrid(x, y)
            # z = numpy.sin(X)
            # note
            # plt.contourf(X, Y, Z, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
            # plt.imshow(Z.T, cmap=cmap, vmin=-.5, vmax=1.5)
            plt.pcolormesh(X, Y, Z, cmap=cmap, vmin=vmin, vmax=vmax, norm=None)#, cmap=cmap)#, norm=norm)
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
        # pressure plots require an offset
        plt.plot(self.lagPoints[:,0], self.lagPoints[:,1], ".k-", alpha=0.8)

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
        for uvp, lag, resid, xf, yf in itertools.izip(uvpFiles,lagFiles,residFiles,xFiles,yFiles):
            dm = DataMuncher(
                uvpFile=uvp,
                lagFile = lag,
                xFile = xf,
                yFile = yf,
                residFile = resid
                )
            jefeList.append(dm)
        self.jefeList = jefeList
        self.lims = {
            "u": (-.5, 1.5),
            "v": (-.5, 1.5),
            "p": (-0.05, 0.05),

        }

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

    def makeColorMaps(self, dm, j, fast=True, tVector=None, saveDir="", zoom=False, short=False):
        if not tVector:
            tVector = dm.tVector
            if short:
                tVector = tVector[:len(tVector)//3]
        for i in range(len(tVector)):
            plotStr = j + " frame(%i) "%i + dm.figSuffix
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_subplot(111,aspect="equal")
            dm.plotColorMap(i, j, figTitle=plotStr, fast=fast, lims=self.lims[j], zoom=zoom)
            plt.savefig(saveDir+plotStr + '.png', format='png')
            plt.close()

    def plotDetector(self, dm, j, xMin, xMax, yMin, yMax, fast=True, tVector=None, saveDir="", zoom=False):
            plotStr = j + " detector region" + dm.figSuffix
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_subplot(111,aspect="equal")
            dm.plotColorMap(-1, j, figTitle=plotStr, fast=fast, lims=self.lims[j], zoom=zoom)
            box = numpy.asarray([
                [xMin, yMax],
                [xMax, yMax],
                [xMax, yMin],
                [xMin, yMin],
                [xMin, yMax],
            ])
            plt.plot(box[:,0], box[:,1], 'k')
            plt.savefig(saveDir+plotStr + '.png', format='png')
            plt.close()

    def movieReynoldsSweep(self, reRange, aspectRatio, u_v_or_p):
        dms = []
        dmsPrelim  = self.filterDataMunchers(aspectRatio=aspectRatio)
        # throw out those outside of reRange
        for dm in dmsPrelim:
            if dm.re in reRange:
                dms.append(dm)
        fig = plt.figure()
        nPanels = len(dms)
        for t in range(len(dms[0].tVector)):
            for i, dm in enumerate(dms):
                ax = fig.add_subplot(nPanels, 1, i, aspect="equal")
                figTitle = "Aspect (%i) Reynolds (%.2f)" % (dm.aspectRatio, dm.re)
                dm.plotColorMap(t, u_v_or_p, figTitle=figTitle, fast=False, lims=self.lims[u_v_or_p], zoom=False)
        plt.savefig("reSweep_%i.png"%t, format="png")
        plt.close()


    def filterDataMunchers(self, gridsize=None, re=None, aspectRatio=None):
        outList = []
        for dm in self.jefeList:
            if gridsize:
                if tuple(dm.gridsize) != tuple(gridsize):
                    continue
            if re:
                if dm.re != re:
                    continue
            if aspectRatio:
                if dm.aspectRatio != aspectRatio:
                    continue
            outList.append(dm)
        return outList[0] if len(outList)==1 else outList


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

    elJefe = elJefe("/Volumes/Boof/AA/_output_lowres")
    elJefe.movieReynoldsSweep(reRange=[70, 100, 125, 150, 175, 200, 300], aspectRatio=2, u_v_or_p="u")
    #for dm in elJefe.jefeList:
    #dm = elJefe.filterDataMunchers(re=150, aspectRatio=1)
    # dm.plotTimeseries(150, 270, 1)
    # dm.plotTimeseries(155, 270, 2)
    # dm.plotTimeseries(160, 270, 3)
    # dm.plotTimeseries(165, 270, 4)
    # dm.plotTimeseries(170, 270, 5)
    # dm.plotTimeseries(175, 270, 6)
    # dm.plotFFTSum()
    # dm.plotTimeseries("p", 270, 160, 1)
    #elJefe.plotDetector(dm, "p", 270,280,150,180, fast=False, zoom=True)
   # elJefe.makeColorMaps(dm, j="p", fast=False, zoom=True, short=False)
    #elJefe.plotProfiles()
    #elJefe.dump2dirs("flow"), 1
