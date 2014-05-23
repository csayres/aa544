# -*- coding: utf-8 -*-
import numpy
import time
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml

# from matplotlib import rc
# rc('text', usetex=True)
# rc('font', family='serif')

class DataMuncher(object):
    def __init__(self, uvpFile, lagFile, forceFile, xFile, yFile, residFile):
        """Object for interacting with fortran code output, various plotting methods, etc.

        @param[in] uvpFile: string, path to [uvpFile].dat, output from fortran routine
        @param[in] lagFile: string, path to file containing lagrangian point poisions
        @param[in] forceFile: string, path to [uvpFile].dat, output from fortran routine
        """
        self.uvpFile2Mat(uvpFile)
        # self.forceMat, foo, foo = self.uvpFile2Mat(forceFile)
        self.lagPoints = self.lagFile2Mat(lagFile)
        self.xPts = numpy.loadtxt(xFile)
        self.yPts = numpy.loadtxt(yFile)
        self.residMat = self.resid2Mat(residFile)
        assert len(self.xPts)==self.gridsize[0]
        assert len(self.yPts)==self.gridsize[1]
        boxLength = uvpFile.split("_")[-1].split(".")[0]
        self.figSuffix = "grid (%ix%i) box length (%s)" % (self.gridsize[0], self.gridsize[1], boxLength)

    def getTimeFromLine(self, line):
        """@param[in] line from uvpFile containing a new time point
            @return float, the time specified in the line

        line looks like this: ZONE T="t =     0.1596000000000005E+00" F=POINT, I=  200 J=   80
        """
        # keep the strin`g between the ""
        return float(line.split('"')[1].split()[-1])

    def getGridSizeFromLine(self, line):
        """From a line determine the correct grid size
        """
        gridsize = [int(x) for x in line.split() if x not in ("Grid", "Size")]
        assert len(gridsize) == 2
        return gridsize

    def resid2Mat(self, residualFile):
        """Convert an residual file (output from fortran code) into a 2D numpy matrix

        @param[in] residualFile: string, path to [residualFile].dat, output from fortran routine
        @return a 2D numpy array of shape n x [n, i, j, resid]
        """
        outArray = []
        with open(residualFile, 'r') as f:
            lines = f.readlines()
        for line in lines:
            # skip first two lines which only contain header info
            # split on whitespace, cast to floats
            lineArray = [float(x) for x in line.split()]
            # append to outArray
            outArray.append(lineArray)
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
            gridsize = self.getGridSizeFromLine(line1)
            line2 = f.readline() # ignored
            line3 = f.readline()
            tVector.append(self.getTimeFromLine(line3))
            uArray, vArray, pArray = [], [], []
            firstArray = True
            while True:
                line = f.readline()
                if not line:
                    # end of file
                    if len(uArray) == gridsize[0] * gridsize[1]:
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
        self.gridsize = gridsize
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

    def makeColorMaps(self, fast=False, tVector=None):
        if not tVector:
            tVector = self.tVector
        for i in len(tVector):
            for j in ["u", "v", "p"]:
                plotStr = j + " timestep(%i) " + self.figSuffix
                fig = plt.figure()
                ax = fig.add_subplot(111,aspect="equal")
                self.plotColorMap(i, j, figTitle=plotStr, fast=fast)
                plt.savefig(figName + '.eps', format='eps')
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
            plt.imshow(Z.T, cmap=cmap)
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
        if figName:
            plt.savefig(figName + '.eps', format='eps')
            plt.close()
        else:
            plt.show()

    def plotQuiver(self, timeStep, figName=None):
        plt.quiver(self.x,self.y,self.u[timeStep],self.v[timeStep])
        if figName:
            plt.savefig(figName + '.eps', format='eps')
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

if __name__ == "__main__":
    # makeFigures()
    x = DataMuncher(
        uvpFile="_output/UVP_ 64_2.dat",
        lagFile="_output/lagrangian_points 64_2.dat",
        forceFile="_output/force_grid 64_2.dat",
        xFile="_output/x_points 64_2.dat",
        yFile="_output/y_points 64_2.dat",
        residFile="_output/residual_ 64_2.dat"
    )
    # x.plotAll(-1)
    x.makeColorMaps()
    x.plotResidSemiLog("resid")

