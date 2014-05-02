import numpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
matplotlib.use('Agg')
# from matplotlib import rc
# rc('text', usetex=True)
# rc('font', family='serif')

class DataMuncher(object):
    def __init__(self, udpFile, residualFile):
        """Object for interacting with fortran code output, various plotting methods, etc.

        @param[in] udpFile: string, path to [udpFile].dat, output from fortran routine
        @param[in] residualFile: string, path
        """
        self.udpMat, self.tVector = self.udpFile2Mat(udpFile)
        self.residMat = self.resid2Mat(residualFile)

    def getTimeFromLine(self, line):
        """@param[in] line from udpFile containing a new time point
            @return float, the time specified in the line

        line looks like this: ZONE T="t =     0.1596000000000005E+00" F=POINT, I=  200 J=   80
        """
        # keep the string between the ""
        return float(line.split('"')[1].split()[-1])

    def udpFile2Mat(self, udpFile):
        """Convert an udp file (output from fortran code) into a 2D numpy matrix

        @param[in] udpFile: string, path to [udpFile].dat, output from fortran routine
        @return outArray, tVector
            out array: a 3D numpy array of shape timeSteps x nPoints x [x, y, u, v, p]
            tVector: time vector of length timeSteps
        """
        tVector = [] # will be 1D
        gridArray = [] # will be 2D
        outArray = [] # will be 3D (will hold an array of gridArrays)
        with open(udpFile, 'r') as f:
            lines = f.readlines()
        # append the starting time to tVector
        tVector.append(self.getTimeFromLine(lines[1]))
        for line in lines[2:]:
            # start from line 2 (already grabbed time)
            if "ZONE" in line:
                # new time step encountered, parse the time, start a  gridArray
                outArray.append(gridArray)
                gridArray = []
                tVector.append(self.getTimeFromLine(line))
            else:
                # split on whitespace, cast to floats
                lineArray = [float(x) for x in line.split()]
                # append to outArray
                gridArray.append(lineArray)
        # return a numpy matrix
        return numpy.asarray(outArray, dtype=float), numpy.asarray(tVector)

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
        # plt.figure()
        plt.plot(self.residMat[:,0], numpy.log(self.residMat[:,-1]))
        plt.xlabel("Step Number")
        plt.ylabel("Log Residual")
        plt.savefig(figName + '.eps', format='eps')

    def plotColorMap(self, timeStep, u_v_or_p, figName, figTitle):
        """Plot a 2D color contour

        @param[in] int, timeStep to use (-1 is last time step)
        @param[in] u_v_or_p, one of 'u', 'v', 'p'
        @param[in] figName: name of the figure
        @param[in] figTitle: title for the figure
        """
        # grab the 2D matrix to plot
        assert u_v_or_p in ["u", "v", "p"]
        # fig = plt.figure()
        uvpIndDict = {
            "u": 2,
            "v": 3,
            "p": 4
        }
        ind = uvpIndDict[u_v_or_p]
        plotMat = self.udpMat[timeStep,:,:]
        x = plotMat[:,0]
        y = plotMat[:,1]
        z = plotMat[:,ind]
        # reduce the array sizes by a factor of 10 to contain
        # x = x[::100]
        # y = y[::100]
        # z = z[::100]
        xi = numpy.linspace(min(x), max(x), 1000)
        yi = numpy.linspace(min(y), max(y), 2000)
        cmap = plt.get_cmap('winter')
        X, Y = numpy.meshgrid(xi, yi)
        Z = ml.griddata(x, y, z, xi, yi)
        # X, Y = numpy.meshgrid(x, y)
        # z = numpy.sin(X)
        # note
        plt.contourf(X, Y, Z, cmap=cmap)#, norm=norm)
        plt.colorbar()
        plt.title(figTitle)
        # im = plt.imshow(value)
        plt.savefig(figName + '.eps', format='eps')
        # plt.show()

def makeFigures():
    """Create all the figures for homework 2.
    """
    for grid in ['coarse', 'med', 'fine']:
        x = DataMuncher(udpFile = "run2/UVP_" + grid + ".dat", residualFile = "run2/residual_" + grid + ".dat")
        x.plotResidSemiLog("resid_" + grid)
        x.plotColorMap(-1, 'u', "u_" + grid, "u " + grid  + " grid")
        x.plotColorMap(-1, 'v', "v_" + grid, "v " + grid  + " grid")
        x.plotColorMap(-1, 'p', "p_" + grid, "p " + grid  + " grid")

if __name__ == "__main__":
    makeFigures()

