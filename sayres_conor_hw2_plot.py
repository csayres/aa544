# -*- coding: utf-8 -*-
import numpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml

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
            if "VARIABLES" in line:
                continue
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
        plt.figure()
        plt.plot(self.residMat[:,0], numpy.log(self.residMat[:,-1]))
        plt.xlabel("Step Number")
        plt.ylabel("Log Residual")
        plt.savefig(figName + '.eps', format='eps')
        plt.close()

    def getXYZ(self, timeStep, u_v_or_p):
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
        return x, y, z

    def plotColorMap(self, timeStep, u_v_or_p, figName, figTitle):
        """Plot a 2D color contour

        @param[in] int, timeStep to use (-1 is last time step)
        @param[in] u_v_or_p, one of 'u', 'v', 'p'
        @param[in] figName: name of the figure
        @param[in] figTitle: title for the figure
        """
        plt.figure()
        # grab the 2D matrix to plot
        x,y,z = self.getXYZ(timeStep, u_v_or_p)

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
        plt.xlabel("x location")
        plt.ylabel("y location")
        # im = plt.imshow(value)
        plt.savefig(figName + '.eps', format='eps')
        plt.close()
        # plt.show()

    def plotQuiver(self, timeStep):
        plt.figure()
        x,y,z1 = self.getXYZ(timeStep, "u")
        x,y,z2 = self.getXYZ(timeStep, "v")
        plt.quiver(x,y,z1,z2)
        plt.show()


    def plotBlasExact(self):
        """Plot an exact blas profile
        values taken directly from Hirsch.
        """
        n = str(u"0.00000E+00 0.14142E+00 0.28284E+00 0.42426E+00 0.56569E+00 0.70711E+00 0.84853E+00 0.98995E+00 0.11314E+01 0.12728E+01 0.14142E+01 0.15556E+01 0.16971E+01 0.18385E+01 0.19799E+01 0.21213E+01 0.22627E+01 0.24042E+01 0.25456E+01 0.26870E+01 0.28284E+01 0.31113E+01 0.33941E+01 0.36770E+01 0.39598E+01 0.42426E+01 0.45255E+01 0.48083E+01 0.50912E+01 0.53740E+01 0.56569E+01 0.59397E+01 0.62225E+01 0.65054E+01 0.67882E+01 0.70711E+01 0.73539E+01 0.76368E+01 0.79196E+01 0.82024E+01 0.84853E+01")
        n = numpy.asarray([float(x.strip()) for x in n.split()], dtype=float)
        df_dn = str(u"0.00000E+00 0.46960E−01 0.93910E−01 0.14081E+00 0.18761E+00 0.23423E+00 0.28058E+00 0.32653E+00 0.37196E+00 0.41672E+00 0.46063E+00 0.50354E+00 0.54525E+00 0.58559E+00 0.62439E+00 0.66147E+00 0.69670E+00 0.72993E+00 0.76106E+00 0.79000E+00 0.81669E+00 0.86330E+00 0.90107E+00 0.93060E+00 0.95288E+00 0.96905E+00 0.98037E+00 0.98797E+00 0.99289E+00 0.99594E+00 0.99777E+00 0.99882E+00 0.99940E+00 0.99970E+00 0.99986E+00 0.99994E+00 0.99997E+00 0.99999E+00 0.99999E+00 0.10000E+01 0.10000E+01")
        df_dn = numpy.asarray([float(x.strip()) for x in df_dn.split()], dtype=float)
        ndf_dn_f = str(u"0.00000E+00 0.33177E−02 0.13282E−01 0.29858E−01 0.53025E−01 0.82696E−01 0.11873E+00 0.16098E+00 0.20916E+00 0.26296E+00 0.32193E+00 0.38563E+00 0.45345E+00 0.52475E+00 0.59881E+00 0.67483E+00 0.75202E+00 0.82955E+00 0.90656E+00 0.98226E+00 0.10558E+01 0.11940E+01 0.13167E+01 0.14209E+01 0.15051E+01 0.15720E+01 0.16215E+01 0.16569E+01 0.16812E+01 0.16972E+01 0.17072E+01 0.17133E+01 0.17168E+01 0.17187E+01 0.17198E+01 0.17203E+01 0.17206E+01 0.17207E+01 0.17207E+01 0.17208E+01 0.17208E+01")
        ndf_dn_f = numpy.asarray([float(x.strip()) for x in ndf_dn_f.split()], dtype=float)
        plt.figure()
        plt.plot(range(n), n)
        plt.close()

    def plotBlas1(self, timeStep, figName):
        # plot u(y)/U from
        plt.figure()
        x, y, z = self.getXYZ(timeStep, "u")
        # choose x = 8
        inds = numpy.nonzero(x==8)[0]
        x = x[inds]
        y = y[inds]
        z = z[inds]
        # perform an upper limit at y=0.5
        inds = numpy.nonzero(y<0.5)[0]
        x = x[inds]
        y = y[inds]
        z = z[inds]
        plt.plot(y, z)
        plt.title("Blasius profile u(y) at x = 8")
        plt.ylabel("u")
        plt.xlabel("y location")
        plt.savefig(figName + '.eps', format='eps')
        plt.close()

    def plotBlas2(self, timeStep, figName):
        # plot v(y)/U from
        plt.figure()
        x, y, z = self.getXYZ(timeStep, "v")
        # choose x = 8
        inds = numpy.nonzero(x==8)[0]
        x = x[inds]
        y = y[inds]
        z = z[inds]
        inds = numpy.nonzero(y<0.5)[0]
        x = x[inds]
        y = y[inds]
        z = z[inds]
        plt.plot(y, z)
        plt.title("Blasius profile v(y) at x = 8")
        plt.ylabel("v")
        plt.xlabel("y location")
        plt.savefig(figName + '.eps', format='eps')
        plt.close()

    def plotBlas3(self, timeStep, figName):
        # plot u(y)/U from
        plt.figure()
        x, y, z = self.getXYZ(timeStep, "u")
        # choose y = 8
        inds = numpy.nonzero(y==y[2300])[0] # about 0.05
        x = x[inds]
        y = y[inds]
        z = z[inds]
        plt.plot(x, z)
        plt.title("Blasius profile u(y) at y = 0.05")
        plt.xlabel("x location")
        plt.ylabel("u")
        plt.savefig(figName + '.eps', format='eps')
        plt.close()

    def plotBlas4(self, timeStep, figName):
        # plot u(y)/U from
        plt.figure()
        x, y, z = self.getXYZ(timeStep, "v")
        # choose y = 8
        inds = numpy.nonzero(y==y[2300])[0] # about 0.05
        x = x[inds]
        y = y[inds]
        z = z[inds]
        plt.plot(x, z)
        plt.title("Blasius profile v(y) at y = 0.05")
        plt.xlabel("x location")
        plt.ylabel("v")
        plt.savefig(figName + '.eps', format='eps')
        plt.close()

    def plotSkinFriction(self, timeStep, figName):
        """calculate the skin friction coefficient (but only on the wall section)
        average across all x greater than half way through the flow at the y = 0 boundary
        """
        x, y, u = self.getXYZ(timeStep, "u")
        # keep only values where x is > than half way, ensuring we're on the wall
        inds = numpy.nonzero(x>numpy.mean(x))[0]
        x = x[inds]
        y = y[inds]
        u = u[inds]
        # next keep only inds of y_i = 1
        inds = numpy.nonzero(y==y[1])[0]
        print 'y[1]', y[1]
        x = x[inds]
        y = y[inds]
        u = u[inds]
        # take the ratio along the boundary and plot it against x
        sf = u / y
        plt.figure()
        plt.plot(x, sf)
        plt.xlabel("x position")
        plt.ylabel("Skin Coeff")
        expected = 0.664*(100)**-0.5
        computed = 2.*numpy.mean(sf)*(1.**-3.)

        print figName, "skin coeff", numpy.mean(sf),  (computed - expected)/expected
        # plt.plot([min(x), max(x)], [analySkValue, analySkValue], 'r')
        plt.savefig(figName + '.eps', format='eps')
        plt.close()
        # plt.show()





def makeFigures():
    """Create all the figures for homework 2.
    """
    for grid in ['coarse', 'med', 'fine']:
        x = DataMuncher(udpFile = "run2/UVP_" + grid + ".dat", residualFile = "run2/residual_" + grid + ".dat")
        # x.plotResidSemiLog("resid_" + grid)
        # x.plotColorMap(-1, 'u', "u_" + grid, "u (" + grid  + " grid)")
        # x.plotColorMap(-1, 'v', "v_" + grid, "v (" + grid  + " grid)")
        # x.plotColorMap(-1, 'p', "p_" + grid, "p (" + grid  + " grid)")
        # x.plotBlas1(-1, "blas1_" + grid)
        # x.plotBlas2(-1, "blas2_" + grid)
        # x.plotBlas3(-1, "blas3_" + grid)
        # x.plotBlas4(-1, "blas4_" + grid)
        x.plotSkinFriction(-1, "sf_" + grid)

if __name__ == "__main__":
    # makeFigures()
    x = DataMuncher(udpFile="_output/UVP.dat", residualFile="run2/residual_coarse.dat")
    x.plotQuiver(0)

