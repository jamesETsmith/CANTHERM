from numpy import *
import pdb
import readGeomFc
import math
#import Gnuplot
#import Gnuplot.funcutils


class Harmonics:
  numFit = 0
  Kcos = []
  Ksin = []
  A = 0.0

  def __init__(self, numFit, Kcos, Ksin):
    self.numFit = numFit
    self.Kcos = Kcos
    self.Ksin = Ksin

  def getPotential(self, angle):
    # assuming the angle is in degrees
    rad = angle * math.pi / 180
    pot = self.A
    for i in range(self.numFit):
      pot = pot + \
          self.Ksin[i] * math.sin((i + 1) * rad) + \
          self.Kcos[i] * math.cos((i + 1) * rad)
    return pot

  def fitPotential(self, file):
    read = open(file, 'r')
    lines = read.readlines()
    angles = []
    potentials = []
    pot0 = 0.0
    # first three lines are comments and last time is repeat
    nfit = len(lines) - 4
    potgiven = []
    for i in range(3, len(lines) - 1):
      tokens = lines[i].split()

      potentials.append(float(tokens[1]))
      if i == 3:
        pot0 = potentials[0]
      potentials[i - 3] = (potentials[i - 3] - pot0) * 627.5095
      potgiven.append([float(i - 3) * 360 / nfit, potentials[i - 3]])

    # now fit the potentials
    Y = transpose(matrix(potentials[:nfit]))
    X = matrix(zeros((nfit, 11), dtype=float))
    for i in range(nfit):
      angle = float(i) * 2 * math.pi / nfit
      X[i, 0] = 1.0
      for j in range(5):
        X[i, j + 1] = math.cos((j + 1) * angle)
        X[i, 6 + j] = math.sin((j + 1) * angle)

    XtX = transpose(X) * X
    XtY = transpose(X) * Y
    b = linalg.inv(XtX) * XtY

    for i in range(5):
      self.Kcos.append(0.0)
      self.Kcos[i] = float(b[i + 1])
      self.Ksin.append(0.0)
      self.Ksin[i] = float(b[i + 6])
    self.A = float(b[0])
    self.numFit = 5

    # print self.Kcos
    # print self.Ksin

    # print "Potential-Read Potential-Fit"
    pot = []
    for i in range(3 * nfit + 1):
      angle = i * 360 / 3 / nfit
      pot.append([angle, self.getPotential(angle)])
    #   print '%14.2f'%potentials[i]+'%14.3f'%pot[i]
    # print
    # g=Gnuplot.Gnuplot()
    #g('set data style linespoints')
    #g('set xrange [0:360]')
    # g.plot(potgiven,pot)
    #raw_input('Please press enter to continue ...\n')

    return
