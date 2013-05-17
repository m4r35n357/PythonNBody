#!/opt/pypy-2.0-src/pypy/goal/pypy-c

import sys
from nbody3d import *

def twoBody ():
	g = 0.05
	ts = 0.001
	errorLimit = -60.0;
	simulationTime = 1.0e3
	bodies = []
	bodies.append(Particle(1.0, 2.0, 0.0, 0.1, 0.1, 0.0, 5.0))
	bodies.append(Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0))
	integratorOrder = 4
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, integratorOrder)

def threeBody ():
	g = 1.0
	ts = 0.001
	errorLimit = -60.0;
	simulationTime = 1.0e3
	bodies = []
	bodies.append(Particle(1.07590, 0.0, 0.0, 0.0, 0.19509, 0.0, 1.0))
	bodies.append(Particle(-0.07095, 0.0, 0.0, -0.2, -1.23187, 0.0, 1.0))
	bodies.append(Particle(-1.00496, 0.0, 0.0, 0.0, 1.03678, 0.0, 1.0))
	integratorOrder = 6
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, integratorOrder)

def fourBody ():
	g = 3.5
	ts = 0.001
	errorLimit = -60.0;
	simulationTime = 1.0e3
	bodies = []
	bodies.append(Particle(1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0))
	bodies.append(Particle(-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0))
	bodies.append(Particle(1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0))
	bodies.append(Particle(-1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0))
	integratorOrder = 6
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, integratorOrder)

def eightBody ():
	g = 0.05
	ts = 0.001
	errorLimit = -60.0;
	simulationTime = 1.0e3
	bodies = []
	bodies.append(Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0))
	bodies.append(Particle(0.0, 4.5, 0.4, -0.2, 0.0, 1.8, 2.0))
	bodies.append(Particle(-6.0, 0.0, -0.4, 0.0, -0.6, 1.0, 3.0))
	bodies.append(Particle(3.0, 0.0, -0.2, 0.0, 5.8, -0.2, 5.0))
	bodies.append(Particle(0.0, -4.0, 0.1, -3.6, 0.0, 0.2, 4.0))
	bodies.append(Particle(-4.0, 0.0, -0.1, 0.0, -0.2, -2.6, 3.0))
	bodies.append(Particle(8.0, 0.0, -0.3, 0.0, 1.2, -0.2, 3.0))
	bodies.append(Particle(0.0, 4.0, -0.2, -4.8, 0.0, -0.2, 4.0))
	integratorOrder = 6
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, integratorOrder)

def stupidPythonMain ():  # need to be inside a function to return . . .
	n = 0
	scenario = threeBody()  # create a symplectic integrator object
	h0 = scenario.hamiltonian()
	hMin = h0
	hMax = h0
	while (n <= scenario.iterations):
		scenario.solveQP()  # perform one integration step
		hNow = scenario.hamiltonian()
		tmp = math.fabs(hNow - h0)  # protect log against negative arguments
		dH = tmp if tmp > 0.0 else 1.0e-18  # protect log against small arguments
		if (hNow < hMin):
			hMin = hNow
		elif (hNow > hMax):
			hMax = hNow
		if ((n % scenario.outputInterval) == 0):
			print scenario.particlesJson()
			dbValue = 10.0 * math.log10(math.fabs(dH / h0))
			print >> sys.stderr, "t: " + str(n * scenario.timeStep) + ", H:" + str(hNow) + ", H0:" + str(h0) + ", H-:" + str(hMin) + ", H+:" + str(hMax) + ", ER:" + str(dbValue) + " dBh0"
			if (dbValue > scenario.errorLimit):
				print >> sys.stderr, "Hamiltonian error, giving up!" 
				return
		n += 1

if __name__ == "__main__":
	stupidPythonMain()
else:
	print >> sys.stderr, __name__ + " module loaded"
