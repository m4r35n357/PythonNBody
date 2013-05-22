import sys
from nbody3d import *
'''
def twoBody ():
	g = 0.05
	ts = 0.01
	errorLimit = -60.0;
	simulationTime = 1.0e3
	bodies = []
	bodies.append(Particle(1.0, 2.0, 0.0, 0.1, 0.1, 0.0, 5.0))
	bodies.append(Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0))
	variant = 0
	integratorOrder = 6
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, variant, integratorOrder)

def threeBody ():
	g = 1.0
	ts = 0.01
	errorLimit = -60.0;
	simulationTime = 1.0e3
	bodies = []
	bodies.append(Particle(1.07590, 0.0, 0.0, 0.0, 0.19509, 0.0, 1.0))
	bodies.append(Particle(-0.07095, 0.0, 0.0, -0.2, -1.23187, 0.0, 1.0))
	bodies.append(Particle(-1.00496, 0.0, 0.0, 0.0, 1.03678, 0.0, 1.0))
	variant = 0
	integratorOrder = 10
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, variant, integratorOrder)

def fourBody ():
	g = 3.5
	ts = 0.01
	errorLimit = -60.0;
	simulationTime = 1.0e3
	bodies = []
	bodies.append(Particle(1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0))
	bodies.append(Particle(-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0))
	bodies.append(Particle(1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0))
	bodies.append(Particle(-1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0))
	variant = 0
	integratorOrder = 6
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, variant, integratorOrder)
'''
def planets ():
	g = 2.95912208286e-4
	ts = 10.0
	errorLimit = -30.0;
	simulationTime = 1.0e5
	bodies = []
	bodies.append(Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0))  # Sun
	mass = 0.000954786104043
	bodies.append(Particle(-3.5025653, -3.8169847, -1.5507963, 0.00565429 * mass, -0.00412490 * mass, -0.00190589 * mass, mass))  # Jupiter
	mass = 0.000285583733151
	bodies.append(Particle(9.0755314, -3.0458353, -1.6483708, 0.00168318 * mass, 0.00483525 * mass, 0.00192462 * mass, mass))  # Saturn
	mass = 0.0000437273164546
	bodies.append(Particle(8.3101420, -16.2901086, -7.2521278, 0.00354178 * mass, 0.00137102 * mass, 0.00055029 * mass, mass))  # Uranus
	mass = 0.0000517759138449
	bodies.append(Particle(11.4707666, -25.7294829, -10.8169456, 0.00288930 * mass, 0.00114527 * mass, 0.00039677 * mass, mass))  # Neptune
	variant = 0
	integratorOrder = 1
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, variant, integratorOrder)

def eightBody ():
	g = 0.05
	ts = 0.01
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
	variant = 0
	integratorOrder = 6
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, variant, integratorOrder)

if __name__ == "__main__":
	print >> sys.stderr, __name__ + " is not an executable"
else:
	print >> sys.stderr, __name__ + " module loaded"
