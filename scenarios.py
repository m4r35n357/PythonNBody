#!/opt/pypy-2.0-src/pypy/goal/pypy-c

import sys
from nbody3d import *

def planets ():  # here for the time being, until I get a JSON version written
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

def icJson (fileName) :
	ic = json.loads(open(fileName, 'r').read())
	bodies = []
	for p in ic['bodies']:
		bodies.append(Particle(p['qX'], p['qY'], p['qZ'], p['pX'], p['pY'], p['pZ'], p['mass']))
	return Symplectic(ic['g'], ic['simulationTime'], ic['ts'], ic['errorLimit'], bodies, ic['variant'], ic['integratorOrder'])

def main ():  # need to be inside a function to return . . .
	n = 0
	if len(sys.argv) > 1:
		scenario = icJson(sys.argv[1])  # create a symplectic integrator object from JSON input
	else:
		scenario = planets()  # create a symplectic integrator object using a function above
	h0 = scenario.hamiltonian()
	hMin = h0
	hMax = h0
	while (n <= scenario.iterations):
		scenario.solve()  # perform one integration step
		print scenario.particlesToJson()  # particle data
		hNow = scenario.hamiltonian()
		tmp = math.fabs(hNow - h0)  # protect log against negative arguments
		dH = tmp if tmp > 0.0 else 1.0e-18  # protect log against small arguments
		if (hNow < hMin):
			hMin = hNow
		elif (hNow > hMax):
			hMax = hNow
		dbValue = 10.0 * math.log10(math.fabs(dH / h0))
		print >> sys.stderr, 't:%.2f, H:%.9e, H0:%.9e, H-:%.9e, H+:%.9e, ER:%.1fdBh0' % (n * scenario.timeStep, hNow, h0, hMin, hMax, dbValue)  # progress
		if (dbValue > scenario.errorLimit):
			print >> sys.stderr, "Hamiltonian error, giving up!" 
			return
		n += 1

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
