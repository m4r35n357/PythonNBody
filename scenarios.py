#!/opt/pypy-2.0-src/pypy/goal/pypy-c

from nbody3d import *

def twoBody ():
	g = 0.05
	ts = 0.001
	outputInterval = 1000
	errorLimit = -60.0;
	simulationTime = 1.0e4
	bodies = []
	bodies.append(Particle(1.0, 2.0, 0.0, 0.1, 0.1, 0.0, 5.0))
	bodies.append(Particle(2.0, 1.0, 0.0, -0.1, -0.1, 0.0, 1.0))
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies)

def threeBody ():
	g = 1.0
	ts = 0.001
	outputInterval = 1000
	errorLimit = -60.0;
	simulationTime = 1.0e4
	bodies = []
	bodies.append(Particle(1.07590, 0.0, 0.0, 0.0, 0.19509, 0.0, 1.0))
	bodies.append(Particle(-0.07095, 0.0, 0.0, -0.2, -1.23187, 0.0, 1.0))
	bodies.append(Particle(-1.00496, 0.0, 0.0, 0.0, 1.03678, 0.0, 1.0))
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies)
	
def fourBody ():
	g = 3.5
	ts = 0.001
	outputInterval = 1000
	errorLimit = -60.0;
	simulationTime = 1.0e4
	bodies = []
	bodies.append(Particle(1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0))
	bodies.append(Particle(-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0))
	bodies.append(Particle(1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0))
	bodies.append(Particle(-1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0))
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies)

def eightBody ():
	g = 0.05
	ts = 0.001
	outputInterval = 1000
	errorLimit = -60.0;
	simulationTime = 1.0e4
	bodies = []
	bodies.append(Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0))
	bodies.append(Particle(0.0, 4.5, 0.4, -0.2, 0.0, 1.8, 2.0))
	bodies.append(Particle(-6.0, 0.0, -0.4, 0.0, -0.6, 1.0, 3.0))
	bodies.append(Particle(3.0, 0.0, -0.2, 0.0, 5.8, -0.2, 5.0))
	bodies.append(Particle(0.0, -4.0, 0.1, -3.6, 0.0, 0.2, 4.0))
	bodies.append(Particle(-4.0, 0.0, -0.1, 0.0, -0.2, -2.6, 3.0))
	bodies.append(Particle(8.0, 0.0, -0.3, 0.0, 1.2, -0.2, 3.0))
	bodies.append(Particle(0.0, 4.0, -0.2, -4.8, 0.0, -0.2, 4.0))
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies)
'''
'''
def stupidPythonMain ():  # need to be inside a function to return . . .
	n = 0
	scenario = threeBody()
	h0 = scenario.hamiltonian()
	hMin = h0
	hMax = h0
	while (n <= scenario.iterations):
		scenario.stormerVerlet4(scenario.updateQ, scenario.updateP)
		scenario.cog()
		hNow = scenario.hamiltonian()
		tmp = math.fabs(hNow - h0)  # protect log against negative arguments
		dH = tmp if tmp > 0.0 else 1.0e-18  # protect log against small arguments
		if (hNow < hMin):
			hMin = hNow
		elif (hNow > hMax):
			hMax = hNow
		if ((n % scenario.outputInterval) == 0):
			l = ["["]
			for i in range(scenario.np):
				p = scenario.particles[i]
				l.append("{\"Qx\":" + str(p.qX) + ",\"Qy\":" + str(p.qY) + ",\"Qz\":" + str(p.qZ) + ",\"Px\":" + str(p.pX) + ",\"Py\":" + str(p.pY) + ",\"Pz\":" + str(p.pZ) + "},")
			print(''.join(l) + "]")
			dbValue = 10.0 * math.log10(math.fabs(dH / h0))
			print("t: " + str(n * scenario.timeStep) + ", H:" + str(hNow) + ", H0:" + str(h0) + ", H-:" + str(hMin) + ", H+:" + str(hMax) + ", ER:" + str(dbValue) + " dBh")
			if (dbValue > scenario.errorLimit):
				print("Hamiltonian error, giving up!")
				return
		n += 1
	
if __name__ == "__main__":
	stupidPythonMain()
else:
	print __name__ + " module loaded"
