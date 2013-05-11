#!/usr/bin/env python

"""
"""
import math

class Symplectic(object):

    def __init__(self, g, simulationTime, timeStep, outputInterval, bodies):
		self.particles = bodies
		self.np = len(bodies)
		self.g = g
		self.timeStep = timeStep
		self.outputInterval = outputInterval
		self.iterations = simulationTime / timeStep

    def __str__(self):
        return 'np: ' + str(self.np) + ', g: ' + str(self.g) + ', ts: ' + str(self.timeStep) + ', n: ' + str(self.iterations) + ', n: ' + str(self.particles)

class Particle(object):

   def __init__(self, qX, qY, qZ, pX, pY, pZ, mass):
		self.qX = qX
		self.qY = qY
		self.qZ = qZ
		self.pX = pX
		self.pY = pY
		self.pZ = pZ
		self.mass = mass

#	def __str__(self):
# 		return 'qX: ' + str(self.qX) + ', qY: ' + str(self.qY) + ', qZ: ' + str(self.qZ) + ', pX: ' + str(self.pX) + ', pY: ' + str(self.pY) + ' pZ: ' + str(self.pZ)

def distance (xA, yA, zA, xB, yB, zB):
	return math.sqrt(math.pow(xB - xA, 2) + math.pow(yB - yA, 2) + math.pow(zB - zA, 2))

def hamiltonian (s):
	energy = 0.0
	for i in range(s.np):
		a = s.particles[i]
		energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass
		for j in range(s.np):
			if (i > j):
				b = s.particles[j]
				energy -= s.g * a.mass * b.mass / distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
	return energy

def updateQ (s, c):
	for i in range(s.np):
		a = s.particles[i]
		tmp = c / a.mass * s.timeStep
		a.qX += a.pX * tmp
		a.qY += a.pY * tmp
		a.qZ += a.pZ * tmp

def updateP (s, c):
	for i in range(s.np):
		a = s.particles[i]
		for j in range(s.np):
			b = s.particles[j]
			if (i > j):
				tmp = - c * s.g * a.mass * b.mass / math.pow(distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3) * s.timeStep
				dPx = (b.qX - a.qX) * tmp
				dPy = (b.qY - a.qY) * tmp
				dPz = (b.qZ - a.qZ) * tmp
				a.pX -= dPx
				a.pY -= dPy
				a.pZ -= dPz
				b.pX += dPx
				b.pY += dPy
				b.pZ += dPz

def euler (s, first, second):
	first(s, 1.0)
	second(s, 1.0)

def stormerVerlet2 (s, first, second):
	first(s, 0.5)
	second(s, 1.0)
	first(s, 0.5)

def stormerVerlet4 (s, first, second):
	first(s, 0.6756035959798289)
	second(s, 1.3512071919596578)
	first(s, -0.17560359597982883)
	second(s, -1.7024143839193153)
	first(s, -0.17560359597982883)
	second(s, 1.3512071919596578)
	first(s, 0.6756035959798289)

def threeBody ():
	g = 1.0
	ts = 0.01
	outputInterval = 1000
	simulationTime = 1.0e4
	bodies = []
	bodies.append(Particle(1.07590, 0.0, 0.0, 0.0, 0.19509, 0.0, 1.0))
	bodies.append(Particle(-0.07095, 0.0, 0.0, -0.2, -1.23187, 0.0, 1.0))
	bodies.append(Particle(-1.00496, 0.0, 0.0, 0.0, 1.03678, 0.0, 1.0))
	return Symplectic(g, simulationTime, ts, outputInterval, bodies)

if __name__ == "__main__":
		n = 0
		s = threeBody()
		h0 = hamiltonian(s)
		hMin = h0
		hMax = h0
		while (n <= s.iterations):
			stormerVerlet4(s, updateQ, updateP)
			hNow = hamiltonian(s)
			dH = hNow - h0
			if (hNow < hMin):
				hMin = hNow
			elif (hNow > hMax):
				hMax = hNow
			if ((n % s.outputInterval) == 0):
				l = ["["]
				for i in range(s.np):
					p = s.particles[i]
					l.append("{\"Qx\":" + str(p.qX) + ",\"Qy\":" + str(p.qY) + ",\"Qz\":" + str(p.qZ) + ",\"Px\":" + str(p.pX) + ",\"Py\":" + str(p.pY) + ",\"Pz\":" + str(p.pZ) + "},")
				print(''.join(l) + "]")
				print("t: " + str(n * s.timeStep) + ", H:" + str(hNow) + ", H0:" + str(h0) + ", H-:" + str(hMin) + ", H+:" + str(hMax) + ", ER:" + str((10.0 * math.log10(math.fabs(dH / h0)))))
			n += 1
else:
    print __name__ + " module loaded"
