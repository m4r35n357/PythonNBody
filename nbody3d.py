#!/opt/pypy-2.0-src/pypy/goal/pypy-c
##!/usr/bin/env python

import math

class Symplectic(object):

	def __init__(self, g, simulationTime, timeStep, errorLimit, outputInterval, bodies):
		self.particles = bodies
		self.np = len(bodies)
		self.g = g
		self.timeStep = timeStep
		self.errorLimit = errorLimit
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

def distance (xA, yA, zA, xB, yB, zB):
	return math.sqrt(math.pow(xB - xA, 2) + math.pow(yB - yA, 2) + math.pow(zB - zA, 2))

def hamiltonian (s):  # Energy
	energy = 0.0
	for i in range(s.np):
		a = s.particles[i]
		energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass
		for j in range(s.np):
			if (i > j):
				b = s.particles[j]
				energy -= s.g * a.mass * b.mass / distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
	return energy

def updateQ (s, c):  # Update Positions
	for i in range(s.np):
		a = s.particles[i]
		tmp = c / a.mass * s.timeStep
		a.qX += a.pX * tmp
		a.qY += a.pY * tmp
		a.qZ += a.pZ * tmp

def updateP (s, c):  # Update Momenta
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

def euler (s, first, second):  # First order
	first(s, 1.0)
	second(s, 1.0)

def stormerVerlet2 (s, first, second):  # Second order
	first(s, 0.5)
	second(s, 1.0)
	first(s, 0.5)

def stormerVerlet4 (s, first, second):  # Fourth order
	first(s, 0.6756035959798289)
	second(s, 1.3512071919596578)
	first(s, -0.17560359597982883)
	second(s, -1.7024143839193153)
	first(s, -0.17560359597982883)
	second(s, 1.3512071919596578)
	first(s, 0.6756035959798289)

if __name__ == "__main__":
	pass
else:
	print __name__ + " module loaded"
