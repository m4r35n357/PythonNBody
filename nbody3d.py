import sys
import math
import json

class Particle(object):

	def __init__(self, qX, qY, qZ, pX, pY, pZ, mass):
		self.qX = qX
		self.qY = qY
		self.qZ = qZ
		self.pX = pX
		self.pY = pY
		self.pZ = pZ
		self.mass = mass

	def __str__(self):
#		return "{\"qX\":" + str(self.qX) + ",\"qY\":" + str(self.qY) + ",\"qZ\":" + str(self.qZ) + ",\"pX\":" + str(self.pX) + ",\"pY\":" + str(self.pY) + ",\"pZ\":" + str(self.pZ) + ",\"mass\":" + str(self.mass) + "}"
		return "{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\": %.3e}" % (self.qX, self.qY, self.qZ, self.pX, self.pY, self.pZ, self.mass)

class Symplectic(object):

	def __init__(self, g, simulationTime, timeStep, errorLimit, bodies, order):
		self.particles = bodies
		self.np = len(bodies)
		self.pRange = range(self.np)
		self.g = g
		self.timeStep = timeStep
		self.errorLimit = errorLimit
		self.outputInterval = math.floor(0.01 / math.fabs(timeStep))
		self.iterations = simulationTime / math.fabs(timeStep)
		if (order == 1):  # First order
			self.integrator = self.euler
		elif (order == 2):  # Second order
			self.integrator = self.stormerVerlet2
		elif (order == 4):  # Fourth order
			self.integrator = self.stormerVerlet4
		elif (order == 6):  # Sixth order
			self.integrator = self.stormerVerlet6
		elif (order == 8):  # Eighth order
			self.integrator = self.stormerVerlet8
		else:  # Wrong value for integrator order
			raise Exception('>>> ERROR! Integrator order must be 1, 2, 4 or 6 <<<')

	def distance (self, xA, yA, zA, xB, yB, zB):  # Euclidean distance between point A and point B
		return math.sqrt(math.pow(xB - xA, 2) + math.pow(yB - yA, 2) + math.pow(zB - zA, 2))

	def hamiltonian (self):  # Energy
		energy = 0.0
		for i in self.pRange:
			a = self.particles[i]
			energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass
			for j in self.pRange:
				if (i > j):
					b = self.particles[j]
					energy -= self.g * a.mass * b.mass / self.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
		return energy

	def updateQ (self, c):  # Update Positions
		for i in self.pRange:
			a = self.particles[i]
			tmp = c / a.mass * self.timeStep
			a.qX += a.pX * tmp
			a.qY += a.pY * tmp
			a.qZ += a.pZ * tmp

	def updateP (self, c):  # Update Momenta
		for i in self.pRange:
			a = self.particles[i]
			for j in self.pRange:
				if (i > j):
					b = self.particles[j]
					tmp = - c * self.g * a.mass * b.mass / math.pow(self.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3) * self.timeStep
					dPx = (b.qX - a.qX) * tmp
					dPy = (b.qY - a.qY) * tmp
					dPz = (b.qZ - a.qZ) * tmp
					a.pX -= dPx
					a.pY -= dPy
					a.pZ -= dPz
					b.pX += dPx
					b.pY += dPy
					b.pZ += dPz

	def euler (self, first, second):  # First order
		first(1.0)
		second(1.0)

	def stormerVerletBase (self, first, second, step):
		first(0.5 * step)
		second(step)
		first(0.5 * step)

	def stormerVerlet2 (self, first, second):  # Second order
		self.stormerVerletBase(first, second, 1.0)

	def stormerVerlet4 (self, first, second):  # Fourth order
		self.stormerVerletBase(first, second, 1.351207191959657)
		self.stormerVerletBase(first, second, -1.702414383919315)
		self.stormerVerletBase(first, second, 1.351207191959657)

	def stormerVerlet6 (self, first, second):  # Sixth order
		self.stormerVerletBase(first, second, 0.784513610477560e0)
		self.stormerVerletBase(first, second, 0.235573213359357e0)
		self.stormerVerletBase(first, second, -1.17767998417887e0)
		self.stormerVerletBase(first, second, 1.31518632068391e0)
		self.stormerVerletBase(first, second, -1.17767998417887e0)
		self.stormerVerletBase(first, second, 0.235573213359357e0)
		self.stormerVerletBase(first, second, 0.784513610477560e0)

	def stormerVerlet8 (self, first, second):  # Eighth order
		self.stormerVerletBase(first, second, 0.104242620869991e1)
		self.stormerVerletBase(first, second, 0.182020630970714e1)
		self.stormerVerletBase(first, second, 0.157739928123617e0)
		self.stormerVerletBase(first, second, 0.244002732616735e1)
		self.stormerVerletBase(first, second, -0.716989419708120e-2)
		self.stormerVerletBase(first, second, -0.244699182370524e1)
		self.stormerVerletBase(first, second, -0.161582374150097e1)
		self.stormerVerletBase(first, second, -0.17808286265894516e1)
		self.stormerVerletBase(first, second, -0.161582374150097e1)
		self.stormerVerletBase(first, second, -0.244699182370524e1)
		self.stormerVerletBase(first, second, -0.716989419708120e-2)
		self.stormerVerletBase(first, second, 0.244002732616735e1)
		self.stormerVerletBase(first, second, 0.157739928123617e0)
		self.stormerVerletBase(first, second, 0.182020630970714e1)
		self.stormerVerletBase(first, second, 0.104242620869991e1)

	def solveQP (self):  # Update positions first
		self.integrator(self.updateQ, self.updateP)

	def solvePQ (self):  # Update momenta first
		self.integrator(self.updateP, self.updateQ)

	def particlesJson (self):
		data = []
		for i in self.pRange:
			data.append(str(self.particles[i]))
		return "[" + ','.join(data) + "]"

if __name__ == "__main__":
	pass
else:
	print >> sys.stderr, __name__ + " module loaded"
