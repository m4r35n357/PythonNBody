import math

class Particle(object):

	def __init__(self, qX, qY, qZ, pX, pY, pZ, mass):
		self.qX = qX
		self.qY = qY
		self.qZ = qZ
		self.pX = pX
		self.pY = pY
		self.pZ = pZ
		self.mass = mass

class Symplectic(object):

	def __init__(self, g, simulationTime, timeStep, errorLimit, outputInterval, bodies):
		self.particles = bodies
		self.np = len(bodies)
		self.g = g
		self.timeStep = timeStep
		self.errorLimit = errorLimit
		self.outputInterval = outputInterval
		self.iterations = simulationTime / timeStep

	def distance (self, xA, yA, zA, xB, yB, zB):
		return math.sqrt(math.pow(xB - xA, 2) + math.pow(yB - yA, 2) + math.pow(zB - zA, 2))

	def hamiltonian (self):  # Energy
		energy = 0.0
		for i in range(self.np):
			a = self.particles[i]
			energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass
			for j in range(self.np):
				if (i > j):
					b = self.particles[j]
					energy -= self.g * a.mass * b.mass / self.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
		return energy

	def updateQ (self, c):  # Update Positions
		for i in range(self.np):
			a = self.particles[i]
			tmp = c / a.mass * self.timeStep
			a.qX += a.pX * tmp
			a.qY += a.pY * tmp
			a.qZ += a.pZ * tmp

	def updateP (self, c):  # Update Momenta
		for i in range(self.np):
			a = self.particles[i]
			for j in range(self.np):
				b = self.particles[j]
				if (i > j):
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

	def stormerVerlet2 (self, first, second):  # Second order
		first(0.5)
		second(1.0)
		first(0.5)

	def stormerVerlet4 (self, first, second):  # Fourth order
		first(0.6756035959798289)
		second(1.3512071919596578)
		first(-0.17560359597982883)
		second(-1.7024143839193153)
		first(-0.17560359597982883)
		second(1.3512071919596578)
		first(0.6756035959798289)

if __name__ == "__main__":
	pass
else:
	print __name__ + " module loaded"
