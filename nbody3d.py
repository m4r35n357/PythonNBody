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
		return "{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\":%.3e}" % (self.qX, self.qY, self.qZ, self.pX, self.pY, self.pZ, self.mass)

class Symplectic(object):

	def __init__(self, g, simulationTime, timeStep, errorLimit, bodies, order):
		self.cubeRootTwo = math.pow(2.0, 1.0 / 3.0)
		self.particles = bodies
		self.np = len(bodies)
		self.pRange = range(self.np)
		self.g = g
		self.timeStep = timeStep
		self.errorLimit = errorLimit
		self.iterations = simulationTime / math.fabs(timeStep)  # we can run backwards too!
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
			raise Exception('>>> ERROR! Integrator order must be 1, 2, 4, 6 or 8 <<<')

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

	def euler (self):  # First order
		self.updateQ(1.0)
		self.updateP(1.0)

	def sympBase (self, c):
		self.updateQ(0.5 * c)
		self.updateP(c)
		self.updateQ(0.5 * c)

	def stormerVerlet2 (self):  # Second order
		self.sympBase(1.0)

	def stormerVerlet4 (self):  # Fourth order
		gamma = 1.0 / (2.0 - self.cubeRootTwo);
		self.sympBase(gamma)
		self.sympBase(- self.cubeRootTwo * gamma)
		self.sympBase(gamma)

	def stormerVerlet6 (self):  # Sixth order
		self.sympBase(0.78451361047755726381949763)
		self.sympBase(0.23557321335935813368479318)
		self.sympBase(-1.17767998417887100694641568)
		self.sympBase(1.31518632068391121888424973)
		self.sympBase(-1.17767998417887100694641568)
		self.sympBase(0.23557321335935813368479318)
		self.sympBase(0.78451361047755726381949763)

	def stormerVerlet8 (self):  # Eighth order
		self.sympBase(0.74167036435061295344822780)
		self.sympBase(-0.40910082580003159399730010)
		self.sympBase(0.19075471029623837995387626)
		self.sympBase(-0.57386247111608226665638773)
		self.sympBase(0.29906418130365592384446354)
		self.sympBase(0.33462491824529818378495798)
		self.sympBase(0.31529309239676659663205666)
		self.sympBase(-0.79688793935291635401978884)
		self.sympBase(0.31529309239676659663205666)
		self.sympBase(0.33462491824529818378495798)
		self.sympBase(0.29906418130365592384446354)
		self.sympBase(-0.57386247111608226665638773)
		self.sympBase(0.19075471029623837995387626)
		self.sympBase(-0.40910082580003159399730010)
		self.sympBase(0.74167036435061295344822780)

	def particlesToJson (self):
		data = []
		for i in self.pRange:
			data.append(str(self.particles[i]))
		return "[" + ','.join(data) + "]"

def icJson (fileName) :
	ic = json.loads(open(fileName, 'r').read())
	bodies = []
	for p in ic['bodies']:
		bodies.append(Particle(p['qX'], p['qY'], p['qZ'], p['pX'], p['pY'], p['pZ'], p['mass']))
	return Symplectic(ic['g'], ic['simulationTime'], ic['ts'], ic['errorLimit'], bodies, ic['integratorOrder'])

if __name__ == "__main__":
	print >> sys.stderr, __name__ + " is not an executable"
else:
	print >> sys.stderr, __name__ + " module loaded"
