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

	def __init__(self, g, simulationTime, timeStep, errorLimit, bodies, variant, order):
		self.particles = bodies
		self.np = len(bodies)
		self.pRange = range(self.np)
		self.g = g
		self.timeStep = timeStep
		self.errorLimit = errorLimit
		self.iterations = simulationTime / math.fabs(timeStep)  # we can run backwards too!
		self.variant = variant
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
#		elif (order == 10):  # Eighth order
#			self.integrator = self.stormerVerlet10
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
		self.stormerVerletBase(first, second, 0.78451361047755726381949763)
		self.stormerVerletBase(first, second, 0.23557321335935813368479318)
		self.stormerVerletBase(first, second, -1.17767998417887100694641568)
		self.stormerVerletBase(first, second, 1.31518632068391121888424973)
		self.stormerVerletBase(first, second, -1.17767998417887100694641568)
		self.stormerVerletBase(first, second, 0.23557321335935813368479318)
		self.stormerVerletBase(first, second, 0.78451361047755726381949763)

	def stormerVerlet8 (self, first, second):  # Eighth order
		self.stormerVerletBase(first, second, 0.74167036435061295344822780)
		self.stormerVerletBase(first, second, -0.40910082580003159399730010)
		self.stormerVerletBase(first, second, 0.19075471029623837995387626)
		self.stormerVerletBase(first, second, -0.57386247111608226665638773)
		self.stormerVerletBase(first, second, 0.29906418130365592384446354)
		self.stormerVerletBase(first, second, 0.33462491824529818378495798)
		self.stormerVerletBase(first, second, 0.31529309239676659663205666)
		self.stormerVerletBase(first, second, -0.79688793935291635401978884)
		self.stormerVerletBase(first, second, 0.31529309239676659663205666)
		self.stormerVerletBase(first, second, 0.33462491824529818378495798)
		self.stormerVerletBase(first, second, 0.29906418130365592384446354)
		self.stormerVerletBase(first, second, -0.57386247111608226665638773)
		self.stormerVerletBase(first, second, 0.19075471029623837995387626)
		self.stormerVerletBase(first, second, -0.40910082580003159399730010)
		self.stormerVerletBase(first, second, 0.74167036435061295344822780)
'''
	def stormerVerlet10 (self, first, second):  # Tenth order
		self.stormerVerletBase(first, second, 0.09040619368607278492161150)
		self.stormerVerletBase(first, second, 0.53591815953030120213784983)
		self.stormerVerletBase(first, second, 0.35123257547493978187517736)
		self.stormerVerletBase(first, second, -0.31116802097815835426086544)
		self.stormerVerletBase(first, second, -0.52556314194263510431065549)
		self.stormerVerletBase(first, second, 0.14447909410225247647345695)
		self.stormerVerletBase(first, second, 0.02983588609748235818064083)
		self.stormerVerletBase(first, second, 0.17786179923739805133592238)
		self.stormerVerletBase(first, second, 0.09826906939341637652532377)
		self.stormerVerletBase(first, second, 0.46179986210411860873242126)
		self.stormerVerletBase(first, second, -0.33377845599881851314531820)
		self.stormerVerletBase(first, second, 0.07095684836524793621031152)
		self.stormerVerletBase(first, second, 0.23666960070126868771909819)
		self.stormerVerletBase(first, second, -0.49725977950660985445028388)
		self.stormerVerletBase(first, second, -0.30399616617237257346546356)
		self.stormerVerletBase(first, second, 0.05246957188100069574521612)
		self.stormerVerletBase(first, second, 0.44373380805019087955111365)
		self.stormerVerletBase(first, second, 0.05246957188100069574521612)
		self.stormerVerletBase(first, second, -0.30399616617237257346546356)
		self.stormerVerletBase(first, second, -0.49725977950660985445028388)
		self.stormerVerletBase(first, second, 0.23666960070126868771909819)
		self.stormerVerletBase(first, second, 0.07095684836524793621031152)
		self.stormerVerletBase(first, second, -0.33377845599881851314531820)
		self.stormerVerletBase(first, second, 0.46179986210411860873242126)
		self.stormerVerletBase(first, second, 0.09826906939341637652532377)
		self.stormerVerletBase(first, second, 0.17786179923739805133592238)
		self.stormerVerletBase(first, second, 0.02983588609748235818064083)
		self.stormerVerletBase(first, second, 0.14447909410225247647345695)
		self.stormerVerletBase(first, second, -0.52556314194263510431065549)
		self.stormerVerletBase(first, second, -0.31116802097815835426086544)
		self.stormerVerletBase(first, second, 0.35123257547493978187517736)
		self.stormerVerletBase(first, second, 0.53591815953030120213784983)
		self.stormerVerletBase(first, second, 0.09040619368607278492161150)
'''
	def solve (self):
		if self.variant == 0:  # Update positions first
			self.integrator(self.updateQ, self.updateP)
		elif self.variant == 1:  # Update momenta first
			self.integrator(self.updateP, self.updateQ)

	def particlesToJson (self):
		data = []
		for i in self.pRange:
			data.append(str(self.particles[i]))
		return "[" + ','.join(data) + "]"

if __name__ == "__main__":
	print >> sys.stderr, __name__ + " is not an executable"
else:
	print >> sys.stderr, __name__ + " module loaded"
