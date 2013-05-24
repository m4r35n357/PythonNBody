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
		self.cubeRt2 = math.pow(2.0, 1.0 / 3.0)
		self.bodies = bodies
		self.np = len(bodies)
		self.pRange = range(self.np)
		self.g = g
		self.ts = timeStep
		self.eMax = errorLimit
		self.n = simulationTime / math.fabs(timeStep)  # we can run backwards too!
		if (order == 1):  # First order
			self.iterate = self.euler
		elif (order == 2):  # Second order
			self.iterate = self.stormerVerlet2
		elif (order == 4):  # Fourth order
			self.iterate = self.stormerVerlet4
		elif (order == 6):  # Sixth order
			self.iterate = self.stormerVerlet6
		elif (order == 8):  # Eighth order
			self.iterate = self.stormerVerlet8
		elif (order == 10):  # Tenth order
			self.iterate = self.stormerVerlet10
		else:  # Wrong value for integrator order
			raise Exception('>>> ERROR! Integrator order must be 1, 2, 4, 6, 8 or 10 <<<')

	def modR (self, xA, yA, zA, xB, yB, zB):  # Euclidean distance between point A and point B
		return math.sqrt(math.pow(xB - xA, 2) + math.pow(yB - yA, 2) + math.pow(zB - zA, 2))

	def h (self):  # Energy
		energy = 0.0
		for i in self.pRange:
			a = self.bodies[i]
			energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass
			for j in self.pRange:
				if (i > j):
					b = self.bodies[j]
					energy -= self.g * a.mass * b.mass / self.modR(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
		return energy

	def updateQ (self, c):  # Update Positions
		for i in self.pRange:
			a = self.bodies[i]
			tmp = c / a.mass * self.ts
			a.qX += a.pX * tmp
			a.qY += a.pY * tmp
			a.qZ += a.pZ * tmp

	def updateP (self, c):  # Update Momenta
		for i in self.pRange:
			a = self.bodies[i]
			for j in self.pRange:
				if (i > j):
					b = self.bodies[j]
					tmp = - c * self.g * a.mass * b.mass / math.pow(self.modR(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3) * self.ts
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

	def sympBase (self, c):  # build higher order integrators by composition
		self.updateQ(0.5 * c)
		self.updateP(c)
		self.updateQ(0.5 * c)

	def stormerVerlet2 (self):  # Second order
		self.sympBase(1.0)

	def stormerVerlet4 (self):  # Fourth order
		y = 1.0 / (2.0 - self.cubeRt2);
		self.sympBase(y)
		self.sympBase(- self.cubeRt2 * y)
		self.sympBase(y)

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

	def stormerVerlet10 (self):  # Tenth order
		self.sympBase(0.09040619368607278492161150)
		self.sympBase(0.53591815953030120213784983)
		self.sympBase(0.35123257547493978187517736)
		self.sympBase(-0.31116802097815835426086544)
		self.sympBase(-0.52556314194263510431065549)
		self.sympBase(0.14447909410225247647345695)
		self.sympBase(0.02983588609748235818064083)
		self.sympBase(0.17786179923739805133592238)
		self.sympBase(0.09826906939341637652532377)
		self.sympBase(0.46179986210411860873242126)
		self.sympBase(-0.33377845599881851314531820)
		self.sympBase(0.07095684836524793621031152)
		self.sympBase(0.23666960070126868771909819)
		self.sympBase(-0.49725977950660985445028388)
		self.sympBase(-0.30399616617237257346546356)
		self.sympBase(0.05246957188100069574521612)
		self.sympBase(0.44373380805019087955111365)
		self.sympBase(0.05246957188100069574521612)
		self.sympBase(-0.30399616617237257346546356)
		self.sympBase(-0.49725977950660985445028388)
		self.sympBase(0.23666960070126868771909819)
		self.sympBase(0.07095684836524793621031152)
		self.sympBase(-0.33377845599881851314531820)
		self.sympBase(0.46179986210411860873242126)
		self.sympBase(0.09826906939341637652532377)
		self.sympBase(0.17786179923739805133592238)
		self.sympBase(0.02983588609748235818064083)
		self.sympBase(0.14447909410225247647345695)
		self.sympBase(-0.52556314194263510431065549)
		self.sympBase(-0.31116802097815835426086544)
		self.sympBase(0.35123257547493978187517736)
		self.sympBase(0.53591815953030120213784983)
		self.sympBase(0.09040619368607278492161150)

	def bodiesJson (self):
		data = []
		for i in self.pRange:
			data.append(str(self.bodies[i]))
		return "[" + ','.join(data) + "]"

def icJson (fileName) :
	ic = json.loads(open(fileName, 'r').read())
	bodies = []
	for p in ic['bodies']:
		bodies.append(Particle(p['qX'], p['qY'], p['qZ'], p['pX'], p['pY'], p['pZ'], p['mass']))
	return Symplectic(ic['g'], ic['simulationTime'], ic['timeStep'], ic['errorLimit'], bodies, ic['integratorOrder'])

if __name__ == "__main__":
	print >> sys.stderr, __name__ + " is not an executable"
else:
	print >> sys.stderr, __name__ + " module loaded"
