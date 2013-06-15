#!/opt/pypy-2.0.2-src/pypy/goal/pypy-c

from sys import argv, stderr
from math import fabs, log10, sqrt
from json import loads

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

	def __init__(self, g, runTime, timeStep, errorLimit, bodies, order):
		self.bodies = bodies
		self.np = len(bodies)
		self.pRange = range(self.np)
		self.g = g
		self.ts = timeStep
		self.eMax = errorLimit
		self.n = runTime / fabs(timeStep)  # We can run backwards too!
		self.coefficients = []
		if (order == 2):  # Second order
			self.coefficients.append(1.0)
		elif (order == 4):  # Fourth order
			cbrt2 = 2.0 ** (1.0 / 3.0)
			y = 1.0 / (2.0 - cbrt2)
			self.coefficients.append(y)
			self.coefficients.append(- y * cbrt2)
		elif (order == 6):  # Sixth order
			self.coefficients.append(0.78451361047755726381949763)
			self.coefficients.append(0.23557321335935813368479318)
			self.coefficients.append(-1.17767998417887100694641568)
			self.coefficients.append(1.31518632068391121888424973)
		elif (order == 8):  # Eighth order
			self.coefficients.append(0.74167036435061295344822780)
			self.coefficients.append(-0.40910082580003159399730010)
			self.coefficients.append(0.19075471029623837995387626)
			self.coefficients.append(-0.57386247111608226665638773)
			self.coefficients.append(0.29906418130365592384446354)
			self.coefficients.append(0.33462491824529818378495798)
			self.coefficients.append(0.31529309239676659663205666)
			self.coefficients.append(-0.79688793935291635401978884)
		elif (order == 10):  # Tenth order
			self.coefficients.append(0.09040619368607278492161150)
			self.coefficients.append(0.53591815953030120213784983)
			self.coefficients.append(0.35123257547493978187517736)
			self.coefficients.append(-0.31116802097815835426086544)
			self.coefficients.append(-0.52556314194263510431065549)
			self.coefficients.append(0.14447909410225247647345695)
			self.coefficients.append(0.02983588609748235818064083)
			self.coefficients.append(0.17786179923739805133592238)
			self.coefficients.append(0.09826906939341637652532377)
			self.coefficients.append(0.46179986210411860873242126)
			self.coefficients.append(-0.33377845599881851314531820)
			self.coefficients.append(0.07095684836524793621031152)
			self.coefficients.append(0.23666960070126868771909819)
			self.coefficients.append(-0.49725977950660985445028388)
			self.coefficients.append(-0.30399616617237257346546356)
			self.coefficients.append(0.05246957188100069574521612)
			self.coefficients.append(0.44373380805019087955111365)
		else:  # Wrong value for integrator order
			raise Exception('>>> ERROR! Integrator order must be 2, 4, 6, 8 or 10 <<<')

	def dist (self, xA, yA, zA, xB, yB, zB):  # Euclidean distance between point A and point B
		return sqrt((xB - xA)**2 + (yB - yA)**2 + (zB - zA)**2)

	def h (self):  # Conserved energy
		energy = 0.0
		for i in self.pRange:
			a = self.bodies[i]
			energy += 0.5 * (a.pX**2 + a.pY**2 + a.pZ**2) / a.mass
			for j in self.pRange:
				if (i > j):
					b = self.bodies[j]
					energy -= self.g * a.mass * b.mass / self.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
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
					tmp = - c * self.g * a.mass * b.mass * self.ts / self.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)**3
					dPx = tmp * (b.qX - a.qX)
					dPy = tmp * (b.qY - a.qY)
					dPz = tmp * (b.qZ - a.qZ)
					a.pX -= dPx
					a.pY -= dPy
					a.pZ -= dPz
					b.pX += dPx
					b.pY += dPy
					b.pZ += dPz

	def sympBase (self, c):  # Build higher order integrators by composition
		self.updateQ(c * 0.5)
		self.updateP(c)
		self.updateQ(c * 0.5)

	def SV2 (self):  # Second order
		self.sympBase(1.0)

	def SV4 (self):  # Fourth order
		self.sympBase(self.y)
		self.sympBase(- self.y * self.cbrt2)
		self.sympBase(self.y)

	def SV6 (self):  # Sixth order
		self.sympBase(0.78451361047755726381949763)
		self.sympBase(0.23557321335935813368479318)
		self.sympBase(-1.17767998417887100694641568)
		self.sympBase(1.31518632068391121888424973)
		self.sympBase(-1.17767998417887100694641568)
		self.sympBase(0.23557321335935813368479318)
		self.sympBase(0.78451361047755726381949763)

	def SV8 (self):  # Eighth order
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

	def SV10 (self):  # Tenth order
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

	def solve (self):
		tmp = len(self.coefficients) - 1
		for i in range(0, tmp):
			self.sympBase(self.coefficients[i])
		for i in range(tmp, -1, -1):
			self.sympBase(self.coefficients[i])

	def bodiesJson (self):
		data = []
		for i in self.pRange:
			data.append(str(self.bodies[i]))
		return "[" + ','.join(data) + "]"

def icJson (fileName):
	ic = loads(open(fileName, 'r').read())
	bodies = []
	for a in ic['bodies']:
		if 'pX' in a and 'pY' in a and 'pZ' in a:  # momenta specified
			bodies.append(Particle(a['qX'], a['qY'], a['qZ'], a['pX'], a['pY'], a['pZ'], a['mass']))
		elif 'vX' in a and 'vY' in a and 'vZ' in a:  # velocities specified, convert to momenta
			mass = a['mass']
			bodies.append(Particle(a['qX'], a['qY'], a['qZ'], mass * a['vX'], mass * a['vY'], mass * a['vZ'], mass))
		else:
			raise Exception('>>> ERROR! Specify either momenta or velocites consistently <<<')
	return Symplectic(ic['g'], ic['simulationTime'], ic['timeStep'], ic['errorLimit'], bodies, ic['integratorOrder'])

def main ():  # Need to be inside a function to return . . .
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply a scenario file name <<<')
	s = icJson(argv[1])  # Create a symplectic integrator object from JSON input
	h0 = hMax = hMin = s.h()  # Set up error reporting
	print s.bodiesJson()  # Log initial particle data
	print >> stderr, '{"t":%.2f, "H":%.9e, "H0":%.9e, "H-":%.9e, "H+":%.9e, "ER":%.1f}' % (0.0, h0, h0, h0, h0, -999.9)  # Log initial progress
	n = 1
	while (n <= s.n):
		s.solve()  # Perform one full integration step
		hNow = s.h()
		tmp = fabs(hNow - h0)  # Protect logarithm against negative arguments
		dH = tmp if tmp > 0.0 else 1.0e-18  # Protect logarithm against small arguments
		if (hNow < hMin):  # Low tide
			hMin = hNow
		elif (hNow > hMax):  # High tide
			hMax = hNow
		dbValue = 10.0 * log10(fabs(dH / h0))
		print s.bodiesJson()  # Log particle data
		print >> stderr, '{"t":%.2f, "H":%.9e, "H0":%.9e, "H-":%.9e, "H+":%.9e, "ER":%.1f}' % (n * s.ts, hNow, h0, hMin, hMax, dbValue)  # Log progress
		if (dbValue > s.eMax):
#			print >> stderr, "Hamiltonian error is %.1fdBh0 (limit: %.1fdBh0), giving up!" % (dbValue, s.eMax)
			return
		n += 1

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
