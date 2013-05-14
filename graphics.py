#!/usr/bin/env python

from __future__ import division
from visual import *
from math import *
from pprint import *
from nbody3d import *
from scenarios import *

def earthSun ():
	massOfEarth = 5.98e24   # kg
	massOfSatellite = 7.36e22    # kg
	radiusOfEarth = 6.38e6  # m
	radiusOfSatellite = 1.74e6   # m
	distanceEarthToSatellite = 3.84e8  # m - the moon
	#distanceEarthToSatellite = radiusOfEarth + 
	dt = 100 # seconds (try 100 for the moon, 1 to 10 for lower satellites)
	totalseconds = 0 
	speedOfSatellite = 500 # change this to how fast we think the satellite is going in m/s
	g = 6.67e-11
	ts = dt
	errorLimit = -60.0;
	simulationTime = 1.0e2
	bodies = []
	bodies.append(Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, massOfEarth))
	bodies.append(Particle(0.0, distanceEarthToSatellite, speedOfSatellite * massOfSatellite, 0.0, 0.0, 0.0, massOfSatellite))
	integratorOrder = 4
	return Symplectic(g, simulationTime, ts, errorLimit, bodies, integratorOrder)

def plotMappings ():
	pass
	
def stupidPythonMain ():  # need to be inside a function to return . . .n = 0
	scene.center = (0,0,0)
	scene.width = 800
	scene.height = 800
	scene.range = (10.0, 10.0, 10.0)
	#scene.range = (1.0e9,1.0e9,1.0e9)
	arrowScale = 0.2
	colours = [ (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.7, 0.7, 0.7), (1.0, 1.0, 1.0) ]
	n = 0
	scenario = fourBody()  # create a symplectic integrator object
	scenario.spheres = []
	for i in scenario.pRange:
		p = scenario.particles[i]
		ball = sphere(pos = (p.qX, p.qY, p.qZ), radius = 0.1 * math.pow(p.mass, 1.0 / 3.0), color = colours[i])
		ball.trail = curve(color = ball.color)
#		ball.varr = arrow(pos = ball.pos, axis = arrowScale * vector(p.pX / p.mass, p.pY / p.mass, p.pZ / p.mass), color = colours[i])
		scenario.spheres.append(ball)
	h0 = scenario.hamiltonian()
	hMin = h0
	hMax = h0
	while (n <= scenario.iterations):
		rate(60)
		scenario.solveQP()  # perform one integration step	
		hNow = scenario.hamiltonian()
		tmp = math.fabs(hNow - h0)  # protect log against negative arguments
		dH = tmp if tmp > 0.0 else 1.0e-18  # protect log against small arguments
		if (hNow < hMin):
			hMin = hNow
		elif (hNow > hMax):
			hMax = hNow
		if ((n % scenario.outputInterval) == 0):
			for i in scenario.pRange:
				p = scenario.particles[i]
				ball = scenario.spheres[i]
				position = (p.qX - scenario.cogX, p.qY - scenario.cogY, p.qZ - scenario.cogZ)
				ball.pos = position
				ball.trail.append(pos = position)
#				ball.varr = arrow(pos = ball.pos, axis = arrowScale * vector(p.pX / p.mass, p.pY / p.mass, p.pZ / p.mass), color = colours[i])
			dbValue = 10.0 * math.log10(math.fabs(dH / h0))
#			print >> sys.stderr, "t: " + str(n * scenario.timeStep) + ", H:" + str(hNow) + ", H0:" + str(h0) + ", H-:" + str(hMin) + ", H+:" + str(hMax) + ", ER:" + str(dbValue) + " dBh"
			if (dbValue > scenario.errorLimit):
				print >> sys.stderr, "Hamiltonian error, giving up!"
				return
		n += 1

stupidPythonMain()
