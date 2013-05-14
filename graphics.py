#!/usr/bin/env python
#!/opt/pypy-2.0-src/pypy/goal/pypy-c

from __future__ import division
from visual import *
from math import *
from pprint import *
from nbody3d import *

# this converts totalseconds to a nice string (days, hours, minutes and seconds)
def make_time_string(t):
    "Accept a number of seconds, return a relative string."
    if t < 60: return "%02i seconds"%t
    mins,secs = divmod(t, 60)
    if mins < 60: return "%02i minutes %02i seconds"%(mins,secs)
    hours, mins = divmod(mins,60)
    if hours < 24: return "%02i hours %02i minutes"%(hours,mins)
    days,hours = divmod(hours,24)
    return "%02i days %02i hours"%(days,hours)

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
	outputInterval = 1
	errorLimit = -60.0;
	simulationTime = 1.0e2
	bodies = []
	bodies.append(Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, massOfEarth))
	bodies.append(Particle(0.0, distanceEarthToSatellite, speedOfSatellite * massOfSatellite, 0.0, 0.0, 0.0, massOfSatellite))
	integratorOrder = 4
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies, integratorOrder)

def threeBody ():
	g = 1.0
	ts = 0.01
	outputInterval = 1
	errorLimit = -60.0;
	simulationTime = 1.0e4
	bodies = []
	bodies.append(Particle(1.07590, 0.0, 0.0, 0.0, 0.19509, 0.0, 1.0))
	bodies.append(Particle(-0.07095, 0.0, 0.0, -0.2, -1.23187, 0.0, 1.0))
	bodies.append(Particle(-1.00496, 0.0, 0.0, 0.0, 1.03678, 0.0, 1.0))
	integratorOrder = 4
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies, integratorOrder)

def fourBody ():
	g = 3.51
	ts = 0.01
	outputInterval = 1
	errorLimit = -60.0;
	simulationTime = 1.0e4
	bodies = []
	bodies.append(Particle(1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0))
	bodies.append(Particle(-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0))
	bodies.append(Particle(1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0))
	bodies.append(Particle(-1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0))
	integratorOrder = 4
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies, integratorOrder)

def eightBody ():
	g = 0.05
	ts = 0.01
	outputInterval = 1
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
	integratorOrder = 4
	return Symplectic(g, simulationTime, ts, errorLimit, outputInterval, bodies, integratorOrder)
	
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
	for i in scenario.indices:
		p = scenario.particles[i]
		ball = sphere(pos = (p.qX, p.qY, p.qZ), radius = 0.05 * math.pow(p.mass, 1.0 / 3.0), color = colours[i])
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
#			l = ["["]
			for i in scenario.indices:
				p = scenario.particles[i]
				ball = scenario.spheres[i]
				position = (p.qX - scenario.cogX, p.qY - scenario.cogY, p.qZ - scenario.cogZ)
				ball.pos = position
				ball.trail.append(pos = position)
#				ball.varr = arrow(pos = ball.pos, axis = arrowScale * vector(p.pX / p.mass, p.pY / p.mass, p.pZ / p.mass), color = colours[i])
#				l.append("{\"Qx\":" + str(p.qX) + ",\"Qy\":" + str(p.qY) + ",\"Qz\":" + str(p.qZ) + ",\"Px\":" + str(p.pX) + ",\"Py\":" + str(p.pY) + ",\"Pz\":" + str(p.pZ) + "},")
#			print(''.join(l) + "]")
			dbValue = 10.0 * math.log10(math.fabs(dH / h0))
			print("t: " + str(n * scenario.timeStep) + ", H:" + str(hNow) + ", H0:" + str(h0) + ", H-:" + str(hMin) + ", H+:" + str(hMax) + ", ER:" + str(dbValue) + " dBh")
			if (dbValue > scenario.errorLimit):
				print("Hamiltonian error, giving up!")
				return
		n += 1

stupidPythonMain()
