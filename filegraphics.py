#!/usr/bin/env python

from __future__ import division
from sys import argv
from visual import scene, sphere, curve, rate
from json import loads

def main():
	if len(argv) > 1:
		dataFile = open(argv[1], 'r')
	else:
		raise Exception('>>> ERROR! Please supply a data file name <<<')
	# scene basics
	scene.center = (0,0,0)
	scene.width = scene.height = 1000
	scene.range = (10.0, 10.0, 10.0)
	# get data dimensions
	line = dataFile.readline()
	bodies = loads(line)
	pRange = range(len(bodies))
	#  set up the balls
	colours = [ (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.5, 0.5, 0.5), (1.0, 1.0, 1.0) ]
	spheres = []
	for j in pRange:
		p = bodies[j]
		ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 0.1 * p['mass']**(1.0 / 3.0), color = colours[j])
		ball.trail = curve(color = ball.color)
		spheres.append(ball)
	while line:
		rate(60)
		bodies = loads(line)
		X = Y = Z = mT = 0.0;
		for j in pRange:  # COG correction
			p = bodies[j]
			X += p['qX'] * p['mass']
			Y += p['qY'] * p['mass']
			Z += p['qZ'] * p['mass']
			mT += p['mass']
		for j in pRange:
			p = bodies[j]
			ball = spheres[j]
			position = (p['qX'] - X / mT, p['qY'] - Y / mT, p['qZ'] - Z / mT)
			ball.pos = position
			ball.trail.append(pos = position)
		line = dataFile.readline()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
