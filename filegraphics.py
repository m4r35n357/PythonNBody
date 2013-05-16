#!/usr/bin/env python

from __future__ import division
from visual import *
import json
from pprint import *
from nbody3d import *

def readJson (filename):
	data = []
	for line in open(filename, 'r'):
		data.append(json.loads(line))
	return data

if __name__ == "__main__":
	scene.center = (0,0,0)
	scene.width = 1000
	scene.height = 1000
	scene.range = (10.0, 10.0, 10.0)

	colours = [ (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.5, 0.5, 0.5), (1.0, 1.0, 1.0) ]
	spheres = []

	data = readJson(sys.argv[1])
	rowRange = range(len(data))
	particleRange = range(len(data[0]))
	for j in particleRange:
		p = data[0][j]
		ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 0.1 * math.pow(p['mass'], 1.0 / 3.0), color = colours[j])
		ball.trail = curve(color = ball.color)
		spheres.append(ball)
	for i in rowRange:
		rate(60)
		X = 0.0;
		Y = 0.0;
		Z = 0.0;
		mT = 0.0;
		for j in particleRange:  # COG correction
			p = data[i][j]
			X += p['qX'] * p['mass']
			Y += p['qY'] * p['mass']
			Z += p['qZ'] * p['mass']
			mT += p['mass']
		for j in particleRange:
			p = data[i][j]
			ball = spheres[j]
			position = (p['qX'] - X / mT, p['qY'] - Y / mT, p['qZ'] - Z / mT)
			ball.pos = position
			ball.trail.append(pos = position)
else:
	print >> sys.stderr, __name__ + " module loaded"
