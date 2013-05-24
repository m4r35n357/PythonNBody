#!/usr/bin/env python

from __future__ import division
from visual import *
import json

def main():
	if len(sys.argv) > 1:
		dataFile = open(sys.argv[1], 'r')
	else:
		raise Exception('>>> ERROR! Please supply a data file name <<<')
	# scene basics
	scene.center = (0,0,0)
	scene.width = scene.height = 1000
	scene.range = (10.0, 10.0, 10.0)
	# get data dimensions
	line = dataFile.readline()
	lineData = json.loads(line)
	particleRange = range(len(lineData))
	#  set up the balls
	colours = [ (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.5, 0.5, 0.5), (1.0, 1.0, 1.0) ]
	spheres = []
	for j in particleRange:
		p = lineData[j]
		ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 0.1 * math.pow(p['mass'], 1.0 / 3.0), color = colours[j])
		ball.trail = curve(color = ball.color)
		spheres.append(ball)
	while line:
		rate(60)
		lineData = json.loads(line)
		X = 0.0;
		Y = 0.0;
		Z = 0.0;
		mT = 0.0;
		for j in particleRange:  # COG correction
			p = lineData[j]
			X += p['qX'] * p['mass']
			Y += p['qY'] * p['mass']
			Z += p['qZ'] * p['mass']
			mT += p['mass']
		for j in particleRange:
			p = lineData[j]
			ball = spheres[j]
			position = (p['qX'] - X / mT, p['qY'] - Y / mT, p['qZ'] - Z / mT)
			ball.pos = position
			ball.trail.append(pos = position)
		line = dataFile.readline()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
