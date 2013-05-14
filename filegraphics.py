#!/usr/bin/env python

from __future__ import division
from visual import *
import json
from pprint import *

def readJson (filename):
	data = []
	for line in open(filename, 'r'):
		data.append(json.loads(line))
	return data

scene.center = (0,0,0)
scene.width = 800
scene.height = 800
scene.range = (10.0, 10.0, 10.0)
colours = [ (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.7, 0.7, 0.7), (1.0, 1.0, 1.0) ]
spheres = []
data = readJson('./aaa')
for i in range(len(data[0])):
	p = data[0][i]
	ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 0.05 * math.pow(p['mass'], 1.0 / 3.0), color = colours[i])
	ball.trail = curve(color = ball.color)
	spheres.append(ball)
for i in range(len(data)):
	rate(60)
	for j in range(len(data[0])):
		p = data[i][j]
		ball = spheres[j]
		position = (p['qX'], p['qY'], p['qZ'])
		ball.pos = position
		ball.trail.append(pos = position)
