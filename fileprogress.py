#!/usr/bin/env python

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from json import loads

def main():
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply a data file name <<<')
	dataFile = open(argv[1], 'r')
	line = dataFile.readline()  # throw away the first line
	line = dataFile.readline()
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.set_xlabel('Time')
	ax1.set_ylabel('Hamiltonian (Energy)', color='r')
	ax2 = ax1.twinx()
	ax2.set_ylabel('Error ratio, dBH0', color='b')
	while line:
		p = loads(line)
		ax1.plot(p['t'], p['H'], 'r.')
		ax2.plot(p['t'], p['ER'], 'b.')
		line = dataFile.readline()
	plt.show()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
