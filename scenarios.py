#!/opt/pypy-2.0-src/pypy/goal/pypy-c

import sys
from nbody3d import *

def main ():  # need to be inside a function to return . . .
	n = 0
	if len(sys.argv) > 1:
		s = icJson(sys.argv[1])  # create a symplectic integrator object from JSON input
	else:
		raise Exception('>>> ERROR! Please supply a s file name <<<')
	print s.bodiesJson()  # log initial particle data
	h0 = s.h()  # set up error reporting
	hMin = h0
	hMax = h0
	while (n <= s.n):
		s.iterate()  # perform one full integration step
		print s.bodiesJson()  # log particle data
		hNow = s.h()
		tmp = math.fabs(hNow - h0)  # protect logarithm against negative arguments
		dH = tmp if tmp > 0.0 else 1.0e-18  # protect logarithm against small arguments
		if (hNow < hMin):  # low tide
			hMin = hNow
		elif (hNow > hMax):  # high tide
			hMax = hNow
		dbValue = 10.0 * math.log10(math.fabs(dH / h0))
		print >> sys.stderr, 't:%.2f, H:%.9e, H0:%.9e, H-:%.9e, H+:%.9e, ER:%.1fdBh0' % (n * s.ts, hNow, h0, hMin, hMax, dbValue)  # log progress
		if (dbValue > s.eMax):
			print >> sys.stderr, "Hamiltonian error is %.1fdBh0 (limit: %.1fdBh0), giving up!" % (dbValue, s.eMax)
			return
		n += 1

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
