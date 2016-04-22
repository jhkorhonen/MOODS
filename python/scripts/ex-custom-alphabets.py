#!/usr/bin/env python

import MOODS.scan
import MOODS.tools
import random

# random sequence of ACGTM
seq = ''.join([random.choice(['A', 'C', 'G', 'T', 'M']) for i in range(10000)])
# a simple matrix
mat = [
	   [1,0,0,0,0], # A
	   [0,1,0,0,0], # C
	   [0,0,1,0,0], # G
	   [0,0,0,1,0], # T
	   [0,0,0,0,4]  # M
	  ]

# this specifes which symbols in the sequence are mapped to the rows of the matrices
alph = ['aA', 'cC', 'gG', 'tT', 'mM']

# we will not bother log-odds transformations etc. this time
bg = MOODS.tools.flat_bg(5)
res = MOODS.scan.scan(seq, [mat], bg, [6], 6, alph)

for x in res[0]:
    print x.pos, x.score, seq[x.pos:x.pos+5]
