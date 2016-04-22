#!/usr/bin/env python

import MOODS.scan
import MOODS.tools

import sys
import time
import os
import random

mat1 = [
	   [1,0,0,0,1], # A
	   [0,1,0,0,1], # C
	   [0,0,1,0,1], # G
	   [0,0,0,1,1], # T
	  ]

# AAACCCT or
# GGGGTTC
mat2 = [
	[1,1,0,0,0,0], #AA
	[0,0,1,0,0,0], #AC
	[0,0,0,0,0,0], #AG
	[0,0,0,0,0,0], #AT
	[0,0,0,0,0,0], #CA
	[0,0,0,1,1,0], #CC
	[0,0,0,0,0,0], #CG
	[0,0,0,0,0,2], #CT
	[0,0,0,0,0,0], #GA
	[0,0,0,0,0,0], #GC
	[1,1,1,0,0,0], #GG
	[0,0,0,1,0,0], #GT
	[0,0,0,0,0,0], #TA
	[0,0,0,0,0,10], #TC
	[0,0,0,0,0,0], #TG
	[0,0,0,0,1,0]  #TT
]
matrices = [mat1, mat2]

bg = MOODS.tools.flat_bg(4)
thresholds = [5,7]

seq = ''.join([random.choice(['A', 'C', 'G', 'T']) for j in range(1000)])
seq = seq + 'AAACCCT'
seq = seq + ''.join([random.choice(['A', 'C', 'G', 'T']) for j in range(1000)])
seq = seq + 'GGGGTTC'
seq = seq + ''.join([random.choice(['A', 'C', 'G', 'T']) for j in range(1000)])

results = MOODS.scan.scan_dna(seq, matrices, bg, thresholds, 2)

for mat, result in zip(matrices, results):
    for r in result:
        hitseq = seq[r.pos:r.pos+len(mat[0])]
        print str(r.pos) + '|' + hitseq + '|'  + str(r.score)
