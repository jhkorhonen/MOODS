import MOODS.scan
import MOODS.tools
import random

seq = ''.join([random.choice(['A', 'C', 'G', 'T', 'M']) for i in range(10000)])
print seq

mat = [
	   [1,0,0,0,0], # A
	   [0,1,0,0,0], # C
	   [0,0,1,0,0], # G
	   [0,0,0,1,0], # T
	   [0,0,0,0,4]  # M
	  ]
alph = ['aA', 'cC', 'gG', 'tT', 'mM']

bg = MOODS.tools.flat_bg(5)
res = MOODS.scan.scan(seq, [mat], bg, [7], 6, alph)

for x in res[0]:
    print x.pos, x.score, seq[x.pos:x.pos+5]
