#!/usr/bin/python

import MOODS.scan
import MOODS.tools


import sys
import time
import os
import math

import timeit
from numpy import median


def load_adm(filename):
    with open(filename, "r") as file_handle:
        lines = [[i for i in line.split()] for line in file_handle]

        mat_1storder = [[float(i) for i in line[:-2]] for line in lines[:16]]
        zero_terms = [float(line[0]) for line in lines[-4:]]

        mat_0order = [[float(i) for i in line[:-1]] for line in lines[-4:]]
        
        return mat_0order, (mat_1storder, zero_terms)


def load_chr(N):
	filename = "paper-test-data/chroms/chr{}.fa".format(N)
	with open(filename, "r") as file_handle:
		file_handle.next()
		return "".join([line.strip() for line in file_handle])


def log_odds_dna_1storder(m, bg):
	mat = m[0]
	zero_terms = m[1]
	los = [[0] * len(row) for row in mat]

	for i in xrange(len(mat[0])):
		for j in xrange(16):
			los[j][i] = math.log(mat[j][i]) - math.log(bg[j & 3])

	for k in xrange(4):
		for j in xrange(4):
			los[(j << 2) ^ k][0] = los[(j << 2) ^ k][0] + math.log(zero_terms[j]) + math.log(bg[j])
	return los


def setup_scanner_0order(threshold):
	matrix_directory = sys.argv[1]
	matrix_names = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.adm']
	raw_matrices = [load_adm(matrix_directory + filename) for filename in matrix_names]
	bg = MOODS.tools.flat_bg(4)

	log_odds_0order = [MOODS.tools.log_odds(m[0], bg, 0.1) for m in raw_matrices]
	thresholds_0order = [MOODS.tools.threshold_from_p(m, bg, threshold) for m in log_odds_0order]

	scanner = MOODS.scan.Scanner(7)
	scanner.set_motifs(log_odds_0order, bg, thresholds_0order)

	return scanner


def setup_scanner_1storder(threshold):
	matrix_directory = sys.argv[1]
	matrix_names = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.adm']
	raw_matrices = [load_adm(matrix_directory + filename) for filename in matrix_names]

	bg = MOODS.tools.flat_bg(4)
	log_odds_1storder = [log_odds_dna_1storder(m[1], bg) for m in raw_matrices]
	thresholds_1storder =  [MOODS.tools.threshold_from_p(m, bg, threshold, 4) for m in log_odds_1storder]

	scanner = MOODS.scan.Scanner(7)
	scanner.set_motifs(log_odds_1storder, bg, thresholds_1storder)
	return scanner

def scan_chr(N, scanner):
	seq = load_chr(N)
	res = scanner.scan(seq)

def chr_len(N):
	seq = load_chr(N)
	return len(seq)

print "-- preprocessing --"
for th in [0.00001,0.0001,0.001]:
	print "0-order", th, timeit.timeit("setup_scanner_0order(th)",setup="from __main__ import setup_scanner_0order, th", number=100)
	print "first-order", th, timeit.timeit("setup_scanner_1storder(th)",setup="from __main__ import setup_scanner_1storder, th", number=100)


print "-- chr lengths --"
lengths = []
for i in map(str, range(1,23)):
	lengths.append(chr_len(i))
print "lengths =", lengths

print "-- scanning --"
for th in [0.00001,0.0001,0.001]:
	print "0-order", th
	scanner = setup_scanner_0order(th)
	times = []
	for i in map(str, range(1,23)):
		ts = timeit.repeat("scan_chr(i, scanner)",setup="from __main__ import scan_chr, i, scanner", repeat=6, number=1)
		times.append(ts)

	print "times =", times
	print "median_times =", [median(ts) for ts in times]


for th in [0.00001,0.0001,0.001]:
	print "1st-order", th
	scanner = setup_scanner_1storder(th)
	times = []
	for i in map(str, range(1,23)):
		ts = timeit.repeat("scan_chr(i, scanner)",setup="from __main__ import scan_chr, i, scanner", repeat=6, number=1)
		times.append(ts)
	print "times =", times
	print "median_times =", [median(ts) for ts in times]





