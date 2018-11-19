import mmap
import pysam
import sys
import os
import re

fn = sys.argv[1]
if os.path.exists(fn):
	match = 0
	mismatch = 0
	count = 0
	with open(fn) as bed:
		for riga in bed:
			lista = re.split(r'\t+', riga)
			search = (lista[3] + '\t' + lista[0])
			with open("id_results.txt", "r") as ins:
				s = mmap.mmap(ins.fileno(), 0, access=mmap.ACCESS_READ)
				if s.find(search) != -1:
					match += 1
				else:
					mismatch +=1
			count += 1
		print "match: ",match, "\n"
		print "mismatch: ", mismatch, "\n" 
		print "total: ",count, "\n"
else:
	print "file doesn't match"
