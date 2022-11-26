#!/usr/bin/python

import glob

files = glob.glob('*.txt')

all_words = []
for filename in files:
	words = open(filename,'r').readlines()
	all_words.extend(words)

all_words = [w.strip() for w in all_words]
all_words = list(set(all_words))
all_words.sort()

output_file = open('all_stop_words.txt','w')

for w in all_words:
	output_file.write(w + '\n')

