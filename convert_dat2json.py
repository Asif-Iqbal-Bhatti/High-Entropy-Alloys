#!/usr/bin/env python
#
# USAGE: Script to parse a Two column .dat file to json file
#
import json
import numpy as np
from tqdm import tqdm
import pandas as pd


with open('VDWID.dat', 'r') as f:
	data = f.read().split('\n')

dd = {}
for line in tqdm(data):
	fout = ' '.join(line.split())
	s = fout.split(" ")
	dd[s[0]] = s[1:]

json.dump(dd, open("VDWID.json", 'w'), indent = 1, sort_keys=False, ensure_ascii=True)

print(len(dd))
