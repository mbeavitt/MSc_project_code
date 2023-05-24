#!/exports/applications/apps/SL7/anaconda/5.3.1/bin/python3

import sys

with open(f'{sys.argv[1]}/files.txt', 'r') as file:
	content = file.read().strip()

content = content.split('\n')[3:-2]

paths = []

for i in content:
	paths.append(i.split()[3])

with open(f'{sys.argv[1]}/paths.txt', 'w') as outfile:
	for i in paths:
		print(i, file = outfile)
