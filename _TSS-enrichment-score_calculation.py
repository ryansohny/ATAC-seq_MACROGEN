## written by Min-Hwan Ryan Sohn on June 5th
import sys
from subprocess import call
bedtools = sys.argv[1]
tss1kb = sys.argv[2]
tssflank = sys.argv[3]
bamtobed = sys.argv[4]

call(bedtools + ' coverage -counts -sorted -a ' + tss1kb + ' -b ' + bamtobed + ' > tmp1.bed', shell=True)
call(bedtools + ' coverage -counts -sorted -a ' + tssflank + ' -b ' + bamtobed + ' > tmp2.bed', shell=True)

db = {}

dfh1 = open('tmp1.bed', 'r')
for i in dfh1:
	line = i.strip().split('\t')
	if line[-1] != '0':
		id = line[-2]
		db[id] = line[-1] + '/'
dfh1.close()
dfh2 = open('tmp2.bed', 'r')
for i in dfh2:
	line = i.strip().split('\t')
	id = line[-2]
	try:
		if 'p' not in id:
			db[id] += line[-1] + '/'
		else:
			db[id.rstrip('p')] += line[-1]
	except KeyError:
		pass
dfh2.close()
signal = [] # TSS signal
noise = [] # TSS noise
for i in db.itervalues():
	line = i.split('/')
	signal.append(int(line[0]))
	noise.append(int(line[1]) + int(line[2]))
rfh = open('tmp3.bed', 'w')
rfh.write(str(float(sum(signal)) / sum(noise)))
rfh.close()

call('/bin/rm -rf tmp1.bed', shell=True)
call('/bin/rm -rf tmp2.bed', shell=True)


