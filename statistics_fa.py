#####################################################################################
def get_N(sortedLength_LIST, N):
    totalLength = sum(sortedLength_LIST)
    sumLength = 0
    for idx, length in enumerate(sortedLength_LIST):
        sumLength += length
        if float(sumLength) / totalLength >= (float(N)/100):
            return length, idx+1
#####################################################################################
from optparse import OptionParser
import sys
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None:
    print('Basic usage')
    print('')
    print('     python statistics_fa.py -i test.fa')
    print('')
    sys.exit()

infile = opt.INPUT
fin = open(infile)

ref_DICT = {}
seqName_LIST = []
for line in fin:
    if line.startswith('>') == True:
        seqName = line[1:].split('\t')[0].split(' ')[0].rstrip('\n')
        seqName_LIST += [seqName]
        ref_DICT[seqName] = []
    else:
        seq = line.rstrip('\n')
        ref_DICT[seqName] += [seq]
fin.close()

seqLength_LIST = []
for seqName in seqName_LIST:
    seq = ''.join(ref_DICT[seqName])
    seqLength = len(seq)
    seqLength_LIST += [seqLength]

seqN = len(seqLength_LIST)
totalLength = sum(seqLength_LIST)
minSeqLength = min(seqLength_LIST)
maxSeqLength = max(seqLength_LIST)

seqLength_sortedLIST = sorted(seqLength_LIST, reverse=True)

print('Number of sequences:' + '\t' + str(seqN))
print('Total sequence length:' + '\t' + str(totalLength))
print('Min sequence length:' + '\t' + str(minSeqLength))
print('Max sequence length:' + '\t' + str(maxSeqLength))
print('sequence N50:' + '\t' + '\t'.join(map(str, get_N(seqLength_sortedLIST, 50))))
print('sequence N90:' + '\t' + '\t'.join(map(str, get_N(seqLength_sortedLIST, 90))))

fout = open(infile + '.' + 'lengthDistribution', 'w')
for N in range(1, 101):
    N_length = get_N(seqLength_sortedLIST, N)
    fout.write(str(N) + '\t' + '\t'.join(map(str,(N_length))) + '\n')
fout.close()

import math

logCount_LIST = [0]*(int(math.log2(maxSeqLength)) + 1)
logSum_LIST   = [0]*(int(math.log2(maxSeqLength)) + 1)
for seqLength in seqLength_sortedLIST:
    logCountIDX = int(math.log2(seqLength))
    logCount_LIST[logCountIDX] += 1
    logSum_LIST[logCountIDX]   += seqLength

fout = open(infile + '.' + 'logCount', 'w')
for logCountIDX, (logCount, logSum) in enumerate(zip(logCount_LIST, logSum_LIST)):
    fout.write('\t'.join(map(str, [logCountIDX, logCount, logSum])) + '\n')
fout.close()