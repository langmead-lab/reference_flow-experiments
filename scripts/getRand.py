#Takes in a file containing basically any string (separated by lines)
#Randomly selects N strings and outputs them to terminal

import re, sys
import random

inputFile = sys.argv[1]
randsize = int(sys.argv[2])


bigList = []
with open(inputFile) as myFile:
    for line in myFile:
        bigList.append(line.strip())

subset=random.sample(bigList, randsize)

for k in subset:
    print k
