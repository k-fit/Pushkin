import dateutil
from collections import Counter

f = open('purch_data1.txt')

idlist = [] 

for line in f:
    items = line.split(',')
    idlist.append(items[0])

idcounter = Counter(idlist)

multiple = [c for c in idcounter if idcounter[c] > 2]
print len(multiple)
