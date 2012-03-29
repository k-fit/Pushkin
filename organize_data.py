import dateutil
from collections import Counter

def clean_ids(name, outname):
    f = open(name)
    g = open(outname, 'w')
    ids = {}
    
    counter = 0

    for line in f:
        items = line.split(',')
        if items[0] in ids:
            items[0] = str(ids[items[0]])
        else:
            ids[items[0]] = counter
            items[0] = str(counter)
            counter += 1
            L = ','.join(items)
        g.write(L)
    f.close()
    g.close()

for i in range(1,4):
    name = "purch_data" + str(i) + ".txt"
    outname = "clean_purch" + str(i) + ".txt"
    print name, outname
    clean_ids(name,outname)
