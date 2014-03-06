import sys
from collections import defaultdict
startdict = defaultdict(list)
enddict = defaultdict(list)
readdict=defaultdict(list)

if __name__ == '__main__':
    with open (sys.argv[1],'rU') as fin:
        fout = open(sys.argv[2],'w')
        for line in fin:
            l = line.strip().split('\t')
            #print(l)
            ID = l[0]
            #print(ID)
            subID = ID [0:ID.find('/')].strip()
            start = l[2].strip()
            end = l[3].strip()
            read = l[16].strip()
            startdict[subID].append(int(start))
            enddict[subID].append(int(end))
            readdict[subID].append(read)
#print(ID,start,startdict)
#print(enddict)
    newlist=[]        
    for key in startdict:
        newlist.append([key,min(startdict[key]),max(enddict[key]),readdict[key]])
#print(newlist)
    for item in newlist:
        values = '\t'.join(str(k) for k in item)+'\n'
    
        
        fout.write(values)
    fout.close()
