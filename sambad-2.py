import sys
import os
getBin = lambda x, n: x >= 0 and str(bin(x))[2:].zfill(n) or "-" + str(bin(x))[3:].zfill(n)
if __name__ == '__main__':
    with open (sys.argv[1]) as fin:
	fh = open(sys.argv[2],'w')
	fh.write('ID'+'\t'+'chromsome number'+'\t'+'start'+'\t'+'end'+'\t'+'supplementary alignment'+'\t'+'read is PCR or optical duplicate'+'\t'+'read fails platform/vendor quality checks'+'\t'+'not primary alignment'+'\t'+'second in pair'+'\t'+'first in pair'+'\t'+'mate reverse strand'+'\t'+'read reverse strand'+'\t'+'mate unmapped'+'\t'+'read unmapped'+'\t'+'read mapped in proper pair'+'\t'+'read paired'+'\t'+'read sequence'+'\n')
        for line in fin:
            l = line.strip().split('\t')
            start = int(l[3])-1
            length = len(l[9])
            end = start+ length
            bedrow = str(l[2])+'\t'+ str(start)+ '\t'+str(end)+'\t'
            bi=getBin(int(l[1]),12)
            de = list(str(bi))
	    fh.write(str(l[0])+'\t'+bedrow+'\t'.join(de)+'\t'+str(l[9])+'\n')
	fh.close()


