#!/opt/homebrew/bin/python3

#thanks to Jan Bernauer for writing the original version of this
import sys
from numpy import array
from math import sqrt

f=open(sys.argv[1])
maxGoodDist=float(sys.argv[2])
#print ("dist threshold = ", maxGoodDist)
#scan for contour
first=True
while True:
    # find contour
    while True:
        l=f.readline()
        if not l:
            exit()
        #look for 'Contour' header.  break when we find it.
        tok=l.split()
        if len(tok)==2 and tok[0]=="Contour":
            contnum=int(tok[1])
            break
    firstline=array([float(v) for v in f.readline().split()[:3]])
    nlines=1
    sumpos=firstline
    while True:
        l=f.readline()
        if not l:
            exit()
        if l.strip()=="":
            break
        lastline=array([float(v) for v in l.split()[:3]])
        #nlines=nlines+1
        #sumpos=sumpos+lastline
    #avepos=sumpos/nlines
    #if contnum==11 or contnum==223:
    #    print(contnum,"center:",avepos)
    dist=sqrt(sum((firstline-lastline)**2))
    if dist>maxGoodDist:
        if first:
            print (sys.argv[1])
            first=False
        print (contnum,"Delta:",dist)
        
