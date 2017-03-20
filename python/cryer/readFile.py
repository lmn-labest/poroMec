# -*- coding: utf-8 -*-

def readCsv(name,x,y,step,cX,cY,scaleX,scaleY):
    
    with open(name,"r") as f:
        data = f.read()    
    
    lines = data.split('\n')

    lines = lines[1:len(lines)-1]

    for line in lines:
        a = line.split(',')
        x.append(float(a[cX])*scaleX)
        y.append(float(a[cY])/scaleY)
     
     
     
def readTxt(name,x,y,cX,cY,scaleX,scaleY):
    
    
    with open(name,"r") as f:
        data = f.read()    
    
    lines = data.split('\n')

    lines = lines[1:len(lines)-1]

    for line in lines:
        a = line.split()
        x.append(float(a[cX])*scaleX)
        y.append(float(a[cY])*scaleY)