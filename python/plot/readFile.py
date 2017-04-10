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
     
def readCsvV2(name,x,y,step,strX,strY,scaleX,scaleY):
    
    with open(name,"r") as f:
        data = f.read()    
    
    lines = data.split('\n')

    nc = 0
    cX = 0
    cY = 0
    for word in  lines[0].split(',') :
        
        if word == strX:
            cX = nc
            
        elif word == strY:   
            cY = nc
            
        nc+=1 
    
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