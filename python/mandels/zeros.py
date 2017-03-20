# -*- coding: utf-8 -*-
import math as m

"""
defincao da funcao

"""
def f(x,A):
    return m.tan(x) - A*x
    

"""
defincao da funcao

"""
def Df(x,A):
    sec = 1.0/m.cos(x)    
    
    return sec**2 - A

def bissecao(x1,x2,A,tol,maxIt):    
    
    it   = 0
    stop = False
    while not stop and it < maxIt:
        f1 = f(x1,A)
        xm = 0.5*(x1+x2)
        fm = f(xm,A)
            
        if fm*f1 > 0 :
            x1 = xm
        else:
            x2 = xm
               
        if m.fabs(f(xm,A)) < tol:
            stop = True
                
        it+=1   
 
    return xm,stop     
    
def newton(x,A,tol,maxIt):    
    
    it   = 0
    stop = False
    while not stop and it < maxIt:
        x -= f(x,A)/Df(x,A)
        if m.fabs(f(x,A)) < tol:
            stop = True
       
        it+=1   
 
    return x,stop     
   
"""
encontra os zeros da funcao f
"""
def zeros(A,root,nRoots):
    
    eps   = 1.e-07
    tol   = 1.e-07  
    maxIt = 1000000
    
    for i in range (0,nRoots): 
           
        if  i == 0 :      
            x1 = m.pi/4
            x2 = m.pi/2 - eps
        else:
            x1   = i*m.pi - m.pi/2 + eps
            x2   = i*m.pi + m.pi/2 - eps
       
        xm, stop   = bissecao(x1,x2,A,tol,maxIt)
#       xm, stop   = newton(x2,A,tol,maxIt)
        
        if  not stop :          
          print ("Numero maximo de it excedido na raiz %i: "%i)
          print ("f(%f) = %e"%(xm,m.fabs(f(xm,A))))          
          exit(0)
        
        
        root.append(xm) 