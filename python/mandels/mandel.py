# -*- coding: utf-8 -*-

import math as m
import myPlot 
import tempoNormalizado as t
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
    
    eps   = 1.e-08
    tol   = 1.e-08  
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

"""
Solucao analica das poros pressoes

entrada:
alphan - 
x      - pontos onde seram plotados (adimensional)
t      - tempo andimensional 

saida:
y      - solucao nos pontos x[] no tempo t

"""  
def analiticoX (y,alpha,x,t,scale):




    for xi in x:
        s = 0.0
        for alphai in alpha:
        
            tmp1 = (m.cos(alphai*xi)-m.cos(alphai))*m.sin(alphai)
            tmp2 = alphai - m.sin(alphai)*m.cos(alphai) 
            tmp3 = m.exp(-t*alphai**2)    
            s += (tmp1/tmp2)*tmp3 
        
        y.append(2*s*scale)

"""
Solucao analica dos delocamentos ux

entrada:
alphan   - 
x        - pontos onde seram plotados (adimensional)
t        - tempo andimensional 
w        - forca aplicada
modE     - modulo de elasticidade
poisson  - coeficiente de poisson 
poissonU - coeficiente de poisson nao drenado
L        - metade do comprimento (2a) 

saida:
y      - solucao nos pontos x[] no tempo t

"""  
        
def analiticoUx (y,alpha,x,t,w,modE,poisson,poissonU,L):

    nu     = 0.5*modE/(1+poisson)  

    for xi in x:
        s1 = 0.0
        s2 = 0.0
        for alphai in alpha:
        
            tmp1 = m.sin(alphai)*m.cos(alphai)
            tmp2 = alphai - m.sin(alphai)*m.cos(alphai) 
            tmp3 = m.exp(-t*alphai**2)  
            tmp4 = m.cos(alphai)*m.sin(alphai*xi)
            s1  += (tmp1/tmp2)*tmp3 
            s2  += (tmp4/tmp2)*tmp3 
        
        tmp5 = ((w*poisson)/(2.0*nu*L) - (w*poissonU)/(nu*L)*s1)*L*xi\
             + (w/nu)*s2
        
        y.append(tmp5)        

def readCsv(x,y,step,cX,cY,scaleX,scaleY):
    
    
    name = "csv/mandel_line."+str(step)+".csv"   
    
    with open(name,"r") as f:
        data = f.read()    
    
    lines = data.split('\n')

    lines = lines[1:len(lines)-1]

    for line in lines:
        a = line.split(',')
        x.append(float(a[cX])*scaleX)
        y.append(float(a[cY])/scaleY)
    
# dimensao
L = 10.0

modE    = 20.0 
poisson =  0.2



# modulo de viscosidade volumetrica
K       = (modE*poisson)/((1.0+poisson)*(1.0-2.0*poisson))\
        + (1.0/3.0)*(modE/(1.0+poisson))

#modulo de biot e coeficiente de biot
bioMod  = 1.0e+04
bioCoef = 1.0   

#fluido
k       = 1.0e-9
g       = 10.0
rof     = 1.0e+3

# calculo de cv


Cv = t.cv(k,rof,g,bioMod,bioCoef,modE,poisson)

A, B, poissonU = t.prop(modE,poisson,bioMod,bioCoef)

root = []

#
x    = [0.0]
y    = []
nDiv = 100 
dx   = 1.0/nDiv

for i in range(1,nDiv+1):
    xNew = dx + x[i-1]    
    x.append(xNew) 
    


zeros(A,root,50)

#
print ('B =',B,'Cv =',Cv,'vu',poissonU )


# escala
scaleP = B*(1.0+poissonU)/3.0
scaleX = 1.0/L

#
yn = []
xn = []

dt = 1.0e-3

print ('dt',t.deltat(L,Cv,dt))

color = ['red','blue','black','green']

# poropressao

myPlot.param['legenda'] = True 
myPlot.param['xrange']  = True 
myPlot.param['yrange']  = True 
xleg  = (1.01,0.95) 
rx=[0.0,1.0]
ry=[0.0,1.4]   
fileName = 'mandel_pressao.png' 
eixos=['Normalized distance x','Normalized pressure p']
myPlot.grafico(eixos,'Poro-pressure',0)

i = 0
for istep in [1,10,100,1000]:
    
    t = istep*dt    
    print (t)
#   tv = Cv/(L**2)*t
    analiticoX(y,root,x,t,1.0)
    
    legenda = 'Analytical(tv='+str(t)+')'
    myPlot.serie(x,y,'','-',legenda,color[i],myPlot.param,rx,ry,xleg)
   
    legenda = 'Numerical (tv='+str(t)+')'
    readCsv(xn,yn,istep,6,3,scaleX,scaleP)
    myPlot.serie(xn,yn,'o','--',legenda,color[i],myPlot.param,rx,ry,xleg)
#    
    i+=1
    xn.clear() 
    yn.clear() 
    y.clear()
 
myPlot.plot(False,fileName) 

   
myPlot.param['legenda'] = True 
myPlot.param['xrange']  = False 
myPlot.param['yrange']  = False 
xleg  = (1.01,0.95) 
rx=[0.0,1.0]
ry=[0.0,1.4]
myPlot.grafico(eixos,'Deslocamentos',2)   
fileName = 'mandel_pressao.png' 
eixos=['Normalized distance x','u']    
    
    
#deslocamento ux
i = 0
for istep in [1,10,100,1000]:
    
    t = istep*dt    
    print (t)
#   tv = Cv/(L**2)*t
    analiticoUx(y,root,x,t,1.0,modE,poisson,poissonU,L)
    
    legenda = 'Analytical(tv='+str(t)+')'
    myPlot.serie(x,y,'','-',legenda,color[i],myPlot.param,rx,ry,xleg)
   
    legenda = 'Numerical (tv='+str(t)+')'
    readCsv(xn,yn,istep,6,0,scaleX,1)
    myPlot.serie(xn,yn,'o','',legenda,color[i],myPlot.param,rx,ry,xleg)
#    
    i+=1
    xn.clear() 
    yn.clear() 
    y.clear()
    
myPlot.plot(False,fileName)   
