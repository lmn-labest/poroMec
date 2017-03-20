# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 17:23:35 2016

@author: henrique
"""
import math as m
import readFile as r 
import myPlot
import tempoNormalizado as tt
    
def tc(h,mBiot,cBiot,modE,possion,k,rof,g):
    scale = 1.0e-6
    gamma = rof*g*scale
    v   = possion    
    tmp = 1.0/mBiot + (cBiot*cBiot*(1+v)*(1-2.0*v))/(modE*(1.0-v))        
    
    return (gamma*h/k)*tmp

"""
Solucao analica das poros pressoes

entrada:
alphan - 
x      - pontos onde seram plotados (adimensional)
t      - tempo andimensional 
n      - numero de termpos no somatorio
scaleY - fator de escala 
saida:
pressao(x,t)

"""  
def analiticoPc (xn,tn,n,scaleY):

    s = 0.0
    for i in range(n):
        tmp1 = (2.0*i+1.0)*m.pi
        tmp2 = m.sin(0.5*tmp1*xn)
        tmp3 = m.exp(-tn*(0.5*tmp1)**2)
        s   += (4.0/tmp1)*tmp2*tmp3
    return s
  
"""
Fator de escala para normalizacao da pressao
(pag. 126 - PoroMechanics - Oliver Coussy)
entrada:
t        - tempo andimensional 
n        - numero de termpos no somatorio
scaleY   - fator de escala 
modE     - modulo de elasticidade
poisson  - coeficiente de possion
bM       - modulo de Biot
bC       - coeficiente de Biot
F        - modulo da forca aplicada
saida:
recalque(t)

"""      
def analiticoS(tn,n,scaleY,modE,poisson,bM  ,bC   ,F):
  
    nu = 0.5*modE/(1.0+poisson)
    la = modE*poisson/((1.0+poisson)*(1.0-2.0*poisson))
  
    s = 0.0
    
    for i in range(n):
        tmp1 = 2.0*i+1.0
        tmp2 = m.exp(-tn*(0.5*tmp1*m.pi)**2)
        s   += tmp2/tmp1**2

    tmp1 = F/(la+2.0*nu)
    tmp2 = (bM*bC**2)/(la+2.0*nu+bM*bC**2)
    tmp3 = 1.0 - (8.0/m.pi**2)*tmp2*s

    return tmp1*tmp3     
    
"""
Fator de escala para normalizacao da pressao
(pag. 126 - PoroMechanics - Oliver Coussy)
entrada:
modE     - modulo de elasticidade
poisson  - coeficiente de possion
bMod     - modulo de Biot
bCoef    - coeficiente de Biot
F        - modulo da forca aplicada
saida:
fator de escala

"""      
    
def scaleP(modE,poisson,bMod,bCoef,F):
    
    nu = 0.5*modE/(1.0+poisson)
    la = modE*poisson/((1.0+poisson)*(1.0-2.0*poisson))    
    
# modulo de viscosidade volumetrica
    K = la + 1.5*nu    

# modulo de viscosidade volumetrica nao drenado
    Ku  = K + bMod*bCoef**2    
#
    tmp = Ku + (4.0/3.0)*nu
    
#   
    return (bCoef*bMod*F)/tmp

def timePlot(modE,poisson,bMod,bCoef,k,rof,g,H,F,cv):
    
    global N    
    
#
    xn   = [0.0]
    nDiv = 20 
    dx   = 1.0/nDiv
    
    for i in range(1,nDiv+1):
        xNew = dx + xn[i-1]    
        xn.append(xNew)
    
    tn   = [0.0]
    nDiv = 6000
    dt   = 6.0/nDiv
    for i in range(1,nDiv+1):
        tNew = dt + tn[i-1]
        tn.append(tNew) 

# fatores de escala 
    sp = scaleP(modE,poisson,bMod,bCoef,F)
    st = H**2/cv 
    sH = 1.0/H
   
    color = ('red','blue','black','green')
# poropressao normalizada
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.95) 
    rx    = (0.0,6.0)
    ry     =(0.0,1.1)   
    fileName = 'pressure_norm.png' 
    eixos    =('t/tc','Normalized p')
    myPlot.grafico(eixos,'Fluid Pressure at Base',0)
  
 
# lendo do arquivo
    name = 'plotTempo/solo2_up_node_441.txt' 
    xNum = []
    yNum = []
    r.readTxt(name,xNum,yNum,1,5,1/st,1/sp)

# solucao analitica
    y = []   
    for t in tn: 
        yn = analiticoPc(1.0,t,N,1.0)
        y.append(yn)    

    leg = "Analytical"
    myPlot.serie(tn,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 

    myPlot.plot(True,fileName)    
    
# poropressao
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = False
    myPlot.param['yrange']  = False
    xleg  = (0.65,0.95) 
    rx    = (0.0,25.0)
    ry     =(0.0,4.0)   
    fileName = 'pressure.png' 
    eixos    =('t (Day)','p (MPa)')
    myPlot.grafico(eixos,'Terzaghi\'s Problem',1)
      
    segToDay = 1.0/(24.0*3600.0)
    
    t1 = []
    for i in range(len(y)) :
        t1.append((segToDay*st)*tn[i])
        y[i]  = sp*y[i] 

# lendo do arquivo    
    name = 'plotTempo/solo2_up_node_441.txt'  
    xNum.clear()
    yNum.clear()
    r.readTxt(name,xNum,yNum,1,5,segToDay,1.0)
      
    leg = "Analytical"
    myPlot.serie(t1,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 

    myPlot.plot(True,fileName)     

# calculo do erro
    print (len(xNum))
    for step in [100,1000,6000]:
        tReal = segToDay*(H**2/cv)*tn[step] 
        error = 100*abs(yNum[step]-y[step])/abs(y[step])
        print ('p tn = %f t = %f num = %f ex = %f erro() = %f'
              %(tn[step],tReal,yNum[step],y[step], error ))
 

# deslocamento
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.25) 
    rx    = (0.0,6.0)
    ry     =(0.0,0.2)   
    fileName = 'desloc_norm.png' 
    eixos    =('t/tc','Normalized |uz|')
    myPlot.grafico(eixos,'Terzaghi\'s Problem',2)

# lendo do arquivo    
    name = 'plotTempo/solo2_up_node_9261.txt' 
    xNum.clear()
    yNum.clear()    
    r.readTxt(name,xNum,yNum,1,4,1/st,-sH)    
    
# solucao analitica   
    y.clear()   
    for t in tn: 
        yn = analiticoS(t,N,1.0,modE,poisson,bMod,bCoef,F)
        y.append(yn)  
    
    leg = "Analytical"
    myPlot.serie(tn,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 
    
    # deslocamento
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = False
    myPlot.param['yrange']  = False
    xleg  = (0.65,0.25) 
    rx    = (0.0,25.0)
    ry     =(0.0,0.2)   
    fileName = 'desloc.png' 
    eixos    =('t (day)','s (m)')
    myPlot.grafico(eixos,'Terzaghi\'s Problem',3)

# lendo do arquivo    
    name = 'plotTempo/solo2_up_node_9261.txt' 
    xNum.clear()
    yNum.clear()    
    r.readTxt(name,xNum,yNum,1,4,segToDay,-sH)    
    
# solucao analitica   
    y.clear()   
    for t in tn: 
        yn = analiticoS(t,N,1.0,modE,poisson,bMod,bCoef,F)
        y.append(yn)  
    
    t1.clear()
    for i in range(len(y)) :
        t1.append((segToDay*st)*tn[i])
        y[i]  = y[i]     
    
    
    leg = "Analytical"
    myPlot.serie(t1,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 


# calculo do erro
    print (len(xNum))
    for step in [100,1000,6000]:
        tReal = segToDay*(H**2/cv)*tn[step] 
        error = 100*abs(yNum[step]-y[step])/abs(y[step])
        print ('u tn = %f t = %f num = %f ex = %f erro() = %f'
              %(tn[step],tReal,yNum[step],y[step], error ))
            
       
    myPlot.plot(True,fileName)   


def xPlot(modE,poisson,bMod,bCoef,k,rof,g,H,F,cv):
    
    global N    
    
    
    y    = []
    xn   = []
    
    xp   = []    
    
    xNum = []
    yNum = []
    
# ... fatores de escala 
    sP = scaleP(modE,poisson,bMod,bCoef,F)
    sH = 1.0/H
# .....................................................................

# ...    
    xn   = [0.0]
    nDiv = 100 
    dx   = 1.0/nDiv
  
    for i in range(1,nDiv+1):
        xNew = dx + xn[i-1]    
        xn.append(xNew)    
# .....................................................................

# ... 
    color = ('red','blue','black','green')
# .....................................................................

# ... 
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True 
    myPlot.param['yrange']  = True 
    xleg  = (1.01,0.70) 
    rx=[0.0,1.0]
    ry=[0.0,1.1]
    eixos=('z/H','Normalized p')
    myPlot.grafico(eixos,'Terzaghi\'s Problem',4)   
# .....................................................................

# ...            
    i  = 0
    dt = 0.001
    for istep in [1,10,100,1000]:
    
        t = istep*dt    
        print (t)
#   
        for xi in xn:
            yi = analiticoPc (xi,t,N,1.0)
            y.append(yi)
        
        
        for xx in xn:
            xp.append(1-xx)                   
   
        legenda = 'Analytical  (tn='+str(t)+')'
        myPlot.serie(xp,y,'','-',legenda,color[i],myPlot.param,rx,ry,xleg)
   
        name = "csv/solo2_line."+str(istep)+".csv"  
        legenda = 'Numerical (tn='+str(t)+')'
        r.readCsv(name,xNum,yNum,istep,30,13,sH,sP)
        
        xp.clear()
        for xx in xNum:
            xp.append(1-xx)    
        
        myPlot.serie(xp,yNum,'','--',legenda,color[i],myPlot.param,rx,ry,xleg)
#    
        i+=1
        xNum.clear() 
        yNum.clear() 
        xp.clear()   
        y.clear()
# .....................................................................
        
    fileName = 'pressao_line.png' 
    myPlot.plot(True,fileName) 


N = 1000

def main():

#prop
    modE    = 16.68
    poisson = 0.3369
    bMod    = 5.91074e+4
    bCoef   = 0.9989
    k       = 1.1e-9
    rof     = 1.0e+3
    g       = 10.0
    H       =  1.0
    F       =  3.9
#    
    dt = 1.0e-3
    Cv = tt.cv(k,rof,g,bMod,bCoef,modE,poisson)
    print ('dt',tt.deltat(H,Cv,dt))
    tc = (H**2)/Cv
    print ('tc',tc,tc/(24*3600))   
    
    print ("timePlot")
    timePlot(modE,poisson,bMod,bCoef,k,rof,g,H,F,Cv) 
    print ("xPlot")
    xPlot(modE,poisson,bMod,bCoef,k,rof,g,H,F,Cv)
    print ("Done")
 
if __name__ == "__main__":
    
    main()	
	