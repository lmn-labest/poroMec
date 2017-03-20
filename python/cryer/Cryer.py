# -*- coding: utf-8 -*-
import math as m
import myPlot 
import readFile as r 
import tempoNormalizado as tt
import zeros as z


def analiticoPc(alpha,v,x,t,scale):

    c1 = 4*(1-2*v)*(1-v) 
    c2 = 2*(1-2*v)*(1+v) 
    c3 = (1-v)**2 

    
    s = 0.0
    if x == 0.0:           
        for alphai in alpha:        
            tmp1 = m.sin(alphai)-alphai
            tmp2 = (c2-c3*alphai**2)*m.sin(alphai)
            tmp3 = m.exp(-t*alphai**2)    
            s   += (tmp1/tmp2)*tmp3 
        
        return c1*s
        
    else:
         for alphai in alpha:        
             tmp1 = 1.0 - m.sin(x*alphai)/(m.sin(alphai)*x)
             tmp2 = (c2-c3*alphai**2)
             tmp3 = m.exp(-t*alphai**2)    
             s   += (tmp1 /tmp2)*tmp3
             
         return c1*s    
      
    
def timePlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,R,Cv,dt): 
    global nZeros

   
    tn   = [0.0]
    nDiv = 1000
    dt   = 1/nDiv
    for i in range(1,nDiv+1):
        tNew = dt + tn[i-1]
        tn.append(tNew) 

    neta = tt.propCryer(modE,poisson)

    root = []
    y    = []
    t1   = []
# .....................................................................

# ... erros
    erros = (1,10,100,1000)
   
    print(neta) 
# ... calculo do coeficientes
    z.zeros2(z.fCryer,neta,root,nZeros)
#   print (root)
#
    print ('neta =',neta)

# escala
    scaleP    = F/bioCoef
    scaleX    = 1.0/R
    st        = R**2/Cv
    segToDay  = 1.0/(24.0*3600.0)
   
    color = ('red','blue','black','green')
# poropressao normalizada
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['logx']    = False
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.95) 
    rx    = (0.0,1.0)
    ry     =(0.0,1.5)   
    fileName = 'pressure_cryer_norm_time.png' 
    ex    = r'$\bar t$'
    ey    = r'$\bar p$'
    eixos =(ex,ey)
    eixos    =(ex,ey)
    myPlot.grafico(eixos,'Cryer\'s problem',4)
  
 
# lendo do arquivo
    name = 'plotTime/cryer_up_node_1530.txt' 
    xNum = []
    yNum = []
    r.readTxt(name,xNum,yNum,1,5,1/st,1/scaleP)

# solucao analitica
    
    for t in tn: 
        yn = analiticoPc(root,poisson,0.0,t,1)
        y.append(yn)    

    leg = "Analytical"
    myPlot.serie(tn,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 

    myPlot.plot(True,fileName)    
   
# poropressao 
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['logx']    = False
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.95) 
    rx    = (0.0,25.0)
    ry     =(0.0,3.0)   
    fileName = 'pressure_cryer_time.png'  
    eixos    =('t (day)','p(MPa)')
    myPlot.grafico(eixos,'Cryer\'s problem',5)
# lendo do arquivo
    name = 'plotTime/cryer_up_node_1530.txt' 
    xNum = []
    yNum = []
    r.readTxt(name,xNum,yNum,1,5,segToDay,1.0)

# solucao analitica
    y.clear()
    for t in tn: 
        yn = analiticoPc(root,poisson,0.0,t,1)
        y.append(yn)    

    t1.clear()
    for i in range(len(y)) :
        t1.append((segToDay*st)*tn[i])
        y[i]  = scaleP*y[i]

    leg = "Analytical"
    myPlot.serie(t1,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 

    myPlot.plot(True,fileName) 

# calculo do erro
    print (len(xNum))
    for step in erros:
        tReal = segToDay*(R**2/Cv)*tn[step] 
        error = 100*abs(yNum[step]-y[step])/abs(y[step])
        print ('uy tn = %f t = %f num = %f ex = %f erro() = %f'
              %(tn[step],tReal,yNum[step],y[step], error ))
    

def xPlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,R,Cv,dt):  
 

    global nZeros

    neta = tt.propCryer(modE,poisson)

    root = []
    y    = []
# ...
    nDiv = 1000
    xn   = [0.0]
    dx   = 1.0/nDiv
    for i in range(1,nDiv+1):
        xNew = dx + xn[i-1]    
        xn.append(xNew) 
# .....................................................................
        
# ... calculo do coeficientes
    z.zeros2(z.fCryer,neta,root,nZeros)

    print ('neta =',neta)

# escala
    scaleP    = F/bioCoef
    scaleX    = 1.0/R
    st        = R**2/Cv
    segToDay  = 1.0/(24.0*3600.0)
   
    color = ('red','blue','black','green')

#
    yNum = []
    xNum = []


# poropressao

    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['logx']    = False
    myPlot.param['yrange']  = True 
    xleg  = (1.01,0.70) 
    rx=(0.0,1.0)
    ry=(0.0,1.4)   
    fileName = 'cryer_pressao_line.png' 
    ex    = r'$\bar x$'
    ey    = r'$\bar p$'
    eixos =(ex,ey)
    myPlot.grafico(eixos,'Cryer\'s problem',0)

    i = 0
    et    = r'$\bar t = $'
    for istep in [1,10,100,500]:
        
        t = istep*dt    
        print ('tc',t)

        for xi in xn:
            yi = analiticoPc(root,poisson,xi,t,1.0)
            y.append(yi)
    
        legenda = 'Analytical ('+et+str(t)+')'
        myPlot.serie(xn,y,'','-',legenda,color[i],myPlot.param,rx,ry,xleg)
   
   
        name = "csv/cryer_line."+str(istep)+".csv"  
        legenda = 'Numerical ('+et+str(t)+')'
        r.readCsv(name,xNum,yNum,istep,30,13,scaleX,scaleP)
        myPlot.serie(xNum,yNum,'','--',legenda,color[i],myPlot.param,rx,ry,xleg)
#    
        i+=1
        xNum.clear() 
        yNum.clear() 
        y.clear()
 
    myPlot.plot(True,fileName) 

          
nZeros =500

def main():
    
# dimensao
    R = 1.5

    modE    =  1.0e+01 
    poisson =  0.10 

#modulo de biot e coeficiente de biot
    bioMod  = 1.0e+04
    bioCoef = 1.0

#fluido
    k       = 1.0e-9
    g       = 10.0
    rof     = 1.0e+3
#
    F       = 2.0000 
#   
    # calculo de cv
    Cv = tt.cv(k,rof,g,bioMod,bioCoef,modE,poisson)
#    
    dt = 1.0e-3
    print ('dt',tt.deltat(R,Cv,dt))
    tc = (R**2)/Cv
    print ('tc',tc,tc/(24*3600))   
#    
    xPlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,R,Cv,dt) 
    timePlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,R,Cv,dt)
  
    print ('Done')  


main()  