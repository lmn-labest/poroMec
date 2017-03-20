# -*- coding: utf-8 -*-

import math as m
import myPlot 
import readFile as r 
import tempoNormalizado as tt
import zeros as z

"""
Solucao analica das poros pressoes

entrada:
alphan - 
x      - pontos onde seram plotados (adimensional)
t      - tempo andimensional 

saida:
pressao(x,t)

"""  
def analiticoPc(alpha,x,t,scale):

    s = 0.0
    for alphai in alpha:
        
        tmp1 = (m.cos(alphai*x)-m.cos(alphai))*m.sin(alphai)
        tmp2 = alphai - m.sin(alphai)*m.cos(alphai) 
        tmp3 = m.exp(-t*alphai**2)    
        s += (tmp1/tmp2)*tmp3 
        
    return 2*s*scale

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
        
def analiticoUx (alpha,x,t,w,modE,poisson,poissonU,L):

    nu     = 0.5*modE/(1+poisson)  

    s1 = 0.0
    s2 = 0.0
    
#    for alphai in alpha:
#        
#        tmp1 = m.sin(alphai)*m.cos(alphai)
#        tmp2 = alphai - m.sin(alphai)*m.cos(alphai) 
#        tmp3 = m.exp(-t*alphai**2)  
#        tmp4 = m.cos(alphai)*m.sin(alphai*x)
#        s1  += (tmp1/tmp2)*tmp3 
#        s2  += (tmp4/tmp2)*tmp3 
#        
#        tmp5 = ((w*poisson)/(2.0*nu) - ((w*poissonU)/nu)*s1)*L*x\
#             + ((w*L)/nu)*s2
#        
#    return tmp5

    for alphai in alpha:
        
        tmp1 = m.sin(alphai)*m.cos(alphai)
        tmp2 = alphai - m.sin(alphai)*m.cos(alphai) 
        tmp3 = m.exp(-t*alphai**2)  
        tmp4 = m.cos(alphai)*m.sin(alphai*x)
        s1  += (tmp1/tmp2)*tmp3 
        s2  += (tmp4/tmp2)*tmp3 
        
        tmp5 = (w/nu)*((0.5*poisson - poissonU*s1)*x+s2)
        
    return tmp5 

"""
Solucao analica dos delocamentos uy

entrada:
alphan   - 
x        - pontos onde seram plotados (adimensional)
t        - tempo andimensional 
w        - forca aplicada
modE     - modulo de elasticidade
poisson  - coeficiente de poisson 
poissonU - coeficiente de poisson nao drenado
L        - metade do comprimento (2a) 
H        - metade do espessura   (2b) 

saida:
y      - solucao nos pontos x[] no tempo t

"""  
def analiticoUy (alpha,y,t,w,modE,poisson,poissonU,L,H):

    nu     = 0.5*modE/(1+poisson)  

    s1 = 0.0
#    for alphai in alpha:
        
#        tmp1 = m.sin(alphai)*m.cos(alphai)
#        tmp2 = alphai - m.sin(alphai)*m.cos(alphai) 
#        tmp3 = m.exp(-t*alphai**2)  
#        s1  += (tmp1/tmp2)*tmp3 
        
#        tmp4 = ( (w*(1-poissonU)/nu)*s1 - w*(1-poisson)/(2*nu) )*H*y
             
        
#    return tmp4 
    for alphai in alpha:
        
        tmp1 = m.sin(alphai)*m.cos(alphai)
        tmp2 = alphai - m.sin(alphai)*m.cos(alphai) 
        tmp3 = m.exp(-t*alphai**2)  
        s1  += (tmp1/tmp2)*tmp3 
        
        tmp4 =  (w/nu)*(((1-poissonU)*s1 - 0.5*(1-poisson))*y)
             
        
    return tmp4


def xPlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,L,H,Cv,dt):  
    global nZeros

    A, B, poissonU = tt.prop(modE,poisson,bioMod,bioCoef)

    root = []
    y    = []
# ...
    nDiv = 100 
    xn   = [0.0]
    dx   = 1.0/nDiv
    for i in range(1,nDiv+1):
        xNew = dx + xn[i-1]    
        xn.append(xNew) 
# .....................................................................
        
# ... calculo do coeficientes
    z.zeros(A,root,nZeros)

#
    print ('B =',B,'Cv =',Cv,'vu',poissonU )


# escala
    scaleP   = B*(1.0+poissonU)/3.0
    scaleX   = 1.0/L

#
    yNum = []
    xNum = []

    color = ['red','blue','black','green']

# poropressao

    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True 
    myPlot.param['yrange']  = True 
    xleg  = (1.01,0.70) 
    rx=(0.0,1.0)
    ry=(0.0,1.4)   
    fileName = 'mandel_pressao_line.png' 
    eixos    =('t/tc','Normalized p')
    myPlot.grafico(eixos,'Mandel\'s problem',0)

    i = 0
    for istep in [1,10,100,1000]:
        
        t = istep*dt    
        print ('tc',t)

        for xi in xn:
            yi = analiticoPc(root,xi,t,1.0)
            y.append(yi)
    
        legenda = 'Analytical (tv='+str(t)+')'
        myPlot.serie(xn,y,'','-',legenda,color[i],myPlot.param,rx,ry,xleg)
   
   
        name = "csv/mandel_line."+str(istep)+".csv"  
        legenda = 'Numerical (tv='+str(t)+')'
        r.readCsv(name,xNum,yNum,istep,6,3,scaleX,scaleP)
        myPlot.serie(xNum,yNum,'','--',legenda,color[i],myPlot.param,rx,ry,xleg)
#    
        i+=1
        xNum.clear() 
        yNum.clear() 
        y.clear()
 
    myPlot.plot(True,fileName) 

#
def timePlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,L,H,Cv,dt): 
    global nZeros

#
    xn   = [0.0]
    nDiv = 20 
    dx   = 1.0/nDiv
    
    for i in range(1,nDiv+1):
        xNew = dx + xn[i-1]    
        xn.append(xNew)
    
    tn   = [0.0]
    nDiv = 3000
    dt   = 3.0/nDiv
    for i in range(1,nDiv+1):
        tNew = dt + tn[i-1]
        tn.append(tNew) 

    A, B, poissonU = tt.prop(modE,poisson,bioMod,bioCoef)

    root = []
    y    = []
    t1   = []
# .....................................................................

# ... erros
    erros = (10,1000,2000,3000)
   
     
# ... calculo do coeficientes
    z.zeros(A,root,nZeros)

#
    print ('B =',B,'Cv =',Cv,'vu',poissonU )

# escala
    scaleP    = B*(1.0+poissonU)/3.0
    scaleX    = 1.0/L
    scaleY    = 1.0/H
    st        = L**2/Cv
    segToDay  = 1.0/(24.0*3600.0)
   
    color = ('red','blue','black','green')
# poropressao normalizada
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.95) 
    rx    = (0.0,3.0)
    ry     =(0.0,1.2)   
    fileName = 'pressure_norm_time.png' 
    eixos    =('t/tc','Normalized p')
    myPlot.grafico(eixos,'Mandel\'s problem',4)
  
 
# lendo do arquivo
    name = 'plotTime/mandel_hexa_up_node_1619.txt' 
    xNum = []
    yNum = []
    r.readTxt(name,xNum,yNum,1,5,1/st,1/scaleP)

# solucao analitica
    for t in tn: 
        yn = analiticoPc(root,0.0,t,1)
        y.append(yn)    

    leg = "Analytical"
    myPlot.serie(tn,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 

    myPlot.plot(True,fileName)    
    
# poropressao 
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.95) 
    rx    = (0.0,1600.0)
    ry     =(0.0,0.6)   
    fileName = 'pressure_time.png'  
    eixos    =('t (day)','p(MPa)')
    myPlot.grafico(eixos,'Mandel\'s problem',5)
# lendo do arquivo
    name = 'plotTime/mandel_hexa_up_node_1619.txt' 
    xNum = []
    yNum = []
    r.readTxt(name,xNum,yNum,1,5,segToDay,1.0)

# solucao analitica
    y.clear()
    for t in tn: 
        yn = analiticoPc(root,0.0,t,1)
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
    for step in erros:
        tReal = segToDay*(L**2/Cv)*tn[step] 
        error = 100*abs(yNum[step]-y[step])/abs(y[step])
        print ('p tn = %f t = %f num = %f ex = %f erro() = %f'
              %(tn[step],tReal,yNum[step],y[step], error ))
    

# deslocamento
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.9) 
    rx    = (0.0,6.0)
    ry     =(0.0,0.2)   
    fileName = 'deslocx_norm_time.png' 
    eixos    =('t/tc','ux/L')
    myPlot.grafico(eixos,'Mandel\'s problem',6)

# lendo do arquivo    
    name = 'plotTime/mandel_hexa_up_node_1366.txt' 
    xNum.clear()
    yNum.clear()    
    r.readTxt(name,xNum,yNum,1,2,1/st,scaleX)    
    
# solucao analitica   
    y.clear()   
    for t in tn: 
        yi = analiticoUx(root,1.0,t,F,modE,poisson,poissonU,L)
        y.append(yi)  

    leg = "Analytical"
    myPlot.serie(tn,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 
 
    myPlot.plot(True,fileName)
    
  
# deslocamento
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.9) 
    rx    = (0.0,0.1)
    ry     =(0.0,0.2)   
    fileName = 'deslocx_time.png' 
    eixos    =('t(day)','ux(m)')
    myPlot.grafico(eixos,'Mandel\'s problem',7)

# lendo do arquivo    
    name = 'plotTime/mandel_hexa_up_node_1366.txt' 
    xNum.clear()
    yNum.clear()    
    r.readTxt(name,xNum,yNum,1,2,segToDay,1.0)    
    
# solucao analitica   
    y.clear()   
    for t in tn: 
        yi = analiticoUx(root,1.0,t,F,modE,poisson,poissonU,L)
        y.append(yi)  

    t1.clear()
    for i in range(len(y)) :
        t1.append((segToDay*st)*tn[i])
        y[i]  = y[i]/scaleX    
   
    leg = "Analytical"
    myPlot.serie(t1,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 
 
    myPlot.plot(True,fileName)
    
# calculo do erro
    print (len(xNum))
    for step in erros:
        tReal = segToDay*(L**2/Cv)*tn[step] 
        error = 100*abs(yNum[step]-y[step])/abs(y[step])
        print ('uz tn = %f t = %f num = %f ex = %f erro() = %f'
              %(tn[step],tReal,yNum[step],y[step], error ))
# deslocamento
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.9) 
    ry    = (-0.025,0.0)
    rx     =(0.0,3.0)   
    fileName = 'deslocy_norm_time.png' 
    eixos    =('t/tc','uy/H')
    myPlot.grafico(eixos,'Mandel\'s problem',8)

# lendo do arquivo    
    name = 'plotTime/mandel_hexa_up_node_1628.txt' 
    xNum.clear()
    yNum.clear()    
    r.readTxt(name,xNum,yNum,1,3,1/st,scaleY)    
    
# solucao analitica   
    y.clear()   
    for t in tn: 
        yi = analiticoUy(root,0.5,t,F,modE,poisson,poissonU,L,H)
        y.append(yi)  

    leg = "Analytical"
    myPlot.serie(tn,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 
 
    myPlot.plot(True,fileName)

  
# deslocamento
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = True
    myPlot.param['yrange']  = True
    xleg  = (0.65,0.9) 
    ry    = (-1.3,0.0)
    rx     =(0.0,1600)   
    fileName = 'deslocy_time.png' 
    eixos    =('t (day)','uy(m)')
    myPlot.grafico(eixos,'Mandel\'s problem',9)

# lendo do arquivo    
    name = 'plotTime/mandel_hexa_up_node_1628.txt' 
    xNum.clear()
    yNum.clear()    
    r.readTxt(name,xNum,yNum,1,3,segToDay,1.0)    
    
# solucao analitica   
    y.clear()   
    for t in tn: 
        yi = analiticoUy(root,0.5,t,F,modE,poisson,poissonU,L,H)
        y.append(yi)  

    t1.clear()
    for i in range(len(y)) :
        t1.append((segToDay*st)*tn[i])
        y[i]  = y[i]/scaleY
   
    leg = "Analytical"
    myPlot.serie(t1,y,'','-',leg,color[0],myPlot.param,rx,ry,xleg) 
    
    leg = "Numerical"
    myPlot.serie(xNum,yNum,'','--',leg,color[2],myPlot.param,rx,ry,xleg) 
 
    myPlot.plot(True,fileName)
    
# calculo do erro
    print (len(xNum))
    for step in erros:
        tReal = segToDay*(L**2/Cv)*tn[step] 
        error = 100*abs(yNum[step]-y[step])/abs(y[step])
        print ('uy tn = %f t = %f num = %f ex = %f erro() = %f'
              %(tn[step],tReal,yNum[step],y[step], error ))

nZeros = 100

def main():
    
# dimensao
    L = 10.0
    H = 50.0

    modE    = 20.0 
    poisson =  0.2

#modulo de biot e coeficiente de biot
    bioMod  = 1.0e+04
    bioCoef = 0.5   

#fluido
    k       = 1.0e-9
    g       = 10.0
    rof     = 1.0e+3
#
    F       = 1.0    
#   
    # calculo de cv
    Cv = tt.cv(k,rof,g,bioMod,bioCoef,modE,poisson)
#    
    dt = 1.0e-3
    print ('dt',tt.deltat(L,Cv,dt))
    tc = (L**2)/Cv
    print ('tc',tc,tc/(24*3600))   
#    
    xPlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,L,H,Cv,dt) 
    timePlot(modE,poisson,bioMod,bioCoef,k,g,rof,F,L,H,Cv,dt)
  
    print ('Done')  


main()  