# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 17:23:35 2016

@author: henrique
"""
import myPlot
import readFile as r 
import tempoNormalizado as tt


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

    name1 =['exemplo1_mesh0_t2_pl_up_node_1.txt'
	       ,'exemplo1_mesh0_t2_v_pl_up_node_1.txt']

 
           
    name2 =['exemplo1_mesh0_t2_pl_up_node_41.txt'
	       ,'exemplo1_mesh0_t2_v_pl_up_node_41.txt']

           
    legenda=['Constante','kc-hs','hs','kc-hs']

    t      = []
    desloc = []
    pres   = []     

# escala
    st      = H**2/cv     
    sx      =-1.0 
    sP      = scaleP(modE,poisson,bMod,bCoef,F)
    
# color
    color =('black','red','blue','green')   
# linhas
    l     =('-','--','-','-')    
    
# deslocamento
    myPlot.param['legenda'] = True 
    myPlot.param['logx']    = False
    myPlot.param['logy']    = False
    myPlot.param['log']     = False
    myPlot.param['xrange']  = True 
    myPlot.param['yrange']  = True 
    xleg  = (0.7,0.25) 
    rx=(0.0,6.0)
    ry=(0.0,0.1)   
    fileName = 'desloc.png' 
    ex   = r'$\bar t$'
    ey   = r'$\Delta h/h$'
    eixos=[ex, ey]
    myPlot.grafico(eixos,'Time evolution of the surface subsidence',0)
    it = 0
    for name in name1:
        with open(name,"r") as f:
            data = f.read()

        lines  = data.split('\n')             

        for i in range(1,len(lines)-1):
            line    = lines[i].split()     
            t.append(float(line[1])/st) 
            desloc.append(float(line[4])*sx)    
 
        myPlot.serie(t,desloc,'',l[it],legenda[it],color[it],myPlot.param,rx,ry,xleg)
      
        it = it + 1  
        t.clear()
        desloc.clear()

    myPlot.plot(False,fileName)    

 # press
    myPlot.param['legenda'] = True 
    myPlot.param['logx']    = False
    myPlot.param['logy']    = False
    myPlot.param['xrange']  = True 
    myPlot.param['yrange']  = False 
    xleg  = (0.7,0.95) 
    rx=(0.0,6.0)
    ry=(0.0,0.5)   
    fileName = 'pres_time.png' 
    ex   = r'$\bar t$'
    ey   = r'$\bar p$'
    fileName = 'pres.png' 
    eixos=[ex,ey]
    myPlot.grafico(eixos,'Time evolution of the pressure in the base',1)
    it = 0
    for name in name2:
        with open(name,"r") as f:
            data = f.read()

        lines  = data.split('\n')          

        lines  = lines[1:len(lines)]

        for i in range(1,len(lines)-1):
            line    = lines[i].split()   
            t.append(float(line[1])/st)   
            pres.append(float(line[5])/sP)   

        myPlot.serie(t,pres,'',l[it],legenda[it],color[it],myPlot.param,rx,ry,xleg)
      
        it = it + 1  
        t.clear()
        pres.clear()

    myPlot.plot(False,fileName)    
    
def xPlot(modE,poisson,bMod,bCoef,k,rof,g,H,F,cv):
    
    global N    
     
    y    = []
    xn   = []
    
    xp   = []    
    
    xNum = []
    yNum = []
    
# ... fatores de escala 
    sP   = scaleP(modE,poisson,bMod,bCoef,F)
    st   = H**2/cv     
    sH   = 1.0/H
    sTd  = 3600.0
# .....................................................................

# ...    
# .....................................................................

# ... 
    color = ('red','blue','black','green')
# .....................................................................

# ... 
    myPlot.param['legenda'] = True 
    myPlot.param['xrange']  = False 
    myPlot.param['yrange']  = False
    xleg  = (1.01,0.60) 
    rx=[0.0,1.0]
    ry=[0.0,1.1]
    ex   = r'$\bar z$'
    ey   = r'$\bar p$'
    eixos=[ex, ey]
    myPlot.grafico(eixos,'Evolution of the pressure in the layer',4)   
# .....................................................................

# ...
    t = [0.0]
    for i in range(1,1300):
        if 0 <= i <= 100:
            dt = 900.0
        else:
            dt = 1800.0    
        
        t.append(t[i-1]+dt)  
# .....................................................................

# ...            
    i  = 0
    for istep in [20,200,600]:
    
        nameE= "csv/solo2_elastic."+str(istep)+".csv" 
        nameP= "csv/solo2_plastic_line."+str(istep)+".csv" 
        
        
        tt      = t[istep]/sTd
        legenda = 'Elastic (t ='+str(tt)+' h )'
        r.readCsvV2(nameE,xNum,yNum,istep,'"arc_length"','"pressao"',sH,sP)
        
        xp.clear()
        for xx in xNum:
            xp.append(1-xx)    
        
        myPlot.serie(xp,yNum,'','-',legenda,color[i],myPlot.param,rx,ry,xleg)
        
        xNum.clear() 
        yNum.clear() 
        xp.clear()   
        y.clear()        
        
        
        legenda = 'Plastic (t ='+str(tt)+' h )'
        r.readCsvV2(nameP,xNum,yNum,istep,'"arc_length"','"pressao"',sH,sP)
        
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
 
def main(): 
 
    modE    = 16.68
    poisson = 0.3369
    bMod    = 5.91074e+4
    bCoef   = 0.9989
    k       = 1.1e-9
    rof     = 1.0e+3
    g       = 10.0
    H       =  1.0
    F       =  1.5
#    
    dt = 1.0e-3
    Cv = tt.cv(k,rof,g,bMod,bCoef,modE,poisson)
    print ('dt',tt.deltat(H,Cv,dt))
    tc = (H**2)/Cv
    print ('tc',tc,tc/(24*3600))   
    
    print ("timePlot")
    timePlot(modE,poisson,bMod,bCoef,k,rof,g,H,F,Cv) 
#    print ("xPlot")
#    xPlot(modE,poisson,bMod,bCoef,k,rof,g,H,F,Cv)
#    print ("Done") 
 
if __name__ == "__main__":
    
    main()	
	