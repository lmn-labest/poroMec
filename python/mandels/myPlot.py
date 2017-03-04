# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

param = {}

param['legenda'] = False
param['logx']    = False
param['logy']    = False
param['log']     = False
param['xrange']  = False
param['yrange']  = False

def grafico(eixos,legenda,icod):

    fig = plt.figure(icod,figsize=(6,6))
    plt.xlabel(eixos[0])
    plt.ylabel(eixos[1])
    plt.title(legenda, fontsize="medium")   
    plt.grid(True)
    
    return fig

def serie(x,y,marker,line,legenda,color,param,rx,ry,xlegenda):
    
    
    if param['xrange']:    
        plt.xlim(rx)
    if param['yrange']:
        plt.ylim(ry)  
   
    if param['logx']:
        plt.semilogx(x,y,linestyle='-',label=legenda)    
    elif param['logy']:
        plt.semilogy(x,y,linestyle='-',label=legenda)    
    elif param['log']:
        plt.loglog(x,y,linestyle='-',label=legenda) 
    else:
        plt.plot(x,y,marker,linestyle=line,label=legenda,color=color)
   
   
    if param['legenda']:   
      plt.legend(bbox_to_anchor=xlegenda, loc=2, borderaxespad=0.
                ,prop={'size':8})


def plot(disco,file):
    
    if disco:        
        plt.savefig(file,ext="png", close=True, verbose=True
               ,bbox_inches='tight',dpi=1000) 
    else:
        plt.show()
