# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 17:25:57 2017

@author: henrique
"""
def prop(modE,poisson,bioMod,bioCoef):

# modulo de viscosidade volumetrica
    K       = (modE*poisson)/((1.0+poisson)*(1.0-2.0*poisson))\
            + (1.0/3.0)*(modE/(1.0+poisson))

# modulo de viscosidade volumetrica nao drenado
    Ku      = K + bioMod*bioCoef**2

    H       = ((1.0-2.0*poisson)*bioMod*(bioCoef**2))/(3.0*Ku)

#poisson nao drenando
    poissonU = (H + poisson)/(1.0-H)

#B
    B = bioMod*bioCoef/Ku 


    A = (1.0 - poisson)/(poissonU-poisson)

    return A, B, poissonU


"""
Calculo de deltaT real considerando um deltaT normalizado        

entrada:
L        - comprimento caracteristico
cv       - coeficente de consolidacao 
dtv      - deltat normalizado                          


saida:
dt - real              

"""

def deltat(L,cv,dtv):
    
    return ((L**2)/cv)*dtv
  
"""
Calcula do coeficiente de adensamento        

entrada:
k        - coeficiente de permeabilidade
rof      - massa especifica da agua   
g        - modulo da gravidade
bMod     - modulo de Biot
bCoef    - coeficiente de Biot
modE     - modulo de elasticidade
poisson  - coeficiente de possion

saida:
retorn cv     

"""  
def cv(k,rof,g,bMod,bCoef,modE,poisson):

# fator de escala para MPa    
    scale = 1.0e-6    

    nu    = 0.5*modE/(1.0+poisson)
    la    = modE*poisson/((1.0+poisson)*(1.0-2.0*poisson))
    
    gamma = rof*g*scale

    tmp1  = la + 2.0*nu

    return (k/gamma)*bMod*(tmp1/(tmp1+bMod*bCoef**2))    
    

    
    