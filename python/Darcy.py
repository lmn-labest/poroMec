#!/bin/python
import sys  

def main():
  
  nArgs=len(sys.argv)

  print (nArgs)
  if nArgs < 2:
      print (sys.argv[0],'[mDarcy]')
      exit(0)

  roAgua = 1.e+3
  g      = 9.81
  muAgua = 1.003e-3
    
# mili darcy
  mD      = float(sys.argv[1])

# permeabilidae intrisica
  ki = mD*9.869233e-16

# permeabilidae hidraulicaa   
  k = ki*roAgua*g/muAgua 

  print ('Ki = %e'%ki)
  print ('k  = %e'%k)

# ...
if __name__ == "__main__":
  main()
# .....................................................................
