#!/bin/python
import sys  
def main():
  
  nArgs=len(sys.argv)

  if nArgs < 6:
    print sys.argv[0],'[Mod de elasticidade] [coef. poison] [force] [mod Biot] [coef Biot]'
    exit(0)

  E       = float(sys.argv[1])
  ni      = float(sys.argv[2])
  force   = float(sys.argv[3])
  modBiot = float(sys.argv[4])
  coefBiot= float(sys.argv[5])
  
  Mi     = 0.5*(E/(1.0+ni)) 
  Lambda = (E*ni)/((1.0+ni)*(1.0-2.0*ni))

  print 'mi    ',Mi
  print 'lambda',Lambda

  print 'z ',force/(Lambda+2.0*Mi)
  print 'dp',(1.e-6)*(modBiot*coefBiot*force)/(Lambda+2.0*Mi)

# ...
if __name__ == "__main__":
  main()
# .....................................................................
