import sys
import pandas as pd
import numpy as np

def giveXplorAcoFromCyana(cyanaacofile):
   '''

   '''

   data = pd.read_csv(cyanaacofile,sep=r"\s+",header=None)

   for i in range(len(data)):
        resi     = int(data.iloc[i,0])
        phiorpsi = data.iloc[i,2]
        #acohi    = float(data.iloc[i,3])  corrected on 23-02-21
        #acolo    = float(data.iloc[i,4])  corrected on 23-02-21

        acolo    = float(data.iloc[i,3])
        acohi    = float(data.iloc[i,4])

        midaco    = (acohi + acolo)/2
        rangeaco  = acohi - midaco

        if phiorpsi.lower() == 'phi':
            print(' ASSIGN  (resid %-4d and name C   ) (resid %-4d and name N   )'%(resi-1,resi))
            print('         (resid %-4d and name CA  ) (resid %-4d and name C   )  1.00 %-5.2f %-5.2f 2\n'%(resi,resi,midaco,rangeaco))
        elif phiorpsi.lower() == 'psi':
            print(' ASSIGN  (resid %-4d and name N   ) (resid %-4d and name CA  )'%(resi,resi))
            print('         (resid %-4d and name C   ) (resid %-4d and name N   )  1.00 %-5.2f %-5.2f 2\n'%(resi,resi+1,midaco,rangeaco))
        elif phiorpsi.lower() == 'chi1':
            print(' ASSIGN  (resid %-4d and name N   ) (resid %-4d and name CA  )'%(resi,resi))
            print('         (resid %-4d and name CB  ) (resid %-4d and name CG  )  1.00 %-5.2f %-5.2f 2\n'%(resi,resi,midaco,rangeaco))
        else:
            print('UNKNOWN: mark %s'%(phiorpsi))



if __name__ == "__main__":
   giveXplorAcoFromCyana(sys.argv[1])