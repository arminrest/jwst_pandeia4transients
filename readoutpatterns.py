#!/usr/bin/env python
import numpy as np
import math
import sys,socket,os,re
import pandas as pd
from pdastro import pdastroclass, AandB
import io


pattern2exptime = {}
pattern2exptime['nircam'] = """
Readout         NGROUPS NINT    tint    NEXP    texp
RAPID           2       1       21.5    4       86.0
BRIGHT2         2       1       42.9    4       171.8
BRIGHT2         3       1       64.4    4       257.6
BRIGHT2         4       1       85.9    4       343.6
BRIGHT1         5       1       96.6    4       386.4
BRIGHT1         6       1       118.1   4       472.4
BRIGHT1         7       1       139.6   4       558.4
BRIGHT1         8       1       161.1   4       644.4
BRIGHT1         9       1       182.5   4       730.0
BRIGHT1         10      1       204.0   4       816.0
SHALLOW4        5       1       257.7   4       1030.8
SHALLOW4        6       1       311.4   4       1245.6
SHALLOW4        7       1       365.1   4       1460.4
SHALLOW4        8       1       418.7   4       1674.8
SHALLOW4        9       1       472.4   4       1889.6
SHALLOW4        10      1       526.1   4       2104.4
MEDIUM8         6       1       622.7   4       2490.8
MEDIUM8         7       1       730.1   4       2920.4
MEDIUM8         8       1       837.5   4       3350.0
MEDIUM8         9       1       944.8   4       3779.2
MEDIUM8         10      1       1052.2  4       4208.8
"""

class readoutpatternclass(pdastroclass):
    def __init__(self,instrument):
        pdastroclass.__init__(self)
        self.allowed_instruments = ['nircam','nirspec','miri']
        self.instrument=instrument.lower()
        if not (self.instrument in self.allowed_instruments):
            raise(RuntimeError,'instrument %s not in %s' % (self.instrument,' '.join(self.allowed_instruments)))    
        self.t = pd.read_csv(io.StringIO(pattern2exptime[self.instrument]),delim_whitespace=True,skipinitialspace=True)

    def getinfo(self,index):
        if index == None:
            return(None)
        if self.instrument == 'nircam':
            info = {'Readout':self.t.at[index,'Readout'],
                    'NGROUPS':self.t.at[index,'NGROUPS'],
                    'NINT':self.t.at[index,'NINT'],
                    'tint':self.t.at[index,'tint'],
                    'NEXP':self.t.at[index,'NEXP'],
                    'texp':self.t.at[index,'texp']}
        else:
            raise (RuntimeError,'instrument %s not yet implemented!' % self.instrument)
        return(info)
                    

    def index4closestexptime(self,exptime):
        indices = self.ix_inrange('texp',exptime)
    #        self.write(indices=indices)
        if len(indices)<1:
            return(self.t.index.size-1)
        elif (indices[0]==0):
            return(0)
        else:
            dtminus=exptime-self.t.at[indices[0]-1,'texp']
            dtplus=self.t.at[indices[0],'texp']-exptime
            if dtplus<dtminus:
                index=indices[0]
            else:
                index=indices[0]-1
            return(index)

    def index4nextbiggestexptime(self,exptime):
        indices = self.ix_inrange('texp',exptime)
        if len(indices)<1:
            return(None)
        else:
            return(indices[0])


    def info4nextbiggestexptime(self,exptime):
        index = self.index4nextbiggestexptime(exptime)
        return(self.getinfo(index))

    def info4closestexptime(self,exptime):
        index = self.index4closestexptime(exptime)
        return(self.getinfo(index))
        

if __name__ == '__main__':
    readoutpattern=readoutpatternclass('nircam')
    
    print(readoutpattern.t)
    index = readoutpattern.index4closestexptime(20000)
    print('bbbbaa1',readoutpattern.t.at[index,'texp'])
    index = readoutpattern.index4nextbiggestexptime(2200)
    print('bbbbaa2',readoutpattern.t.at[index,'texp'])
    info = readoutpattern.info4nextbiggestexptime(2200)
    print('bbbbaa3',info)
    info = readoutpattern.info4closestexptime(880)
    print('bbbbaa4',info)
