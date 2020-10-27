#!/usr/bin/env python3
import numpy as np
import math
import sys,socket,os,re
import pandas as pd
from pdastro import pdastroclass, AandB
import io


pattern2exptime = {}
pattern2exptime['nircam'] = """
Readout         NGROUP  NINT    tint    NEXP    texp
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
        
        # MIRI: group limits 
        # https://jwst-docs.stsci.edu/mid-infrared-instrument/miri-observing-strategies/miri-imaging-recommended-strategies#MIRIImagingRecommendedStrategies-Dwelltimelimit
        # recommended, but not required is min=40 and max = 360 groups
        self.MIRI_Ng_min = 40
        self.MIRI_Ng_max = 360
        # absolute minimum number of groups, in case of very short exposures with Ngroups<self.MIRI_Ng_min
        self.MIRI_Ng_absmin = 10
        # time in seconds per group in MIRI images
        self.MIRI_tgroup_sec = {}
        self.MIRI_tgroup_sec['fast']=2.775
        
        # if MIRI_Ngroups_modval != None: only keep entries for which 
        # Ngroups % MIRI_Ngroups_modval == 0
        # This reduces the # of entries in the table.
        self.MIRI_Ngroups_modval = 4
        
        # maximum # of exposures
        self.MIRI_max_Nexp=4
        
        # to what total exposure time should 
        self.MIRI_max_tint = 10000.0
        
        self.allowed_instruments = ['nircam','niriss','nirspec','miri']
        self.set_instrument(instrument)
        self.loadreadoutpatterntable(instrument)

    def set_instrument(self,instrument):
        if not(instrument is None):
            self.instrument = instrument.lower()
        if self.instrument is None:
            raise(RuntimeError,'An instrument needs to be specified!!')    
        if not (self.instrument in self.allowed_instruments):
            raise(RuntimeError,'instrument %s not in %s' % (self.instrument,' '.join(self.allowed_instruments))) 
        
        return(0)

    def calc_t_MIRI(self,readout,Ngroups,Nint,Nexp):
        tgroup = self.MIRI_tgroup_sec[readout.lower()]
        tint = Ngroups*tgroup*Nint
        texp = tint*Nexp
        return(tint,texp)
    
    def calcMIRIexptimes(self,readout='FAST',tmin=None,tmax=10000.0):
        # header of table
        pattern2exptime['miri']='Readout\tNGROUP\tNINT\ttint\tNEXP\ttexp\n'

        #initialize
        Ngroups_tot = self.MIRI_Ng_absmin        
        Nint=1

        # set to zero for 'while' condition
        texp=0
        while texp<tmax:
            # how many exposures assuming desired minimum # of groups?
            Nexp = int(Ngroups_tot/self.MIRI_Ng_min)
            # If the # of exposures is 1 or smaller, then use self.MIRI_Ng_absmin
            if Nexp<2:
                Nexp = int(Ngroups_tot/self.MIRI_Ng_absmin)
                # We want to stick with 2 exposures, no reason to go to more 
                # exposures with such a small # of groups per integration
                if Nexp>2: Nexp=2
            else:
                # avoid Nexp=3, bad dither pattern
                if Nexp==3: Nexp=4
            
            # If it's more than self.MIRI_max_Nexp exposures, stick with self.MIRI_max_Nexp
            if Nexp>self.MIRI_max_Nexp:
                Nexp = self.MIRI_max_Nexp
               
            # get the number of groups per integration
            Ngroups = int(Ngroups_tot/(Nint*Nexp))

            if Ngroups>self.MIRI_Ng_max:
                Nint+=1
                Ngroups = int(Ngroups_tot/(Nint*Nexp))
            
            # only go for integer steps!
            if Ngroups*Nint*Nexp!=Ngroups_tot:
                Ngroups_tot+=1
                continue
            
            # as default, use this entry
            useflag=True

            # check for 'mod' if specified.
            if not(self.MIRI_Ngroups_modval is None):
                if (Ngroups % self. MIRI_Ngroups_modval) != 0:
                    useflag=False
                    
            #calculate integration and exposure times
            (tint,texp)=self.calc_t_MIRI(readout,Ngroups,Nint,Nexp)
            
            # check for limits
            if not(tmin is None) and (texp<tmin):useflag=False
            if not(tmax is None) and (texp>tmax):useflag=False

            # add the line to the string
            if useflag: pattern2exptime['miri']+='%s\t%d\t%d\t%.1f\t%d\t%.1f\n' % (readout,Ngroups,Nint,tint,Nexp,texp)

            Ngroups_tot+=(Nexp*Nint)
        
        # parse the string into pdastro
        self.t = pd.read_csv(io.StringIO(pattern2exptime[self.instrument]),delim_whitespace=True,skipinitialspace=True)
        self.write()
                                  
        
        #while texp < self.MIRI_max_tint 

    def loadreadoutpatterntable(self,instrument):
        if self.instrument=='nircam':     
            self.t = pd.read_csv(io.StringIO(pattern2exptime[self.instrument]),delim_whitespace=True,skipinitialspace=True)
        elif self.instrument=='miri':
            self.t = self.calcMIRIexptimes()
        elif self.instrument == 'nirspec':
            path = os.path.dirname(os.path.abspath(__file__))
            file = os.path.join(path,'data/nirspecpattern2exptime.csv')
            self.t = pd.read_csv(path)
        else:
            raise RuntimeError("instrument %s not yet implemented!" % self.instrument)
        return(0)

    def getinfo(self,index):
        if index == None:
            return(None)
        if (self.instrument == 'nircam') | (self.instrument == 'nirspec'):
            info = {'readout_pattern':self.t.at[index,'Readout'].lower(),
                    'NGROUP':self.t.at[index,'NGROUP'],
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

    def index4nextsmallesexptime(self,exptime):
        indices = self.ix_inrange('texp',None,exptime)
        if len(indices)<1:
            return(None)
        else:
            return(indices[-1])

    def info4nextbiggestexptime(self,exptime):
        index = self.index4nextbiggestexptime(exptime)
        return(self.getinfo(index))

    def info4nextsmallesexptime(self,exptime):
        index = self.index4nextsmallesexptime(exptime)
        return(self.getinfo(index))

    def info4closestexptime(self,exptime):
        index = self.index4closestexptime(exptime)
        return(self.getinfo(index))
        
    def nextbiggestexptime(self,exptime):
        index = self.index4nextbiggestexptime(exptime)
        if index is None: 
            return(None)
        else:
            return(self.t.at[index,'texp'])

    def nextsmallestexptime(self,exptime):
        index = self.index4nextsmallesexptime(exptime)
        if index is None: 
            return(None)
        else:
            return(self.t.at[index,'texp'])

    def closestexptime(self,exptime):
        index = self.index4closestexptime(exptime)
        if index is None: 
            return(None)
        else:
            return(self.t.at[index,'texp'])
        

if __name__ == '__main__':
    readoutpattern=readoutpatternclass('miri')
    
    sys.exit(0)
    
    print(readoutpattern.t)
    index = readoutpattern.index4closestexptime(20000)
    print('bbbbaa1',readoutpattern.t.at[index,'texp'])
    index = readoutpattern.index4nextbiggestexptime(2200)
    print('bbbbaa2',readoutpattern.t.at[index,'texp'])
    info = readoutpattern.info4nextbiggestexptime(2200)
    print('bbbbaa3',info)
    info = readoutpattern.info4closestexptime(880)
    print('bbbbaa4',info)
