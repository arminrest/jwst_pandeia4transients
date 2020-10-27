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
pattern2exptime['miri'] = """
Readout  NGROUP  NINT   tint  NEXP    texp
   FAST      12     1   33.3     1    33.3
   FAST      16     1   44.4     1    44.4
   FAST      12     1   33.3     2    66.6
   FAST      16     1   44.4     2    88.8
   FAST      20     1   55.5     2   111.0
   FAST      24     1   66.6     2   133.2
   FAST      28     1   77.7     2   155.4
   FAST      32     1   88.8     2   177.6
   FAST      36     1   99.9     2   199.8
   FAST      40     1  111.0     2   222.0
   FAST      44     1  122.1     2   244.2
   FAST      48     1  133.2     2   266.4
   FAST      52     1  144.3     2   288.6
   FAST      56     1  155.4     2   310.8
   FAST      32     1   88.8     4   355.2
   FAST      36     1   99.9     4   399.6
   FAST      40     1  111.0     4   444.0
   FAST      44     1  122.1     4   488.4
   FAST      48     1  133.2     4   532.8
   FAST      52     1  144.3     4   577.2
   FAST      56     1  155.4     4   621.6
   FAST      60     1  166.5     4   666.0
   FAST      64     1  177.6     4   710.4
   FAST      68     1  188.7     4   754.8
   FAST      72     1  199.8     4   799.2
   FAST      76     1  210.9     4   843.6
   FAST      80     1  222.0     4   888.0
   FAST      84     1  233.1     4   932.4
   FAST      88     1  244.2     4   976.8
   FAST      92     1  255.3     4  1021.2
   FAST      96     1  266.4     4  1065.6
   FAST     100     1  277.5     4  1110.0
   FAST     104     1  288.6     4  1154.4
   FAST     108     1  299.7     4  1198.8
   FAST     112     1  310.8     4  1243.2
   FAST     116     1  321.9     4  1287.6
   FAST     120     1  333.0     4  1332.0
   FAST     124     1  344.1     4  1376.4
   FAST     128     1  355.2     4  1420.8
   FAST     132     1  366.3     4  1465.2
   FAST     136     1  377.4     4  1509.6
   FAST     140     1  388.5     4  1554.0
   FAST     144     1  399.6     4  1598.4
   FAST     148     1  410.7     4  1642.8
   FAST     152     1  421.8     4  1687.2
   FAST     156     1  432.9     4  1731.6
   FAST     160     1  444.0     4  1776.0
   FAST     164     1  455.1     4  1820.4
   FAST     168     1  466.2     4  1864.8
   FAST     172     1  477.3     4  1909.2
   FAST     176     1  488.4     4  1953.6
   FAST     180     1  499.5     4  1998.0
   FAST     184     1  510.6     4  2042.4
   FAST     188     1  521.7     4  2086.8
   FAST     192     1  532.8     4  2131.2
   FAST     196     1  543.9     4  2175.6
   FAST     200     1  555.0     4  2220.0
   FAST     204     1  566.1     4  2264.4
   FAST     208     1  577.2     4  2308.8
   FAST     212     1  588.3     4  2353.2
   FAST     216     1  599.4     4  2397.6
   FAST     220     1  610.5     4  2442.0
   FAST     224     1  621.6     4  2486.4
   FAST     228     1  632.7     4  2530.8
   FAST     232     1  643.8     4  2575.2
   FAST     236     1  654.9     4  2619.6
   FAST     240     1  666.0     4  2664.0
   FAST     244     1  677.1     4  2708.4
   FAST     248     1  688.2     4  2752.8
   FAST     252     1  699.3     4  2797.2
   FAST     256     1  710.4     4  2841.6
   FAST     260     1  721.5     4  2886.0
   FAST     264     1  732.6     4  2930.4
   FAST     268     1  743.7     4  2974.8
   FAST     272     1  754.8     4  3019.2
   FAST     276     1  765.9     4  3063.6
   FAST     280     1  777.0     4  3108.0
   FAST     284     1  788.1     4  3152.4
   FAST     288     1  799.2     4  3196.8
   FAST     292     1  810.3     4  3241.2
   FAST     296     1  821.4     4  3285.6
   FAST     300     1  832.5     4  3330.0
   FAST     304     1  843.6     4  3374.4
   FAST     308     1  854.7     4  3418.8
   FAST     312     1  865.8     4  3463.2
   FAST     316     1  876.9     4  3507.6
   FAST     320     1  888.0     4  3552.0
   FAST     324     1  899.1     4  3596.4
   FAST     328     1  910.2     4  3640.8
   FAST     332     1  921.3     4  3685.2
   FAST     336     1  932.4     4  3729.6
   FAST     340     1  943.5     4  3774.0
   FAST     344     1  954.6     4  3818.4
   FAST     348     1  965.7     4  3862.8
   FAST     352     1  976.8     4  3907.2
   FAST     356     1  987.9     4  3951.6
   FAST     360     1  999.0     4  3996.0
   FAST     184     2  510.6     4  4084.8
   FAST     188     2  521.7     4  4173.6
   FAST     192     2  532.8     4  4262.4
   FAST     196     2  543.9     4  4351.2
   FAST     200     2  555.0     4  4440.0
   FAST     204     2  566.1     4  4528.8
   FAST     208     2  577.2     4  4617.6
   FAST     212     2  588.3     4  4706.4
   FAST     216     2  599.4     4  4795.2
   FAST     220     2  610.5     4  4884.0
   FAST     224     2  621.6     4  4972.8
   FAST     228     2  632.7     4  5061.6
   FAST     232     2  643.8     4  5150.4
   FAST     236     2  654.9     4  5239.2
   FAST     240     2  666.0     4  5328.0
   FAST     244     2  677.1     4  5416.8
   FAST     248     2  688.2     4  5505.6
   FAST     252     2  699.3     4  5594.4
   FAST     256     2  710.4     4  5683.2
   FAST     260     2  721.5     4  5772.0
   FAST     264     2  732.6     4  5860.8
   FAST     268     2  743.7     4  5949.6
   FAST     272     2  754.8     4  6038.4
   FAST     276     2  765.9     4  6127.2
   FAST     280     2  777.0     4  6216.0
   FAST     284     2  788.1     4  6304.8
   FAST     288     2  799.2     4  6393.6
   FAST     292     2  810.3     4  6482.4
   FAST     296     2  821.4     4  6571.2
   FAST     300     2  832.5     4  6660.0
   FAST     304     2  843.6     4  6748.8
   FAST     308     2  854.7     4  6837.6
   FAST     312     2  865.8     4  6926.4
   FAST     316     2  876.9     4  7015.2
   FAST     320     2  888.0     4  7104.0
   FAST     324     2  899.1     4  7192.8
   FAST     328     2  910.2     4  7281.6
   FAST     332     2  921.3     4  7370.4
   FAST     336     2  932.4     4  7459.2
   FAST     340     2  943.5     4  7548.0
   FAST     344     2  954.6     4  7636.8
   FAST     348     2  965.7     4  7725.6
   FAST     352     2  976.8     4  7814.4
   FAST     356     2  987.9     4  7903.2
   FAST     360     2  999.0     4  7992.0
   FAST     244     3  677.1     4  8125.2
   FAST     248     3  688.2     4  8258.4
   FAST     252     3  699.3     4  8391.6
   FAST     256     3  710.4     4  8524.8
   FAST     260     3  721.5     4  8658.0
   FAST     264     3  732.6     4  8791.2
   FAST     268     3  743.7     4  8924.4
   FAST     272     3  754.8     4  9057.6
   FAST     276     3  765.9     4  9190.8
   FAST     280     3  777.0     4  9324.0
   FAST     284     3  788.1     4  9457.2
   FAST     288     3  799.2     4  9590.4
   FAST     292     3  810.3     4  9723.6
   FAST     296     3  821.4     4  9856.8
   FAST     300     3  832.5     4  9990.0
"""

class readoutpatternclass(pdastroclass):
    def __init__(self,instrument):
        pdastroclass.__init__(self)
        
        
        # time in seconds per group for different readout pattern
        self.tgroup_sec = {}
        # MIRI FAST readout tgroup
        self.tgroup_sec['fast']=2.775
        
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
    
    def get_tgroup_for_readoutpattern(self,readout):
        readout=readout.lower()
        if not (readout in self.tgroup_sec):
            raise RuntimeError('readout %s is not yet supported, update self.tgroup_sec!' % readout)
        return(self.tgroup_sec[readout])

    def calc_t(self,Ngroups,Nint,Nexp,tgroup=None,readout=None):
        if tgroup is None:
            if readout is None:
                raise RuntimeError("Neither tgroup nor readout specified, cannot calculate t")    
            tgroup = self.get_tgroup_for_readoutpattern(readout)
        tint = Ngroups*tgroup
        texp = tint*Nint*Nexp
        return(tint,texp)
    
    def calc_MIRI_exptimes(self,
                           # https://jwst-docs.stsci.edu/mid-infrared-instrument/miri-observing-strategies/miri-imaging-recommended-strategies#MIRIImagingRecommendedStrategies-Dwelltimelimit
                           # recommended, but not required is min=40 and max = 360 groups
                           Ng_min=40,Ng_max=360, 
                           Ng_absmin=10,
                           # 4 point dither desired
                           Nexp_max=4,
                           # if Ngroups_modval != None: only keep entries for which 
                           # Ngroups % Ngroups_modval == 0
                           # This reduces the # of entries in the table.
                           Ngroups_modval = 4,
                           tmin=None,tmax=10000.0):
        """
        same as calc_exptimes() but with good default values for MIRI

        https : //jwst-docs.stsci.edu/mid-infrared-instrument/miri-observing-strategies/miri-imaging-recommended-strategies#MIRIImagingRecommendedStrategies-Dwelltimelimit                           # recommended
        #  recommended, but not required is min = 40 and max = 360 groups
 
        Returns
        -------
        None.
        self.t contains the table

        """
        result = self.calc_exptimes('miri','FAST',Ng_min=Ng_min,Ng_max=Ng_max,Ng_absmin=Ng_absmin,
                                    Nexp_max=Nexp_max,Ngroups_modval=Ngroups_modval,tmin=tmin,tmax=tmax)
        return(result)
        

    def calc_exptimes(self,instrument,readout,
                      filename=None,
                      Ng_min=5,Ng_absmin=5,Ng_max=20,
                      Nexp_max=4,Ngroups_modval=None,
                      tmin=None,tmax=10000.0):
        """
        Calculate a list of exposure times with reasonable readout pattern

        Parameters
        ----------
        instrument : string
            one of NIRCam, MIRI, NIRSpec, NIRISS (case-insensitive)
        readout : string
            readout pattern names, e.g., 'FAST' for MIRI, 'SHALLOW4' for NIRCam.
        filename : string, optional
            Save table to this filename if not None. The default is None.
        Ng_min : int, optional
            Minimum desired number of groups pre integration 
            (recommended but not required). The default is 5.
        Ng_absmin : int, optional
            absolute minimum number of groups, in case of very short exposures 
            with Ngroups<self.Ng_min. The default is 5.
        Ng_max : int, optional
            Maximum number of groups per integration. The default is 20.
        Nexp_max : int, optional
            Maximum number of exposures. The default is 4.
        Ngroups_modval : int, optional
            if Ngroups_modval != None: only keep entries for which 
            Ngroups % Ngroups_modval == 0
            This reduces the # of entries in the table. The default is None.
        tmin : float, optional
            minimum total exposure times. The default is None.
        tmax : float, optional
            maximum total exposure time. The default is 10000.0.

        Returns
        -------
        None.
        self.t contains the table
        """
        self.set_instrument(instrument)
        # get the time per group
        tgroup = self.get_tgroup_for_readoutpattern(readout)
        
        # header of table
        pattern2exptime[self.instrument]='Readout\tNGROUP\tNINT\ttint\tNEXP\ttexp\n'

        #initialize
        Ngroups_tot = Ng_absmin        
        Nint=1

        # set to zero for 'while' condition
        texp=0
        while texp<tmax:
            # how many exposures assuming desired minimum # of groups?
            Nexp = int(Ngroups_tot/Ng_min)
            # If the # of exposures is 1 or smaller, then use self.MIRI_Ng_absmin
            if Nexp<2:
                Nexp = int(Ngroups_tot/Ng_absmin)
                # We want to stick with 2 exposures, no reason to go to more 
                # exposures with such a small # of groups per integration
                if Nexp>2: Nexp=2
            else:
                # avoid Nexp=3, bad dither pattern
                if Nexp==3: Nexp=4
            
            # If it's more than self.MIRI_max_Nexp exposures, stick with self.MIRI_max_Nexp
            if Nexp>Nexp_max:
                Nexp = Nexp_max
               
            # get the number of groups per integration
            Ngroups = int(Ngroups_tot/(Nint*Nexp))

            if Ngroups>Ng_max:
                Nint+=1
                Ngroups = int(Ngroups_tot/(Nint*Nexp))
            
            # only go for integer steps!
            if Ngroups*Nint*Nexp!=Ngroups_tot:
                Ngroups_tot+=1
                continue
            
            # as default, use this entry
            useflag=True

            # check for 'mod' if specified.
            if not(Ngroups_modval is None):
                if (Ngroups % Ngroups_modval) != 0:
                    useflag=False
                    
            #calculate integration and exposure times
            (tint,texp)=self.calc_t(Ngroups,Nint,Nexp,tgroup=tgroup)
            
            # check for limits
            if not(tmin is None) and (texp<tmin):useflag=False
            if not(tmax is None) and (texp>tmax):useflag=False

            # add the line to the string
            if useflag: pattern2exptime[self.instrument]+='%s\t%d\t%d\t%.1f\t%d\t%.1f\n' % (readout,Ngroups,Nint,tint,Nexp,texp)

            Ngroups_tot+=(Nexp*Nint)
        
        # parse the string into pdastro
        self.t = pd.read_csv(io.StringIO(pattern2exptime[self.instrument]),delim_whitespace=True,skipinitialspace=True)
        self.write()
        return(0)                                  


    # OLD!!! DELETE!!!!
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
        if self.instrument in ['nircam','miri']:     
            self.t = pd.read_csv(io.StringIO(pattern2exptime[self.instrument]),delim_whitespace=True,skipinitialspace=True)
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
        if self.instrument in ['nircam','nirspec','miri']:
            info = {'readout_pattern':self.t.at[index,'Readout'].lower(),
                    'NGROUP':self.t.at[index,'NGROUP'],
                    'NINT':self.t.at[index,'NINT'],
                    'tint':self.t.at[index,'tint'],
                    'NEXP':self.t.at[index,'NEXP'],
                    'texp':self.t.at[index,'texp']}
        else:
            raise RuntimeError('instrument %s not yet implemented!' % self.instrument)
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
    
    print(readoutpattern.t)
    index = readoutpattern.index4closestexptime(20000)
    print('bbbbaa1',readoutpattern.t.at[index,'texp'])
    index = readoutpattern.index4nextbiggestexptime(2200)
    print('bbbbaa2',readoutpattern.t.at[index,'texp'])
    info = readoutpattern.info4nextbiggestexptime(2200)
    print('bbbbaa3',info)
    info = readoutpattern.info4closestexptime(880)
    print('bbbbaa4',info)
