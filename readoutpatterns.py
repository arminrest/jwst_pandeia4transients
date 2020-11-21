#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:18:15 2020

@author: arest
"""
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
   FAST      10     1   27.8     1    27.8
   FAST      10     1   27.8     2    55.5
   FAST      20     1   55.5     2   111.0
   FAST      30     1   83.2     2   166.5
   FAST      40     1  111.0     2   222.0
   FAST      50     1  138.8     2   277.5
   FAST      30     1   83.2     4   333.0
   FAST      40     1  111.0     4   444.0
   FAST      50     1  138.8     4   555.0
   FAST      60     1  166.5     4   666.0
   FAST      70     1  194.2     4   777.0
   FAST      80     1  222.0     4   888.0
   FAST      90     1  249.8     4   999.0
   FAST     100     1  277.5     4  1110.0
   FAST     110     1  305.2     4  1221.0
   FAST     120     1  333.0     4  1332.0
   FAST     130     1  360.8     4  1443.0
   FAST     140     1  388.5     4  1554.0
   FAST     150     1  416.2     4  1665.0
   FAST     160     1  444.0     4  1776.0
   FAST     170     1  471.8     4  1887.0
   FAST     180     1  499.5     4  1998.0
   FAST     190     1  527.2     4  2109.0
   FAST     200     1  555.0     4  2220.0
   FAST     210     1  582.8     4  2331.0
   FAST     220     1  610.5     4  2442.0
   FAST     230     1  638.2     4  2553.0
   FAST     240     1  666.0     4  2664.0
   FAST     250     1  693.8     4  2775.0
   FAST     260     1  721.5     4  2886.0
   FAST     270     1  749.2     4  2997.0
   FAST     280     1  777.0     4  3108.0
   FAST     290     1  804.8     4  3219.0
   FAST     300     1  832.5     4  3330.0
   FAST     310     1  860.2     4  3441.0
   FAST     320     1  888.0     4  3552.0
   FAST     330     1  915.8     4  3663.0
   FAST     340     1  943.5     4  3774.0
   FAST     350     1  971.2     4  3885.0
   FAST     360     1  999.0     4  3996.0
   FAST     190     2  527.2     4  4218.0
   FAST     200     2  555.0     4  4440.0
   FAST     210     2  582.8     4  4662.0
   FAST     220     2  610.5     4  4884.0
   FAST     230     2  638.2     4  5106.0
   FAST     240     2  666.0     4  5328.0
   FAST     250     2  693.8     4  5550.0
   FAST     260     2  721.5     4  5772.0
   FAST     270     2  749.2     4  5994.0
   FAST     280     2  777.0     4  6216.0
   FAST     290     2  804.8     4  6438.0
   FAST     300     2  832.5     4  6660.0
   FAST     310     2  860.2     4  6882.0
   FAST     320     2  888.0     4  7104.0
   FAST     330     2  915.8     4  7326.0
   FAST     340     2  943.5     4  7548.0
   FAST     350     2  971.2     4  7770.0
   FAST     360     2  999.0     4  7992.0
   FAST     250     3  693.8     4  8325.0
   FAST     260     3  721.5     4  8658.0
   FAST     270     3  749.2     4  8991.0
   FAST     280     3  777.0     4  9324.0
   FAST     290     3  804.8     4  9657.0
   FAST     300     3  832.5     4  9990.0
"""

pattern2exptime['nirspec'] ="""
Readout       NGROUP  NINT     tint    NEXP     texp
nrsirs2rapid       4     1      43.7      2    145.9
nrsirs2rapid       4     1      58.4      3    218.8
nrsirs2rapid       8     1     116.7      3    393.9
nrsirs2rapid      12     1     175.1      3    569.0
nrsirs2rapid      16     1     233.4      3    744.0
nrsirs2rapid      20     1     291.8      3    919.1
nrsirs2rapid      24     1     350.1      3   1094.2
nrsirs2rapid      28     1     408.5      3   1269.2
nrsirs2rapid      32     1     466.8      3   1444.3
nrsirs2rapid      36     1     525.2      3   1619.4
nrsirs2rapid      40     1     583.6      3   1794.4
nrsirs2rapid      44     1     641.9      3   1969.5
nrsirs2rapid      48     1     700.3      3   2144.6
nrsirs2rapid      52     1     758.6      3   2319.6
nrsirs2rapid      56     1     817.0      3   2494.7
nrsirs2rapid      60     1     875.3      3   2669.8
nrsirs2rapid      64     1     933.7      3   2844.8
nrsirs2rapid      68     1     992.1      3   3019.9
nrsirs2rapid      72     1    1050.4      3   3195.0
nrsirs2rapid      76     1    1108.8      3   3370.0
nrsirs2rapid      80     1    1167.1      3   3545.1
nrsirs2rapid      84     1    1225.5      3   3720.2
nrsirs2rapid      88     1    1283.8      3   3895.2
nrsirs2rapid      48     2     700.3      3   4289.1
nrsirs2rapid      52     2     758.6      3   4639.3
nrsirs2rapid      56     2     817.0      3   4989.4
nrsirs2rapid      60     2     875.3      3   5339.5
nrsirs2rapid      64     2     933.7      3   5689.7
nrsirs2rapid      68     2     992.1      3   6039.8
nrsirs2rapid      72     2    1050.4      3   6389.9
nrsirs2rapid      76     2    1108.8      3   6740.1
nrsirs2rapid      80     2    1167.1      3   7090.2
nrsirs2rapid      84     2    1225.5      3   7440.3
nrsirs2rapid      88     2    1283.8      3   7790.5
nrsirs2rapid      64     3     933.7      3   8534.5
nrsirs2rapid      68     3     992.1      3   9059.7
nrsirs2rapid      72     3    1050.4      3   9584.9
"""
class readoutpatternclass(pdastroclass):
    def __init__(self,instrument):
        pdastroclass.__init__(self)
        
        
        # time in seconds per group for different readout pattern
        self.tgroup_sec = {}
        # MIRI FAST readout tgroup
        self.tgroup_sec['fast']=2.775
        # NIRSpec nrsirs2rapid readout tgroup
        self.tgroup_sec['nrsirs2rapid']=14.589
        
        self.allowed_instruments = ['nircam','niriss','nirspec','miri']
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

    def calc_t(self,Ngroups,Nint,Nexp,tgroup=None,readout=None,tadd=0.0):
        if tgroup is None:
            if readout is None:
                raise RuntimeError("Neither tgroup nor readout specified, cannot calculate t")    
            tgroup = self.get_tgroup_for_readoutpattern(readout)
        tint = Ngroups*tgroup
        texp = tint*Nint*Nexp+tadd
        return(tint,texp)
    
    def calc_MIRI_exptimes(self,
                           filename=None,
                           # https://jwst-docs.stsci.edu/mid-infrared-instrument/miri-observing-strategies/miri-imaging-recommended-strategies#MIRIImagingRecommendedStrategies-Dwelltimelimit
                           # recommended, but not required is min=40 and max = 360 groups
                           Ng_min=40,Ng_max=360, 
                           Ng_absmin=10,
                           # 4 point dither desired
                           Nexp_max=4,
                           # if Ngroups_modval != None: only keep entries for which 
                           # Ngroups % Ngroups_modval == 0
                           # This reduces the # of entries in the table.
                           Ngroups_modval = 10,
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
        result = self.calc_exptimes('miri','FAST',filename=filename,Ng_min=Ng_min,Ng_max=Ng_max,Ng_absmin=Ng_absmin,
                                    Nexp_max=Nexp_max,Ngroups_modval=Ngroups_modval,tmin=tmin,tmax=tmax)
        return(result)
        
    def calc_NIRSpec_exptimes(self,
                           filename=None,
                           # https://jwst-docs.stsci.edu/near-infrared-spectrograph/nirspec-observing-strategies/nirspec-detector-recommended-strategies#NIRSpecDetectorRecommendedStrategies-groups
                           # recommended, but not required is min=40 and max = 360 groups
                           Ng_min=3,Ng_max=90, 
                           Ng_absmin=10,
                           # 4 point dither desired
                           Nexp_max=4,
                           # if Ngroups_modval != None: only keep entries for which 
                           # Ngroups % Ngroups_modval == 0
                           # This reduces the # of entries in the table.
                           Ngroups_modval = 4,
                           tmin=None,tmax=10000.0,
                           # additive term: for NIRSPEC
                           tadd=58.35):
        """
        same as calc_exptimes() but with good default values for MIRI

        https : //jwst-docs.stsci.edu/mid-infrared-instrument/miri-observing-strategies/miri-imaging-recommended-strategies#MIRIImagingRecommendedStrategies-Dwelltimelimit                           # recommended
        #  recommended, but not required is min = 40 and max = 360 groups
 
        Returns
        -------
        None.
        self.t contains the table

        """
        print('VVV',tadd)
        result = self.calc_exptimes('nirspec','nrsirs2rapid',filename=filename,Ng_min=Ng_min,Ng_max=Ng_max,Ng_absmin=Ng_absmin,
                                    Nexp_max=Nexp_max,Ngroups_modval=Ngroups_modval,tmin=tmin,tmax=tmax,tadd=tadd)
        return(result)

    def calc_exptimes(self,instrument,readout,
                      filename=None,
                      Ng_min=5,Ng_absmin=5,Ng_max=20,
                      Nexp_max=4,Ngroups_modval=None,
                      tmin=None,tmax=10000.0,tadd=0.0):
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
        tadd : float, optional
            For NIRSpec, 58.35sec needs to be added
            The default is 0.0.

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
             
            if self.instrument.lower() == 'nirspec':
                # for nodding, only exposures of 2, 3, and 5 are valid.
                valid = np.array([2,3,5])
                ind = np.argmin(abs(valid-Nexp))
                Nexp = valid[ind]
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
            (tint,texp)=self.calc_t(Ngroups,Nint,Nexp,tgroup=tgroup,tadd=tadd*Nint)
            
            # check for limits
            if not(tmin is None) and (texp<tmin):useflag=False
            if not(tmax is None) and (texp>tmax):useflag=False

            # add the line to the string
            if useflag: pattern2exptime[self.instrument]+='%s\t%d\t%d\t%.1f\t%d\t%.1f\n' % (readout,Ngroups,Nint,tint,Nexp,texp)

            Ngroups_tot+=(Nexp*Nint)
        
        # parse the string into pdastro
        self.t = pd.read_csv(io.StringIO(pattern2exptime[self.instrument]),delim_whitespace=True,skipinitialspace=True)

        if not(filename is None):
            print('Writing table to %s' % filename)
            self.write(filename)
        return(0)                                  

    def loadreadoutpatterntable(self,instrument):
        self.set_instrument(instrument)
        if self.instrument in ['nircam','miri','nirspec']:     
            self.t = pd.read_csv(io.StringIO(pattern2exptime[self.instrument]),delim_whitespace=True,skipinitialspace=True)
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
