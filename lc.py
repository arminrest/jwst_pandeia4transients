#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:18:15 2020

@author: arest
"""

import math,sys,socket,os,re
if 'TOOLS_PANDEIA_SOURCEDIR' in os.environ:
    print('NNN')
    sys.path.append(os.environ['TOOLS_PANDEIA_SOURCEDIR'])

import numpy as np
import pandas as pd
from pdastro import pdastroclass
import io
import argparse

class lcclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)
        self.lcrootdir = None
        self.lcmodel = None
        self.modeldir = None
        self.verbose=0
        
        self.NIRCamfilters = ['F070W','F090W','F115W','F150W','F150W2','F200W','F277W','F322W2','F356W','F444W']
        self.MIRIfilters = ['F560W','F770W','F1000W','F1130W','F1280W','F1500W','F1800W','F2100W','F2550W']
        
        # subdirs for the different models. Makes it easier
        self.subdirs4models = {}
        self.subdirs4models['metzger']='metzger_20201007'
        
        self.imaging_exptime = pdastroclass()

        self.reftime = None
        self.reffilter = None
        self.refmag = None
        self.dt_discovery = None
        self.dt_trigger = None

        self.lc_maxmag = 30.0
        
        self.phasecolname = 'time'
        
    def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
        
        lcrootdir = None
        if 'LC_DATA' in os.environ:
            lcrootdir = os.environ['LC_DATA']
        
        parser.add_argument("reftime", type=float, help=("time difference between KN and to GW discovery in days"))
        parser.add_argument("reffilter",  help=("reference filter"))
        parser.add_argument("refmag", type=float, help=("reference mag in reference filter"))
                            
        parser.add_argument('-v','--verbose', default=0, action='count')
        parser.add_argument('--lcrootdir',  type=str, default=lcrootdir, help=('specify the rootdir for the light curves (default=%(default)s)'))
        parser.add_argument('--lcmodel',  type=str, default='metzger', help=('specify the model for the light curve. supersedes lcsubdir (default=%(default)s)'))
        parser.add_argument('--lcparams',nargs='+', default=[0.058653,0.143448,0.300000], help=('specify the parameters for the model light curves (default=%(default)s)'))
        parser.add_argument('--lcfilename', default=None, help=('specify the lc filename, overrides --lcmodel and --lcparams (default=%(default)s)'))
        #parser.add_argument('--lcsubdir',  type=str, default=None, help=('specify the subdir for the light curves. Is superseded by lcmodel (default=%(default)s)'))
        
        return(parser)
    
    def getlcdir(self,lcrootdir,lcmodel):
        if lcmodel in self.subdirs4models:
            lcsubdir = self.subdirs4models[lcmodel]
        else:
            lcsubdir = lcmodel
        lcdir="%s/%s" % (lcrootdir,lcsubdir)
        if self.verbose>1: print('lc dir: %s' % lcdir)
        return(lcdir)
        
    def findlcfilename(self,lcrootdir,lcmodel,lcparams):
        lcdir = self.getlcdir(lcrootdir,lcmodel)
        if lcmodel == 'metzger':
            if len(lcparams)!=3: raise RuntimeError("for metzger models: lc params should have 3 entries: vej, mej, and Ye")
            (mej, vej, Ye)=lcparams
            filename = '%s/Metzger_%.6f_%.6f_%.6f.dat' % (lcdir,float(mej),float(vej),float(Ye))
            self.phasecolname = 'time'
            if self.verbose>1: print('lc filename: %s' % filename)
        else:
            raise RuntimeError('lcmodel %s not yet implemented, cannot find filename' % lcmodel)
        return(filename)
    
    def loadmodellc(self,lcfilename,jwstfilter2caps=True):
        if self.verbose: print('Loading lc: %s' % lcfilename)
        self.load_spacesep(lcfilename,comment='#')

        # rename JWST filters to capital letters
        if jwstfilter2caps:
            namesMapping={}
            for col in self.t.columns:
                print(col)
                if re.search('^f\d+',col):
                    namesMapping[col]=col.upper()
            if len(namesMapping)>0:
                print(namesMapping)
                self.t = self.t.rename(columns=namesMapping)
        
    def get_normalized_lc(self,reftime,reffilter,refmag,phaserange,filters,maxmag=None):
        self.initspline(self.phasecolname,reffilter)
        m = self.getspline(reftime,reffilter)
        offset = refmag-m
        if self.verbose: print('%s=%.2f at %f: thus offset=%.2f=%.2f - %.2f' % (reffilter,m,reftime,offset,refmag,m)) 

        cols=[self.phasecolname]
        cols.extend(filters)
        lcnorm = pdastroclass(columns=cols)
        lcnorm.t[self.phasecolname]=phaserange

        #Make sure filters is a list and not a string
        if isinstance(filters,str):filters=[filters]

        for filt in filters:
            self.initspline(self.phasecolname,filt)
            mags = [self.getspline(phase,filt)+offset for phase in phaserange]
            lcnorm.t[filt]=mags
            
        lcnorm.write()
        if not (maxmag is None):
            for filt in filters:
                lcnorm.t.loc[lcnorm.ix_inrange(filt,maxmag),filt]=np.nan
            
            
        lcnorm.write()
        return(lcnorm)
        
if __name__ == '__main__':
    lc=lcclass()
    parser = lc.define_args()
    args = parser.parse_args()

    lc.verbose=args.verbose
    filename = lc.findlcfilename(args.lcrootdir,args.lcmodel,args.lcparams)
    lc.loadmodellc(filename)
    
    indices = lc.ix_remove_null(colnames=['time','f070w'])
    lc.initspline('time','f070w',indices = indices)
    lc.write(indices=indices)
    for x in range(0,20):
        print(x,lc.getspline(x,'f070w'))