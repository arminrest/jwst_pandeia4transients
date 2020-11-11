#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:15:53 2020

@author: arest
"""

import argparse,sys
import numpy as np
from jwst_SNR import jwst_SNRclass
from lc import lcclass

class lc2texpclass(jwst_SNRclass):
    def __init__(self,**kwargs):
        jwst_SNRclass.__init__(self,**kwargs)
        #lcclass.__init__(self)
        self.lc = lcclass()
        self.lcfilters = []
        self.phasecol = 'time'
        
    def lc2texp_table(self, filters, SNR, lambkg4ETC=None,texp_type='best',
                      SNR_tolerance_in_percent=10,saveSNRflag=False,
                      texpformat=None,SNRformat=None,magformat=None,
                      spec=None,**kwargs):
        """

        Parameters
        ----------
        filters : list or string
            Pass (a list of) filter(s) for which to calculate the SNR.
        phaserange : list/tuple of magnitudes
            Pass a list of phases for which to calculate the SNR..
        SNR : float 
            S/N  for which the exposure time should be calculated for.
        lambkg4ETC : list/tuple of two arrays, optional
            the first array is the wavelength, the second the background. 
            The default is None. If none, then the background data in 
            self.lambkg4ETC is used.
        texp_type : string, optional
            'minus', 'plus', or 'best'. The default is 'best'. This is used
            to select which of the exposure times from self.texp4SNRatmag is 
            used
        SNR_tolerance_in_percent : float, optional
            An exposure time with a SNRminus<SNR can be accepted as 'best' if
            SNRminus is closer to SNR than SNRplus (or SNRplus==None), 
            AND if (SNR0-SNRminus)/SNR0<SNR_tolerance_in_percent/100.0 
            (i.e. if the SNRminus is within the tolerance)
        saveSNRflag : True/False, optional
            Save columns <filter>_SN in which the S/N for the given exposure 
            time is saved. The default is False.
        texpformat : string formatter, optional
            String formatter for the exposure time columns. The default is 
            None. If None, then the default formatter self.texpformat is used
        SNRformat : string formatter, optional
            String formatter for the SNR columns. The default is None. If None,
            then the default formatter self.SNRformat is used

        Returns
        -------
        None. the table with the exposure times is saved as a pdastro object 
        in self.texp. The panda table is in self.texp.t. 
        The formatters for the table are in self.formatters4texptable
        Saving the table:
        self.texp.write('myfilename.txt',formatters=self.formatters4texptable)

        """

        #Make sure filters is a list and not a string
        if isinstance(filters,str):filters=[filters]

        # initialize columns
        #cols=[]
        for col in filters: 
            #cols.append(col+'_t')
            self.texp.t[col+'_t']=None
        if saveSNRflag:
            for col in filters: 
                #cols.append(col+'_SN')
                self.texp.t[col+'_SN']=None

        # get the formats for the columns
        if SNRformat == None: SNRformat = self.SNRformat
        if texpformat == None: texpformat = self.texpformat
        if magformat == None: magformat = self.magformat

        self.formatters4texptable = {}
        
        for filt in filters:
            
            self.formatters4texptable[filt+'_t']=texpformat
            self.formatters4texptable[filt]=magformat
            if saveSNRflag:
                self.formatters4texptable[filt+'_SN']=SNRformat
            
            texps=[]
            SNRs=[]
            mags = self.texp.t[filt]
            print('VVVV',filt,mags)
            self.texp.write()
            for mag in mags:
                texp_results=self.texp4SNRatmag(filt,mag,SNR,lambkg4ETC=lambkg4ETC,
                                                SNR_tolerance_in_percent=SNR_tolerance_in_percent,
                                                spec=spec,**kwargs)
                texps.append(texp_results[texp_type][0])
                if saveSNRflag:
                    SNRs.append(texp_results[texp_type][1])
                    
            self.texp.t[filt+'_t']=np.array(texps)
            if saveSNRflag:
                self.texp.t[filt+'_SN']=np.array(SNRs)


if __name__ == '__main__':
    # initialize with instrument and mode
    lc2texp=lc2texpclass()

    parser = argparse.ArgumentParser(usage="create exposure time table for a given set of filters, mags, and target S/N")
    parser.add_argument('SNR', type=float, help=('specify target SNR'))
    #parser.add_argument("--refmag", nargs=3, help=("reference mag, reference filter, reference time"))
    #parser.add_argument("--distance", default=None, type=float, help=("distance to source in Mpc"))
    #parser.add_argument("--redshift", default=None, type=float, help=("redshift to source"))
    parser.add_argument('-i','--instrument', default='nircam', choices=['nircam','miri','niriss'], help=('specify instrument (default=%(default)s)'))
    parser.add_argument('--mode', default=None, choices=['imaging','sw_imaging','lw_imaging'], help=('specify mode. If None, then the default mode for a given instrument is chosen (default=%(default)s)'))
    parser.add_argument('-f','--filters', default=['F200W'],nargs='+', help=('specify filters'))
#    parser.add_argument('-m','--magrange', nargs=3, type=float, default=[24,28,1], help=('specify the magnitude range magmin magmax dm (default=%(default)s)'))
    parser.add_argument('-s','--save', nargs='*',  type=str, default=None, help=('save the table. If no filename specified, then SNR_<exptime>sec.txt is used (default=%(default)s)'))
    parser.add_argument('--bkg_target',  type=str, default='CDF-S', help=('specify the background target (default=%(default)s)'))
    parser.add_argument('--bkg_targetposition', nargs='+', help=('specify the background position in RA and Dec, overwriting --bkg_target. optional add name of position (default=%(default)s)'))
    parser.add_argument('--bkg_percentile',  type=float, default=50.0, help=('specify the background percentile at the given position (default=%(default)s)'))
    parser.add_argument('--bkg_lam4percentile', type=float, default=4.5, help=('specify the background percentile at the given position (default=%(default)s)'))
    parser.add_argument('--bkg_lam',  type=float, default=4.5, help=('specify the wavelength passed to pandeia jwst_backgrounds.background routine (default=%(default)s)'))
    parser.add_argument('--bkg_thresh',  type=float, default=1.1, help=('specify the threshold passed to pandeia jwst_backgrounds.background routine (default=%(default)s)'))
    parser.add_argument('--SNR_tolerance_in_percent', type=float, default=10.0, help=('specify the tolerance in target S/N (default=%(default)s)'))
    parser.add_argument('--saveSNR', default=False, action='store_true', help=('Save the S/N as well in the table'))

    parser.add_argument("--dt_discovery", type=float, default=None, help=("time difference between KN and GW discovery in days. If not specified, then it is assumed that dt_KN=reftime"))
    parser.add_argument("--dt_trigger", type=float, default=0.1, help=("time difference between HST trigger and KN discovery in days"))
    parser.add_argument("--phasemin", type=float, default=0.5, help=("minimum phase since GW discovery in days (default=%(default)s)"))
    parser.add_argument("--phasemax", type=float, default=14.1, help=("maximum phase since GW discovery in days (default=%(default)s)"))
    parser.add_argument("--phasestep", type=float, default=0.5, help=("phase stepsize in days (default=%(default)s)"))

    parser = lc2texp.lc.define_args(parser=parser)

                        
    args = parser.parse_args()
    
    # initialize with instrument and mode
    lc2texp.initialize_pandeia(instrument=args.instrument,mode=args.mode)
    lc2texp.verbose=args.verbose
    lc2texp.lc.verbose=args.verbose
    
    lcfilename = lc2texp.lc.findlcfilename(args.lcrootdir,args.lcmodel,args.lcparams)
    lc2texp.lc.loadmodellc(lcfilename)
    lc2texp.lc.write()

    filters = args.filters
    phaserange = np.arange(args.phasemin,args.phasemax,args.phasestep)
    print(phaserange)
    
    # set the background
    # Note: targetpos overwrites target.
    lc2texp.set_background4jwst(args.bkg_percentile,
                                lam=args.bkg_lam,thresh=args.bkg_thresh,
                                lam4percentile=args.bkg_lam4percentile,
                                target=args.bkg_target,targetpos=args.bkg_targetposition)
        
    
    # get the splined normalized mag table
    if isinstance(args.refmag,list):
        lc2texp.texp = lc2texp.lc.get_normalized_lc(args.refmag[2],args.refmag[1],args.refmag[0],
                                                    phaserange,filters,maxmag=30.0)
    elif isinstance(args.distance,float) | isinstance(args.redshift,float):
        lc2texp.texp = lc2texp.lc.distance_scalling(phaserange,filters,args.distance,
                                                    args.redshift,maxmag=30.0)
    else:
        raise ValueError('Either a distance or reference magnitude must be specified!')

    # exposure time panda table is in lc2texp.texp.t
    lc2texp.lc2texp_table(filters,args.SNR,
                          SNR_tolerance_in_percent=args.SNR_tolerance_in_percent,
                          saveSNRflag=args.saveSNR)
   
    # save the file if wanted
    if not(args.save is None):
        # if verbose, also write it to screen
        if lc2texp.verbose: lc2texp.texp.write(formatters=lc2texp.formatters4texptable)

        # get the filename
        if args.save ==[]:
            if isinstance(args.reffilter,str):
                filename = 'lc2texp_%s_phase%.1f_%s_%.1f_SNR%.0f.txt' % (lc2texp.instrument,
                                                                     args.refmag[2],args.refmag[1],args.refmag[0],
                                                                     args.SNR)
            if isinstance(args.distance,float):
                filename = 'lc2texp_%s_%fMpc_SNR%.0f.txt' % (lc2texp.instrument,args.distance,args.SNR)
            elif isinstance(args.redshift,float):
                filename = 'lc2texp_%s_%fz_SNR%.0f.txt' % (lc2texp.instrument,args.redshift,args.SNR)
        else:
            filename = args.save[0]

        # save the table
        print('Saving table into %s' % filename)
        lc2texp.texp.write(filename,formatters=lc2texp.formatters4texptable)
    else:
        # if not saved, write it to console
        lc2texp.texp.write(formatters=lc2texp.formatters4texptable)