#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:48:36 2020

@author: arest
"""

import argparse,sys
import numpy as np
from jwst_SNR import jwst_SNRclass
import pysynphot as S 
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="create exposure time table for a given set of filters, mags, and target S/N")
    parser.add_argument('SNR', type=float, help=('specify target SNR'))
    parser.add_argument('spec', type=str, help=('input spectrum'))
    parser.add_argument('lam', type=float, default=[2], nargs='+', help=('reference wavelength in micron'))
    parser.add_argument('width', type=float, default=[.5],nargs='+', help=('width to average around reference wavelength fo SNR.'))
    parser.add_argument('-i','--instrument', default='nirspec', choices=['nirspec'], help=('specify instrument (default=%(default)s)'))
    parser.add_argument('--mode', default='fixed_slit', choices=['fixed_slit','ifu'], help=('specify mode. If None, then the default mode for a given instrument is chosen (default=%(default)s)'))
    parser.add_argument('-g','--grating', default=['prism'],nargs='+', help=('specify nirspec disperser'))
    parser.add_argument('-f','--reffilter', default='F200W', help=('specify normalisation filter'))
    parser.add_argument('-m','--magrange', nargs=3, type=float, default=[24,28,1], help=('specify the magnitude range magmin magmax dm (default=%(default)s)'))
    parser.add_argument('-s','--save', nargs='*',  type=str, default=None, help=('save the table. If no filename specified, then SNR_<exptime>sec.txt is used (default=%(default)s)'))
    parser.add_argument('--bkg_target',  type=str, default='CDF-S', help=('specify the background target (default=%(default)s)'))
    parser.add_argument('--bkg_targetposition', nargs='+', help=('specify the background position in RA and Dec, overwriting --bkg_target. optional add name of position (default=%(default)s)'))
    parser.add_argument('--bkg_percentile',  type=float, default=50.0, help=('specify the background percentile at the given position (default=%(default)s)'))
    parser.add_argument('--bkg_lam4percentile', type=float, default=4.5, help=('specify the background percentile at the given position (default=%(default)s)'))
    parser.add_argument('--bkg_lam',  type=float, default=4.5, help=('specify the wavelength passed to pandeia jwst_backgrounds.background routine (default=%(default)s)'))
    parser.add_argument('--bkg_thresh',  type=float, default=1.1, help=('specify the threshold passed to pandeia jwst_backgrounds.background routine (default=%(default)s)'))
    parser.add_argument('-v','--verbose', default=0, action='count')
    parser.add_argument('--SNR_tolerance_in_percent', type=float, default=10.0, help=('specify the tolerance in target S/N (default=%(default)s)'))
    parser.add_argument('--saveSNR', default=False, action='store_true', help=('Save the S/N as well in the table'))
                        
    args = parser.parse_args()
    spec_file - args.spec
    if spec_file is None:
        raise(ValueError('Spectrum must be provided'))
    else:
        spec = pd.read_csv(spec_file)
        if np.nanmedian(spec['wave'].values) < 100:
            waveunits  = 'micron'
        elif np.nanmedian(spec['wave'].values) > 100:
            waveunits = 'angstrom'
        spec = S.ArraySpectrum(spec['wave'].values,spec['flux'].values,waveunits=waveunits,
               fluxunits='flam')
    # initialize with instrument and mode
    jwst_SNR=jwst_SNRclass(instrument=args.instrument,mode=args.mode,grating=args.grating)
    jwst_SNR.verbose=args.verbose

    # set the background
    # Note: targetpos overwrites target.
    jwst_SNR.set_background4jwst(args.bkg_percentile,
                                 lam=args.bkg_lam,thresh=args.bkg_thresh,
                                 lam4percentile=args.bkg_lam4percentile,
                                 target=args.bkg_target,targetpos=args.bkg_targetposition)
        
    
    
    magrange = np.arange(args.magrange[0],args.magrange[1],args.magrange[2])

    # exposure time panda table is in jwst_SNR.texp.t
    jwst_SNR.Spec_texp_table(arg.grating,arg.reffilter,magrange,args.SNR,
                                SNR_tolerance_in_percent=args.SNR_tolerance_in_percent,
                                saveSNRflag=args.saveSNR,spec=spec,wave=args.lam,width=args.width)
   
    # save the file if wanted
    if not(args.save is None):
        # if verbose, also write it to screen
        if jwst_SNR.verbose: jwst_SNR.texp.write(formatters=jwst_SNR.formatters4texptable)

        # get the filename
        if args.save ==[]:
            filename = 'texp_SNR%.0f.txt' % (args.SNR)
        else:
            filename = args.save[0]

        # save the table
        print('Saving table into %s' % filename)
        jwst_SNR.texp.write(filename,formatters=jwst_SNR.formatters4texptable)
    else:
        # if not saved, write it to console
        jwst_SNR.texp.write(formatters=jwst_SNR.formatters4texptable)