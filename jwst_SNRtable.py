#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:48:36 2020

@author: arest
"""

import argparse,sys
import numpy as np
from jwst_SNR import jwst_SNRclass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="create S/N table for a given set of filters, mags, and exposure time")
    parser.add_argument('exptime', type=float, help=('specify exposure time (approximate ok, it will identify the closest readout pattern'))
    parser.add_argument('-i','--instrument', default='nircam', choices=['nircam','miri','niriss'], help=('specify instrument (default=%(default)s)'))
    parser.add_argument('-f','--filters', nargs='+', help=('specify filters'))
    parser.add_argument('-m','--magrange', nargs=3, type=float, default=[24,28,1], help=('specify the magnitude range magmin magmax dm (default=%(default)s)'))
    parser.add_argument('-s','--save', nargs='*',  type=str, default=None, help=('save the table. If no filename specified, then SNR_<exptime>sec.txt is used (default=%(default)s)'))
    parser.add_argument('--bkg_target',  type=str, default='EmptyERS', help=('specify the background target (default=%(default)s)'))
    parser.add_argument('--bkg_position', nargs='+', help=('specify the background position in RA and Dec, overwriting --bkg_target. optional add name of position (default=%(default)s)'))
    parser.add_argument('--bkg_percentile',  type=float, default=50.0, help=('specify the background percentile at the given position (default=%(default)s)'))
    parser.add_argument('--bkg_lam4percentile', type=float, default=4.5, help=('specify the background percentile at the given position (default=%(default)s)'))
    parser.add_argument('--bkg_lam',  type=float, default=4.5, help=('specify the wavelength passed to pandeia jwst_backgrounds.background routine (default=%(default)s)'))
    parser.add_argument('--bkg_thresh',  type=float, default=1.1, help=('specify the threshold passed to pandeia jwst_backgrounds.background routine (default=%(default)s)'))
    parser.add_argument('-v','--verbose', default=0, action='count')
                        
    args = parser.parse_args()
    
    if args.instrument=='nircam':
        mode='sw_imaging'
    else:
        mode='imaging'
        
    # initialize with instrument and mode
    jwst_SNR=jwst_SNRclass(instrument=args.instrument,mode=mode)
    jwst_SNR.verbose=args.verbose

    # set the background
    if args.bkg_position is None:
        # use the pre-defined background
        jwst_SNR.set_background4jwst(args.bkg_percentile,
                                     lam=args.bkg_lam,thresh=args.bkg_thresh,
                                     lam4percentile=args.bkg_lam4percentile,target=args.bkg_target)
    else:
        # set the position ....
        if len(args.bkg_position)<2 or len(args.bkg_position)>3:
            raise RuntimeError('only 2 or 3 arguments allowed for --bkg_position')
        if len(args.bkg_position)==3:
            name = args.bkg_position[-1]
        else:
            name = 'usertarget'
        jwst_SNR.background4jwst.set_position(args.bkg_position[0],args.bkg_position[1],name)
        # then calculate teh background
        jwst_SNR.lambkg4ETC = jwst_SNR.background4jwst.lambkg4ETC(args.bkg_percentile,
                                                                  lam=args.bkg_lam,thresh=args.bkg_thresh,
                                                                  lam4percentile=args.bkg_lam4percentile)
        
    
    filters=['F200W'] # just have a default filter for now
    if not(args.filters is None):
        filters = args.filters
        
    magrange = np.arange(args.magrange[0],args.magrange[1],args.magrange[2])
    exptime = jwst_SNR.Imaging_SNR_table(filters,magrange,args.exptime)
    if not(args.save is None):
        if jwst_SNR.verbose: jwst_SNR.SNR.write(formatters=jwst_SNR.formatters4SNRtable)
        if args.save ==[]:
            filename = 'SNR_%.0fsec.txt' % (exptime)
        else:
            filename = args.save[0]
        print('Saving table into %s' % filename)
        jwst_SNR.SNR.write(filename,formatters=jwst_SNR.formatters4SNRtable)
    else:
        jwst_SNR.SNR.write(formatters=jwst_SNR.formatters4SNRtable)