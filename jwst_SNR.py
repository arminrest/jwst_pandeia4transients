#!/usr/bin/env python3
"""
Created on Thu Oct 15 13:07:57 2020

@author: arest

some definitions based on Dan Coe's jupyter notebooks:
    https://github.com/dancoe/pandeia-imaging/blob/master/depth/NIRCam%20Imaging%20Depth.ipynb
    https://github.com/spacetelescope/nircam_calib/blob/master/nircam_calib/training_notebooks/NIRCam%20Imaging%20Depth%20vs.%20Background.ipynb

"""
import numpy as np
import math
import sys,socket,os,re
import pandas as pd
from pdastro import pdastroclass
import io
from readoutpatterns import readoutpatternclass
from pandeia.engine.perform_calculation import perform_calculation
from pandeia.engine.calc_utils import build_default_calc  # or alternatively, load ETC calculation parameters
from jwst_backgrounds import jbt  # To calculate background
import astropy.io.fits as pyfits  # To load background
import astropy
from background4jwst import background4jwstclass

def nJytoAB(F_nJy):
    return (F_nJy * u.nJy).to(u.ABmag).value


            
class jwst_SNRclass:
    def __init__(self,instrument='nircam',mode='sw_imaging',ETCjsonfile=None):
        self.verbose = 0
        self.allowed_instruments = ['nircam','nirspec','niriss','miri']

        # initialization
        self.initialize_pandeia(instrument,mode,ETCjsonfile=ETCjsonfile)

        # set apertures
        # pandeia default values: 0.1" NIRCam, 0.3" MIRI
        self.aperture_radii = {}
        self.set_apertures({'nircam_sw_imaging':0.08, 'nircam_lw_imaging':0.16, 'niriss':0.16, 'miri':0.3})
        
        # set sky_annulus
        # pandeia default values: 0.22" - 0.4"
        self.sky_annulii = {}
        self.set_sky_annulus({'nircam_sw_imaging':(0.6, 0.99), 'nircam_lw_imaging':(0.6, 0.99), 'niriss':(0.6, 0.99), 'miri':(0.6, 0.99)})
        

        self.background4jwst = background4jwstclass()
        self.lambkg4ETC=None
        
        self.ETCresults = None
       
        
    def get_val4dict(self,d,instrument=None,mode=None):
        if instrument is None: instrument=self.get_instrument()
        if mode is None: mode=self.get_mode()
        if instrument=='nircam':
            s = '%s_%s' % (instrument,mode)
            if not(s in d):
                raise RuntimeError('Could not find %s as a key' % (s))
            val = d[s]
        elif instrument in ['niriss','miri']:
            if not(instrument in d):
                raise RuntimeError('Could not find instrument=%s as a key' % (instrument))
            val = d[instrument]
        else:
            raise RuntimeError("instrument %s is not allowed, only nircam, niriss, and miri" % (instrument))
        return(val)

    def get_aperture(self,instrument=None,mode=None):
        return(self.get_val4dict(self.aperture_radii,instrument=instrument,mode=mode))   
        
    def get_sky_annulus(self,instrument=None,mode=None):
        return(self.get_val4dict(self.sky_annulii,instrument=instrument,mode=mode))           


    def set_apertures(self,ap_dict):
        for key in ap_dict:
            self.aperture_radii[key]=ap_dict[key]
        return(0)


    def set_sky_annulus(self,sky_dict):
        for key in sky_dict:
            self.sky_annulii[key]=sky_dict[key]
        return(0)
    
    def initialize_pandeia(self,instrument,mode,ETCjsonfile=None):
        print('Initializing pandeia with %s, %s' % (instrument,mode))
        
        instrument = instrument.lower()
        mode = mode.lower()
        if not (instrument in self.allowed_instruments):
            raise(RuntimeError,'instrument %s not in %s' % (instrument,' '.join(self.allowed_instruments)))   

        self.readoutpattern=readoutpatternclass(instrument)
 
        if ETCjsonfile is None:
            self.pandeiacfg=build_default_calc('jwst',instrument,mode)
        else:
            with open(ETCjsonfile) as f:  # use a json file you have
                 self.pandeiacfg = json.load(f)
        return(0)
    
    def get_instrument(self):
        if self.pandeiacfg is None:
            raise RuntimeError("pandeia not yet initialized, cannot get the instrument")
        return(self.pandeiacfg['configuration']['instrument']['instrument'])
    
    def get_mode(self):
        if self.pandeiacfg is None:
            raise RuntimeError("pandeia not yet initialized, cannot get the instrument")
        return(self.pandeiacfg['configuration']['instrument']['mode'])

    def determine_imaging_mode_and_aperture(self,filt,instrument=None):
        if instrument is None:
            instrument = self.get_instrument()
            
        if instrument.lower() in ['niriss','miri']:
            mode = ('imaging',None)
            raise RuntimeError('Need to figure out what aperture value is for niriss and miri!!!! fix me!!!')
        elif instrument.lower() == 'nircam':
            lam = int(filt[1:4]) / 100.
            ch = 'sw lw'.split()[lam > 2.4]
            mode = ch+'_imaging'
        else:
             raise RuntimeError("instrument %s is not allowed, only nircam, niriss, and miri" % (instrument))
           
        return(mode,ch)    

    def set_background4jwst(self,percentile,**kwargs):
        """

        Parameters
        ----------
        percentile : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            optional arguments that get passed to background4jwst.lambkg4ETC:
                lam : float, optional
                    Input wavelength for the pandeia 'background' routine. If None, the default value of self.lam is used
                thresh : float, optional
                    Input threshold for the pandeia 'background' routine. If None, the default value of self.thresh is used
                lam4percentile : float, optional
                    wavelength for which the percentiles are calculated for. If None, the default value of self.lam4percentile is used
                target: string, optional
                    if specified, the RA,Dec is set to the target. target must be part of self.defaulttargets
        Returns
        -------
        None.

        """
        self.lambkg4ETC = self.background4jwst.lambkg4ETC(percentile,**kwargs)
        return(0)

    def Imaging_SNR(self,filt, mag, exptime, lambkg4ETC=None):

        if lambkg4ETC is None:
            lambkg4ETC=self.lambkg4ETC
        if lambkg4ETC is None:
            print('!!!!!!!!!!!!!!!\n!!!WARNING!!!!!\n!!!!!!!!!!!!!!!\nNo background specified, calculating SNR WITHOUT background!!')
            lambkg4ETC=[]

        self.pandeiacfg['background'] = lambkg4ETC

        # Assign filter and magnitude
        self.pandeiacfg['configuration']['instrument']['filter'] = filt.lower()
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_flux'] = mag
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_fluxunit'] = 'abmag'

        # mode and aperture
        (mode,aperture)=self.determine_imaging_mode_and_aperture(filt)
        self.pandeiacfg['configuration']['instrument']['aperture'] = aperture
        self.pandeiacfg['configuration']['instrument']['mode'] = mode

        # Photometric aperture and background sky annulus
        self.pandeiacfg['strategy']['aperture_size'] = self.get_aperture()
        self.pandeiacfg['strategy']['sky_annulus'] = self.get_sky_annulus()
        
        # Exposure specifications
        info = self.readoutpattern.info4closestexptime(exptime)
        #print(info)
        for key in ['NEXP','NINT','NGROUP','readout_pattern']:
            self.pandeiacfg['configuration']['detector'][key.lower()] = info[key]
            #self.pandeiacfg['configuration']['detector']['ngroup'] = info['NGROUPS']
        #self.pandeiacfg['configuration']['detector']['nexp'] = 
        #self.pandeiacfg['configuration']['detector']['nint'] = 
        #self.pandeiacfg['configuration']['detector']['ngroup'] = 
        #self.pandeiacfg['configuration']['detector']['readout_pattern'] = .lower()
        
    
    
        # RUN CALCULATION
        self.ETCresults = perform_calculation(self.pandeiacfg)  
        SNR = self.ETCresults['scalar']['sn']
        total_exposure_time = self.ETCresults['scalar']['total_exposure_time']
        print(SNR, total_exposure_time)
        return SNR, total_exposure_time

    def Imaging_SNR_table(self, filters, magrange, exptime, lambkg4ETC=None):
        cols = ['mag']
        cols.extend(filters)
        self.SNR = pdastroclass(columns=cols)
        self.SNR.t['mag']=magrange
        for filt in filters:
            SNRs=[]
            for mag in magrange:
                (SNRval,total_exposure_time)=self.Imaging_SNR(filt,mag,exptime,lambkg4ETC=lambkg4ETC)
                SNRs.append(SNRval)
            self.SNR.t[filt]=np.array(SNRs)
        #self.SNR.t.format({'mag': '{:.2f}'.format, 'F200W': '{:.2f}'.format})
        #self.SNR.write(formatters={'mag': '{:.2f}'.format, 'F200W': '{:.2f}'.format})

if __name__ == '__main__':
    print('hello')
    jwst_SNR=jwst_SNRclass(instrument='nircam',mode='sw_imaging')
    print(jwst_SNR.aperture_radii)
    print(jwst_SNR.sky_annulii)
    
    jwst_SNR.set_background4jwst(50,target='EmptyERS')

    jwst_SNR.Imaging_SNR_table(['F200W'],np.arange(28.0,29.0,0.5),1200)
    
