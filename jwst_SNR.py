#!/usr/bin/env python3
"""
Created on Thu Oct 15 13:07:57 2020

@author: arest

some definitions based on Dan Coe's jupyter notebooks:
    https://github.com/dancoe/pandeia-imaging/blob/master/depth/NIRCam%20Imaging%20Depth.ipynb
    https://github.com/spacetelescope/nircam_calib/blob/master/nircam_calib/training_notebooks/NIRCam%20Imaging%20Depth%20vs.%20Background.ipynb

"""
import numpy as np
import matplotlib.pyplot as plt
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
    def __init__(self,instrument='nircam',mode=None,ETCjsonfile=None):
        self.verbose = 0
        self.allowed_instruments = ['nircam','nirspec','niriss','miri']
        self.instrument = instrument
        # initialization
        self.initialize_pandeia(instrument,mode=mode,ETCjsonfile=ETCjsonfile)

        # set apertures
        # pandeia default values: 0.1" NIRCam, 0.3" MIRI
        self.aperture_radii = {}
        self.set_apertures({'nircam_sw_imaging':0.08, 'nircam_lw_imaging':0.16, 
                            'niriss':0.16, 'miri':0.3,'nirspec':0.15})
        
        # set sky_annulus
        # pandeia default values: 0.22" - 0.4"
        self.sky_annulii = {}
        self.set_sky_annulus({'nircam_sw_imaging':(0.6, 0.99), 'nircam_lw_imaging':(0.6, 0.99), 
                                'niriss':(0.6, 0.99), 'miri':(0.6, 0.99),'nirspec':(0.3,0.5)})
        

        self.background4jwst = background4jwstclass()
        self.lambkg4ETC=None
        
        self.ETCresults = None
        
        self.SNRformat = '{:.2f}'.format
        self.texpformat = '{:.1f}'.format
        self.formatters4SNRtable = None
        self.formatters4texptable = None

       
        
    def get_val4dict(self,d,instrument=None,mode=None):
        if instrument is None: instrument=self.get_instrument()
        if mode is None: mode=self.get_mode()
        if instrument=='nircam':
            s = '%s_%s' % (instrument,mode)
            if not(s in d):
                raise RuntimeError('Could not find %s as a key' % (s))
            val = d[s]
        elif instrument in ['niriss','miri','nirspec']:
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
    
    def initialize_pandeia(self,instrument,mode=None,ETCjsonfile=None):
        print('Initializing pandeia with %s, %s' % (instrument,mode))
        
        # make sure the instrument is correct
        instrument = instrument.lower()
        if not (instrument in self.allowed_instruments):
            raise(RuntimeError,'instrument %s not in %s' % (instrument,' '.join(self.allowed_instruments)))   

        # if no mode is given, choose the default one.
        if mode is None:
            if instrument=='nircam':
                mode='sw_imaging'
            elif instrument in ['niriss','miri']:
                mode='imaging'
            elif instrument=='nirspec':
                mode='fixed_slit'
            else:
                raise RuntimeError('unknown instrument %s!' % instrument)
        mode = mode.lower()

        self.readoutpattern=readoutpatternclass(instrument)
 
        if ETCjsonfile is None:
            print('Initializing',instrument,mode)
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
            mode = 'imaging'
            aperture = 'imager'
            if instrument.lower() in ['niriss']:
                raise RuntimeError('Need to figure out what aperture value is for niriss and miri!!!! fix me!!!')
        elif instrument.lower() == 'nircam':
            lam = int(filt[1:4]) / 100.
            aperture = 'sw lw'.split()[lam > 2.4]
            mode = aperture+'_imaging'
        elif instrument.lower() == 'nirspec':
            mode = 'fixed_slit'
            aperture = 's200a1'
        else:
             raise RuntimeError("instrument %s is not allowed, only nircam, niriss, miri, and nirspec" % (instrument))
           
        return(mode,aperture)    

    def set_background4jwst(self,percentile,**kwargs):
        """

        Parameters
        ----------
        percentile : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            optional arguments that get passed to background4jwst.lambkg4ETC:
                lam : float, optional
                    Input wavelength for the pandeia 'background' routine. 
                    If None, the default value of self.lam is used
                thresh : float, optional
                    Input threshold for the pandeia 'background' routine. 
                    If None, the default value of self.thresh is used
                lam4percentile : float, optional
                    wavelength for which the percentiles are calculated for. 
                    If None, the default value of self.lam4percentile is used
                target: string, optional
                    if specified, the RA,Dec is set to the target. 
                    Target must be part of self.defaulttargets
        Returns
        -------
        None.

        """
        self.lambkg4ETC = self.background4jwst.lambkg4ETC(percentile,**kwargs)
        return(0)

    def Imaging_SNR(self, filt, mag, exptime, lambkg4ETC=None,spec=None,**kwargs):
        """
        Parameters
        ----------
        filt : string
            filter for which the S/N should be calculated for.
        mag : float
            magnitude  for which the S/N should be calculated for.
        exptime : float
             exposure time for which the S/N should be calculated for. The
             readout pattern that closest matches the exposure time is used.
        lambkg4ETC : list/tuple of two arrays, optional
            the first array is the wavelength, the second the background. 
            The default is None. If none, then the background data in 
            self.lambkg4ETC is used.

        Returns
        -------
        (SNR,exptime)
        Note: the returned exptime can differ from the passed exptime. 
        It is the true exposure time associated with the readout pattern
        used for the calculation.

        """
        
        if self.verbose>2: print('Calculating SNR for filter:%s mag:%f, exptime:%f \n' % (filt, mag, exptime))

        if lambkg4ETC is None:
            if self.verbose>2: print('Using saved lambkg4ETC')
            lambkg4ETC=self.lambkg4ETC
        if lambkg4ETC is None:
            print('!!!!!!!!!!!!!!!\n!!!WARNING!!!!!\n!!!!!!!!!!!!!!!\nNo background specified, calculating SNR WITHOUT background!!')
            lambkg4ETC=[]

        self.pandeiacfg['background'] = lambkg4ETC

        # Assign filter and magnitude
        if self.instrument != 'nirspec':
            self.pandeiacfg['configuration']['instrument']['filter'] = filt.lower()
        else:
            self.pandeiacfg['configuration']['detector']['subarray'] = 'full'
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['bandpass'] = 'nircam,sw_imaging,' + filt.lower()
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_flux'] = mag
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_fluxunit'] = 'abmag'
        if not(spec is None):
            print('reference spec')
            # spectrum needs to be 
            spec.convert('micron')
            spec.convert('mJy')
            self.pandeiacfg['scene'][0]['spectrum']['sed']['sed_type'] = 'input'
            self.pandeiacfg['scene'][0]['spectrum']['sed']['spectrum'] = [spec.wave,spec.flux]
            self.pandeiacfg['scene'][0]['spectrum']['sed']['unit'] = 'flam'

        # mode and aperture
        (mode,aperture)=self.determine_imaging_mode_and_aperture(filt)
        if not(aperture is None):
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
        

        # calculate SNR with pandeia
        self.ETCresults = perform_calculation(self.pandeiacfg)  
        
        if self.instrument=='nirspec':
            SNR = self.Av_spec_SNR(**kwargs)
        else:
            SNR = self.ETCresults['scalar']['sn']
        total_exposure_time = self.ETCresults['scalar']['total_exposure_time']
        
        if self.verbose>1: print('filter:%s mag:%.2f, target exptime:%.1f  ==> SNR=%.2f exptime=%.1f' % (filt, mag, exptime,SNR,total_exposure_time))
        return(SNR,total_exposure_time)

    def Av_spec_SNR(self,wave,width=0):
        lam, snr = self.ETCresults['1d']['sn']
        #plt.figure()
        #plt.axvspan(wave-width,wave+width,alpha=.3,color='orange')
        #plt.plot(lam,snr)

        low = wave - width
        high = wave + width
        ind = (lam >= low) & (lam <= high)
        av = np.nanmean(snr[ind])
        return av

    def texp4SNRatmag(self,filt,mag,SNR,lambkg4ETC=None,SNR_tolerance_in_percent=10.0,
                        spec=None,**kwargs):
        """
        Parameters
        ----------
        filt : string
            filter for which the exposure time should be calculated for.
        mag : float
            magnitude  for which the exposure time should be calculated for.
        SNR : float 
            S/N  for which the exposure time should be calculated for.
        lambkg4ETC : list/tuple of two arrays, optional
            the first array is the wavelength, the second the background. 
            The default is None. If none, then the background data in 
            self.lambkg4ETC is used.
        SNR_tolerance_in_percent : float, optional
            An exposure time with a SNRminus<SNR can be accepted as 'best' if
            SNRminus is closer to SNR than SNRplus (or SNRplus==None), 
            AND if (SNR0-SNRminus)/SNR0<SNR_tolerance_in_percent/100.0 
            (i.e. if the SNRminus is within the tolerance)

        Raises
        ------
        RuntimeError
            DESCRIPTION.

        Returns
        -------
        Returns a dictionary with 3 keys: 'plus', 'minus' and 'best'.
        'plus' contains the (texpplus,SNRplus) with SNR>=SNR0
        'minus' contains the (texpminus,SNRminus) with SNR<SNR0
        if the lowest or highest exptime is hit, then they are set to (None,None) accordingly
        
        'best' contains the  (texpbest,SNRbest) for the following criteria:
             if SNRminus is closer to SNR0 than SNRplus (or SNRplus==None), 
             AND if (SNR0-SNRminus)/SNR0<SNR_tolerance_in_percent/100.0 
             (i.e. if the SNRminus is within the tolerance), then best is 
             set to minus. Otherwise best=plus
        """
        if self.verbose: print('#############################\n#### Filter %s, mag %.2f for S/N=%.f \n#############################' % (filt,mag,SNR))

        texp0=1000
        (SNR0, texp0) = self.Imaging_SNR(filt,mag,texp0,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)        
        if self.verbose>1: print('SNR=%6.2f for starting texp=%6.1f' % (SNR0,texp0))
        
        instrument = self.get_instrument()
        if instrument=='nircam':
            pwlindex = 9/8
        elif instrument=='miri':
            pwlindex = 16/8
        elif instrument == 'nirspec':
            pwlindex = 11/8 
        else:
            raise RuntimeError('instrment %s not yet implemented!' % instrument)
            
      
        # guesstimate the best exposure time. The SNR theoretically goes with sqrt(t), 
        # thus t should go with SNR^2. However, it looks like the exponent is smaller than 2
        
        texp_guess = texp0 * math.pow(SNR/SNR0,pwlindex)
        if self.verbose>1: print('texp guess: %.1f' % texp_guess)
        
        tnext = self.readoutpattern.nextbiggestexptime(texp_guess)
        if tnext is None:
            tnext = self.readoutpattern.nextsmallestexptime(texp_guess)
        if tnext is None:
            raise RuntimeError("BUG??? Cannot find any exposure time...")
            
            
        (SNRnext,tnext) = self.Imaging_SNR(filt,mag,tnext,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)        
        if self.verbose>1: print('SNR=%6.2f for next texp=%6.1f' % (SNRnext,tnext))
        
        (tlast,SNRlast)=(tnext,SNRnext)
        if SNRnext<=SNR:
            while (SNRnext<SNR):
                if self.verbose>1: print('SNR=%6.2f<%.2f for texp=%6.1f, checking the next larger exptime...' % (SNRnext,SNR,tnext))
                #print('SNR=%6.2f for texp=%6.1f, SNR=%.2f wanted...' % (SNRnext,tnext,SNR))
                  
                (tlast,SNRlast)=(tnext,SNRnext)
                 
                # get the next bigger exposure time
                tnext = self.readoutpattern.nextbiggestexptime(tnext+1.0)
                
                # If this is the last index of the readoutpattern, and the SNR is still not big enough, return None and the last SNR possible
                if (tnext is None):
                    SNRnext=None
                    print('Warning: could not find an exposure time that is long enough to reach SNR=%.2f, only %.2f' % (SNR,SNRlast))
                    break

                (SNRnext,tnext) = self.Imaging_SNR(filt,mag,tnext,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)        
            if not(tnext is None):
                if self.verbose: print('SNR=%6.2f>=%.2f for texp=%6.1f!! SUCCESS!' % (SNRnext,SNR,tnext))
            (tplus,SNRplus)=(tnext,SNRnext)
            (tminus,SNRminus)=(tlast,SNRlast)
            #return(tnext,SNRnext)
        else:
            while (SNRnext>SNR):
                if self.verbose>1: print('SNR=%6.2f>%.2f for texp=%6.1f, checking the next lower exptime...' % (SNRnext,SNR,tnext))
                 
                (tlast,SNRlast)=(tnext,SNRnext)
                 
                # get the next smaller exposure time
                tnext = self.readoutpattern.nextsmallestexptime(tnext-1.0)
                
                # If this is the last index of the readoutpattern, and the SNR is still not big enough, return None and the last SNR possible
                if (tnext is None):
                    SNRnext=None
                    print('Warning: could not find an exposure time that is short enough to be below SNR=%.2f, thus keeping texp=%.1f and SNR=%.2f' % (SNR,tlast,SNRlast))
                    break
                    
                (SNRnext,tnext) = self.Imaging_SNR(filt,mag,tnext,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)
            if not(tnext is None):
                if self.verbose: print('SNR=%6.2f>=%.2f for texp=%6.1f, and SNR=%.2f<%.2f for texp=%6.1f!! SUCCESS!' % (SNRlast,SNR,tlast,SNRnext,SNR,tnext))
            (tplus,SNRplus)=(tlast,SNRlast)
            (tminus,SNRminus)=(tnext,SNRnext)

        # find the exposure time 
        if tplus is None:
            # not good, didn't find an exposure time that gives enough SNR.
            # checking if minus is within tolerance, if not setting best to None
            if (SNR-SNRminus)/SNR<SNR_tolerance_in_percent/100.0:
                (tbest,SNRbest) = (tminus,SNRminus)
            else:   
                (tbest,SNRbest) = (None,None)
        elif tminus is None:
            # the shortest exposure time gives enough SNR! We can set best to plus
            (tbest,SNRbest) = (tplus,SNRplus)
        else:
            # find the exposure times for which the SNR is best to the desired SNR. If it is the one with
            # a SNR that is smaller than the desired SNR, make sure it is within the tolerance
            if ((SNRplus-SNR)>(SNR-SNRminus)) and ((SNR-SNRminus)/SNR<SNR_tolerance_in_percent/100.0):
                (tbest,SNRbest) = (tminus,SNRminus)
            else:
                (tbest,SNRbest) = (tplus,SNRplus)
               
        results = {'plus':(tplus,SNRplus),
                   'minus':(tminus,SNRminus),
                   'best':(tbest,SNRbest)}    
        
        return(results)
                                

    def Imaging_SNR_table(self, filters, magrange, exptime, lambkg4ETC=None,SNRformat=None,
                            spec=None,**kwargs):
        """
        Parameters
        ----------
        filters : list or string
            Pass (a list of) filter(s) for which to calculate the SNR.
        magrange : list/tuple of magnitudes
            Pass a list of magnitudes for which to calculate the SNR..
        exptime : float
             exposure time for which the S/N should be calculated for. The
             readout pattern that closest matches the exposure time is used.
        lambkg4ETC : list/tuple of two arrays, optional
            the first array is the wavelength, the second the background. 
            The default is None. If none, then the background data in 
            self.lambkg4ETC is used.
        SNRformat : string formatter, optional
            String formatter for the SNR table. The default is None. If None,
            then the default formatter self.SNRformat is used

        Returns
        -------
        None. the table with the SNRs is saved as a pdastro object 
        in self.SNR. The panda table is in self.SNR.t. 
        The formatters for the table are in self.formatters4SNRtable
        Saving the table:
        self.SNR.write('myfilename.txt',formatters=self.formatters4SNRtable)
        """
        
        #Make sure filters is a list and not a string
        if isinstance(filters,str):filters=[filters]
        
        cols = ['mag']
        #cols.extend(filters)
        for col in filters: cols.append(col+'_SN')
       
        if SNRformat == None: SNRformat = self.SNRformat
        
        self.SNR = pdastroclass(columns=cols)
        self.SNR.t['mag']=magrange
        self.formatters4SNRtable = {}
        
        exptime1=None
        for filt in filters:
            self.formatters4SNRtable[filt+'_SN']=SNRformat
            
            SNRs=[]
            for mag in magrange:
                (SNRval,exptime1)=self.Imaging_SNR(filt,mag,exptime,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)
                SNRs.append(SNRval)
                
            self.SNR.t[filt+'_SN']=np.array(SNRs)
            
        return(exptime1)

    def Imaging_texp_table(self, filters, magrange, SNR, lambkg4ETC=None,texp_type='best',
                           SNR_tolerance_in_percent=10,saveSNRflag=False,texpformat=None,SNRformat=None,
                           spec=None,**kwargs):
        """

        Parameters
        ----------
        filters : list or string
            Pass (a list of) filter(s) for which to calculate the SNR.
        magrange : list/tuple of magnitudes
            Pass a list of magnitudes for which to calculate the SNR..
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

        cols = ['mag']
        for col in filters: cols.append(col+'_t')
        if saveSNRflag:
            for col in filters: cols.append(col+'_SN')

        if SNRformat == None: SNRformat = self.SNRformat
        if texpformat == None: texpformat = self.texpformat

        self.texp = pdastroclass(columns=cols)
        self.texp.t['mag']=magrange
        self.formatters4texptable = {}
        
        for filt in filters:
            
            self.formatters4texptable[filt+'_t']=texpformat
            if saveSNRflag:
                self.formatters4texptable[filt+'_SN']=SNRformat
            
            texps=[]
            SNRs=[]
            for mag in magrange:
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
    print('hello')
    jwst_SNR=jwst_SNRclass(instrument='nircam',mode='sw_imaging')
    print(jwst_SNR.aperture_radii)
    print(jwst_SNR.sky_annulii)
    
    jwst_SNR.verbose=1
    jwst_SNR.set_background4jwst(50,target='EmptyERS')

    #jwst_SNR.Imaging_SNR_table(['F200W'],np.arange(28.0,29.0,0.5),1200)
    #jwst_SNR.texp4SNRatmag('F200W',28.0,20.0)
    jwst_SNR.Imaging_texp_table(['F200W'],np.arange(28.0,29.0,0.5),10)
    
