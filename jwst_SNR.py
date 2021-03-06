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

grating_range = {'g140m':[0.97,1.84],
                 'g235m':[1.66,3.07],
                 'g395m':[2.87,5.10],
                 #'g140h':[0.97,1.82],
                 #'g235h':[1.66,3.05],
                 #'g395h':[2.87,5.14],
                 'prism':[0.60,5.30]}

            
class jwst_SNRclass:
    def __init__(self,instrument='nircam',mode=None,ETCjsonfile=None):
        self.verbose = 0
        self.allowed_instruments = ['nircam','nirspec','niriss','miri']
        # book keeping variables
        self.instrument = instrument
        self.spectrum = False
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
        self.magformat = '{:.2f}'.format
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
    
    def set_grating_filter(self,grating):
        
        if (grating == 'g140h') | (grating == 'g140m'):
            message = ('2 filters available: f070lp, or f100lf \n' +
                        'Assigning f070lp')
            #warnings.warn(message)
            filt  = 'f100lp'

        elif (grating == 'g235h') | (grating == 'g235m'):
            filt = 'f170lp'
        elif (grating == 'g395h') | (grating == 'g395m'):
            filt = 'f290lp'
        elif grating == 'prism':
            filt = 'clear'
        return filt

    def Check_grating(self,grating):
        """
        Check to see if the input grating value is a valid option
        """
        allowed = ['prism', 'g140h', 'g140m', 
                   'g235h', 'g235m', 'g395h', 'g395m']
        for i in range(len(allowed)):
            if grating == allowed[i]:
                return grating
        message = ('No such grating available, please choose from:\n '
                   + 'prism\n g140h\n g140m\n g235h\n g235m\n g395h\n g395m')
        raise(ValueError(message))

    def initialize_pandeia(self,instrument,mode=None,grating='prism',ETCjsonfile=None):
        print('Initializing pandeia with %s, %s' % (instrument,mode))
        
        # make sure the instrument is correct
        instrument = instrument.lower()
        self.instrument = instrument
        if not (instrument in self.allowed_instruments):
            raise(RuntimeError,'instrument %s not in %s' % (instrument,' '.join(self.allowed_instruments)))   

        # if no mode is given, choose the default one.
        if mode is None:
            print('no modes')
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

            if instrument == 'nirspec':
                self.pandeiacfg['configuration']['detector']['subarray'] = 'full'
                self.pandeiacfg['configuration']['instrument']['disperser'] = self.Check_grating(grating)
                self.pandeiacfg['configuration']['instrument']['filter'] = self.set_grating_filter(grating)
            #elif instrument == 'miri':
            #    self.pandeiacfg['strategy']['method'] = 'ifunodinscene'
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
            
        if (instrument.lower() in ['niriss','miri']) & ~self.spectrum:
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
        elif (instrument.lower() == 'miri') & self.spectrum:
            mode = 'lrsslit'
            aperture = 'lrsslit'
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

    def setup_for_norm(self,filt):

        sw_imaging = np.array(['f070w','f090w','f115w','f140m','f150w',
                               'f150w2','f164n+f150w2','f164n+f150w2',
                              'f182m','f187n','f200w','f210m','f212n'])
        lw_imaging = np.array(['f250m','f277w','f300m','f322w2','f323n+f322w2',
                              'f335m','f356w','f360','f405n+f444w','f410m',
                              'f430m','f444w','f460m','f466n+f444w',
                              'f470n+f444w','f480m'])
        miri_imaging = np.array(['f560w','f770w','f1000w','f1130w','f1280w',
                                 'f1500w','f1800w','f2100w','f2550w'])

        if (filt == sw_imaging).any():
            i = 'nircam'
            m = 'sw_imaging'
        elif (filt == lw_imaging).any():
            i = 'nircam'
            m = 'lw_imaging'
        
        elif (filt == miri_imaging).any():
            i = 'miri'
            m = 'imaging'
        else:
            message = 'no such imaging filter'
            raise(ValueError(message))

        norm_mode = '{},{},{}'.format(i,m,filt)
        return norm_mode


    def Calculate_SNR(self, filt, mag, exptime, lambkg4ETC=None,spec=None,**kwargs):
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
        if self.spectrum:
            grating = self.pandeiacfg['configuration']['instrument']['disperser']
            if self.verbose>2: print('Calculating SNR for filter:%s mag:%f, %f, exptime:%f \n' % (filt, mag, grating, exptime))
        else:
            if self.verbose>2: print('Calculating SNR for filter:%s mag:%f, exptime:%f \n' % (filt, mag, exptime))

        if lambkg4ETC is None:
            if self.verbose>2: print('Using saved lambkg4ETC')
            lambkg4ETC=self.lambkg4ETC
        if lambkg4ETC is None:
            raise RuntimeError('No background specified!!')

        self.pandeiacfg['background'] = lambkg4ETC

        # Assign filter and magnitude
        if self.instrument != 'nirspec':
            self.pandeiacfg['configuration']['instrument']['filter'] = filt.lower()
            
        if not(filt is None):
            self.pandeiacfg['scene'][0]['spectrum']['normalization']['bandpass'] = self.setup_for_norm(filt.lower())
        else:
            self.pandeiacfg['scene'][0]['spectrum']['normalization']['bandpass'] = None
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_flux'] = mag
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_fluxunit'] = 'abmag'
        self.pandeiacfg['scene'][0]['spectrum']['normalization']['type'] = 'jwst'
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

    def Av_spec_SNR(self,wave,width=0,av_elements=2):
        #lam, snr = self.ETCresults['1d']['sn']
        lam, flux = self.ETCresults['1d']['extracted_flux']
        err = self.ETCresults['1d']['extracted_noise'][1]
        #plt.figure()
        #plt.axvspan(wave-width,wave+width,alpha=.3,color='orange')
        #plt.plot(lam,snr)

        low = wave - width
        high = wave + width
        ind = (lam >= low) & (lam <= high)
        flux = flux[ind]
        err = err[ind]
        ind = np.isfinite(flux) & np.isfinite(err)
        flux = flux[ind]
        err = err[ind]
        #snr = snr[ind]
        #snr = np.nansum(flux) / np.sqrt(np.nansum(err**2))
        #av = np.nanmean(snr)
        
        av_flux = np.nanmedian(flux)
        av_err = np.nanmedian(err)
        snr = np.sqrt(av_elements) * (av_flux/av_err)
        
        if np.isnan(snr): snr = 0
        return snr

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

        if self.instrument == 'nirspec':
            grating = self.pandeiacfg['configuration']['instrument']['disperser']

            if self.verbose: print('#############################\n#### Filter:%s, mag %f, grating: %s, for S/N=%.f \n#############################' % (filt, mag, grating, SNR))
        else:
            if self.verbose: print('#############################\n#### Filter %s, mag %.2f for S/N=%.f \n#############################' % (filt,mag,SNR))

        
        if (pd.isnull(mag)):
            print('#############################\n#### Filter %s, mag is not a number, returning NaNs\n#############################' % filt)
            results = {'plus':(None,None),
                       'minus':(None,None),
                       'best':(None,None)}  
            return(results)

        if self.verbose: print('#############################\n#### Filter %s, mag %.2f for S/N=%.f \n#############################' % (filt,mag,SNR))

           

        texp0=1000
        (SNR0, texp0) = self.Calculate_SNR(filt,mag,texp0,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)        
        if self.verbose>1: print('SNR=%6.2f for starting texp=%6.1f' % (SNR0,texp0))
        if SNR0==0.0 and mag<25:
            texp0=1 # get shortest exposure time
            (SNR0, texp0) = self.Calculate_SNR(filt,mag,texp0,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)        
            if self.verbose>1: print('SNR=%6.2f for starting texp=%6.1f' % (SNR0,texp0))
            if SNR0==0.0:
                #print('#############################\n#### Filter %s, mag %.2f is saturated, returning NaNs\n#############################' % filt)
                results = {'plus':(None,None),
                           'minus':(None,None),
                           'best':(None,None)}  
                return(results)
                
            
        
        
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
            
            
        (SNRnext,tnext) = self.Calculate_SNR(filt,mag,tnext,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)        
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

                (SNRnext,tnext) = self.Calculate_SNR(filt,mag,tnext,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)        
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
                    
                (SNRnext,tnext) = self.Calculate_SNR(filt,mag,tnext,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)
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
                                

    def Calculate_SNR_table(self, filters, magrange, exptime, lambkg4ETC=None,SNRformat=None,
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
                (SNRval,exptime1)=self.Calculate_SNR(filt,mag,exptime,lambkg4ETC=lambkg4ETC,spec=spec,**kwargs)
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

    


    def check_ref_lam(self,grating,lam):
        keys = list(grating_range.keys())
        grating = list(grating)
        g2 = []
        for i in range(len(grating)):
            ind = np.where(np.array(keys) == grating[i])[0]
            if (lam < grating_range[keys[ind]][0]) | (lam > grating_range[keys[ind]][1]):
                m = 'Reference wavelength is outside of wavelength range for {}'.format
                print(m)
            else:
                g2 += [grating[i]]
        return g2

    def get_grating_4_ref_lam(self,lam):
        keys = list(grating_range.keys())
        gratings = []
        for i in range(len(keys)):
            if (lam > grating_range[keys[i]][0]) & (lam < grating_range[keys[i]][1]):
                #m = ('Reference wavelength is outside of grating range, setting '+
                #     'new reference wavelength to be middle of grating')
                gratings += [keys[i]]
                #lam = (grating_range[keys[ind]][0] + grating_range[keys[ind]][1]) / 2 
        return gratings

    def Spec_texp_table(self,wave, reffilter, magrange, SNR,spec_av_width=.1, av_elements=2, gratings=None, lambkg4ETC=None,texp_type='best',
                           SNR_tolerance_in_percent=10,saveSNRflag=False,texpformat=None,SNRformat=None,
                           spec=None):
        """

        Parameters
        ----------
        wave : float/list
            Pass wavelengths to sample for the signal to noise.
        reffilter : string
            Pass filter for which to normalise the spectrum.
        magrange : list/tuple of magnitudes
            Pass a list of magnitudes for which to calculate the SNR..
        SNR : float 
            S/N  for which the exposure time should be calculated for.
        lambkg4ETC : list/tuple of two arrays, optional
            the first array is the wavelength, the second the background. 
            The default is None. If none, then the background data in 
            self.lambkg4ETC is used.
        spec_av_width : float, optional
            width in microns to average the spectrum around the central wavelength.
        av_element : int, optional 
            number of resolution elements to combine to calculate the SNR
        gratings : str or list, optional
            The grating or list of gratings to calculate the SNR with.
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
        spec : pysynphot spectrum
            spectrum to calculate the SNR from.

        Returns
        -------
        None. the table with the exposure times is saved as a pdastro object 
        in self.texp. The panda table is in self.texp.t. 
        The formatters for the table are in self.formatters4texptable
        Saving the table:
        self.texp.write('myfilename.txt',formatters=self.formatters4texptable)

        """
        self.spectrum=True
        #Make sure filters is a list and not a string
        
        g = gratings
        if isinstance(wave,float):wave=[wave]
        if isinstance(wave,int):wave=[wave]

        cols = [reffilter + ' mag']

        if SNRformat == None: SNRformat = self.SNRformat
        if texpformat == None: texpformat = self.texpformat

        self.texp = pdastroclass(columns=cols)
        self.texp.t[reffilter + ' mag']=magrange
        self.formatters4texptable = {}

        for w in wave:
            if self.instrument == 'nirspec':
                
                if g is None:
                    gratings = self.get_grating_4_ref_lam(w)
                else:
                    #gratings = self.check_ref_lam(g,w)
                    if isinstance(g,str):gratings=[g]
                    
            elif self.instrument == 'miri':
                gratings = ['p750l']
            

            for col in gratings: cols.append(str(w)+'micron_'+col+'_t')
            if saveSNRflag:
                for col in gratings: cols.append(str(w)+'micron_'+col+'_SN')

            
            
            for grat in gratings:
                
                self.formatters4texptable[str(w)+'micron_'+grat+'_t']=texpformat
                if saveSNRflag:
                    self.formatters4texptable[str(w)+'micron_'+grat+'_SN']=SNRformat
                if self.instrument == 'nirspec':
                    self.pandeiacfg['configuration']['instrument']['disperser'] = self.Check_grating(grat)
                    self.pandeiacfg['configuration']['instrument']['filter'] = self.set_grating_filter(grat)

                
                texps=[]
                SNRs=[]
                for mag in magrange:
                    texp_results=self.texp4SNRatmag(reffilter,mag,SNR,lambkg4ETC=lambkg4ETC,
                                                    SNR_tolerance_in_percent=SNR_tolerance_in_percent,
                                                    spec=spec,wave=w,width=spec_av_width,av_elements=av_elements)
                    texps.append(texp_results[texp_type][0])
                    if saveSNRflag:
                        SNRs.append(texp_results[texp_type][1])
                        
                self.texp.t[str(w)+'micron_'+grat+'_t']=np.array(texps)
                if saveSNRflag:
                    self.texp.t[str(w)+'micron_'+grat+'_SN']=np.array(SNRs)
                

if __name__ == '__main__':
    print('hello')
    jwst_SNR=jwst_SNRclass(instrument='nircam',mode='sw_imaging')
    print(jwst_SNR.aperture_radii)
    print(jwst_SNR.sky_annulii)
    
    jwst_SNR.verbose=1
    jwst_SNR.set_background4jwst(50,target='EmptyERS')

    #jwst_SNR.Calculate_SNR_table(['F200W'],np.arange(28.0,29.0,0.5),1200)
    #jwst_SNR.texp4SNRatmag('F200W',28.0,20.0)
    jwst_SNR.Imaging_texp_table(['F200W'],np.arange(28.0,29.0,0.5),10)
    
