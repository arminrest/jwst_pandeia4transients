import argparse
from astropy.io import ascii
from pandeia.engine.perform_calculation import perform_calculation
import json
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import minimize
import numpy as np
import warnings

import pysynphot as S # 
from extinction import fm07, fitzpatrick99, apply # https://pypi.org/project/extinction/

from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo 
from astropy.coordinates import Distance
from astropy import units as u
from background4jwst import background4jwstclass
from pandeia.engine.calc_utils import build_default_calc  # or alternatively, load ETC calculation parameters

package_directory = os.path.dirname(os.path.abspath(__file__))

class NIRSpec_SNR(object):
    """
    Wapper for the Pandeia engine to calculate the signal to noise and exposure times for 
    NIRSpec observations. Numerous options are available for readout mode, grating, filter,
    and observing mode.

    Important:
        ref_wave must be provided, wavelength the SN is calculated at
        A spectrum must be provided.
        The distance is measured in Mpc

    Note:
        The spectrum input units should be microns for wavelength and mJy for flux.
        The normalising magnitude is as measured by NIRCam with mode = sw_imaging filter = f115w
    
    Key functions:
        Calculate_SNR()
            Calculates the expected signal to noise for a given configuration and exposuretime.
            The configuration is automatically determined from the input exposure time 
        
        Calculate_exposure_for_SN()
            Calculates the exposure time needed for the desired signal to noise. 
            The instrument configuration is set to the values that achieve the shortest exposure
            time for the same signal to noise.

    Missing components:
        Input spectra can be normalised by a magnitude in a given NIRCam configuration.
        The default is set to "nircam,sw_imaging,f115w" and for now can only be changed by editing the 
        template "nirspec_fs_default.json". 

        It needs a more efficient way to calculate the optimal setup in Calculate_exposure_for_SN(), 
        currently it loops through group values. Scipy minimize didn't work so well here.
        
        Different NIRSpec observing modes. Currently only long_slit is implemented, need to find 
        documentation on the other modes and their calls in Pandeia.

        Commandline interface?

        TBD


    """
    def __init__(self,ref_wave=None,av=None,dist=None,
                 z=None,
                 grating='prism',filt=None,
                 read_mode='nrsirs2rapid',mode='fixed_slit',
                 spec = None,exptime=None):

        self.grating = grating.lower() # grating name
        self.Check_grating()
        if filt is None:
            self.Assign_filter() # filter name
        else:
            self.filter = filt # filter name
        self.ngroups = None # number of frame groups
        self.nint = 1 # Number of integrations
        self.nexp = 4 # Number of exposures
        self.mode = mode # Spectrograph observation mode
        self.read_mode = read_mode # readout mode
        self.exp_time = exptime # Desired exposure time
        self.pandeiacfg = None 
        self.mag = None
        self.ref_filter = 'f200w'
        self.ref_sn = None # signal to nosie 

        if ref_wave is not None:
            self.ref_wave = ref_wave # reference wavelength in mirons
        else:
            self.ref_wave = 2 # reference wavelength in miron
        self.ref_filter
        self.av = av # v band extinction
        self.dist = dist
        self.z = z

        self.spec = spec
        self.background4jwst = background4jwstclass()
        self.lambkg4ETC = None
        # calculated by pandeia
        self.calc_flux = None # units e-/s
        self.calc_wave = None # units micron
        self.calc_err = None  # units e-/s
        self.calc_sn = None   # singal to noise of full spec
        self.calc_exp = None
    
    def Assign_filter(self):
        grating = self.grating
        if (grating == 'g140h') | (grating == 'g140m'):
            message = ('2 filters available: f070lp, or f100lf \n' +
                        'Assigning f070lp')
            warnings.warn(message)
            filt  = 'f070lp'

        elif (grating == 'g235h') | (grating == 'g235m'):
            filt = 'f170lp'
        elif (grating == 'g395h') | (grating == 'g395m'):
            filt = 'f290lp'
        elif grating == 'prism':
            filt = 'clear'
        self.filter = filt



    def Check_grating(self):
        """
        Check to see if the input grating value is a valid option
        """
        allowed = ['prism', 'g140h', 'g140m', 
                   'g235h', 'g235m', 'g395h', 'g395m']
        for i in range(len(allowed)):
            if self.grating == allowed[i]:
                return 
        message = ('No such grating available, please choose from:\n '
                   + 'prism\n g140h\n g140m\n g235h\n g235m\n g395h\n g395m')
        raise(ValueError(message))
        
    def Check_filter(self):
        """
        Check to see if the input filter value is a valid option
        """
        allowed = ['clear', 'f070lp', 'f100lp', 
                   'f110w', 'f140x', 'f170lp', 'f290lp']
        for i in range(len(allowed)):
            if self.filter == allowed[i]:
                return 
        message = ('No such filter available, please choose from:\n '
                   + 'clear\n f070lp\n f100lp\n f110w\n f140x\n f170lp\n f290lp')
        raise(ValueError(message))
    
    def Check_read_mode(self):
        """
        Check to see if the input read mode value is a valid option
        """
        allowed = ['nrsrapid', 'nrsrapidd6','nrs',
                   'nrsirs2rapid','nrsirs2']
        for i in range(len(allowed)):
            if self.read_mode == allowed[i]:
                return 
        message = ('No such readout mode available, please choose from:\n '
                   + 'nrsrapid\n nrsrapidd6\n nrs\n nrsirs2rapid\n nrsirs2')
        raise(ValueError(message))
    
    def Check_mode(self):
        """
        Check to see if the input mode value is a valid option
        """
        allowed = ['fixed_slit','ifu']
        for i in range(len(allowed)):
            if self.mode == allowed[i]:
                return 
        message = ('No such mode available, please choose from:\n '
                   + 'fixed_slit\n ifu\n')
        raise(ValueError(message))

    def Check_exposure_time(self):
        """
        Check to see if the input exposure time value is a valid option
        """
        if self.exp_time is None:
            message = "Exposure time can't be None"
            raise(ValueError(message))
        #if (type(self.exp_time) != float) & (type(self.exp_time) !=int):
        #   print(type(self.exp_time))
        #   message = "Exposure time must be either a float or int"
        #   raise ValueError(message)
        return
    
    def NIRSpec_exp_2_groups(self):
        """
        Calculate the number of groups that can be used recorded for a 
        given exposure time and configuration.
        """
        self.Check_read_mode()
        self.Check_exposure_time()
        mode = self.read_mode
        exp = self.exp_time
        if mode == 'nrsrapid':
            ng = (10.737)
        elif mode == 'nrsrapidd6':
            ng = (75.159)
        elif mode == 'nrs':
            ng = (42.947)
        elif mode == 'nrsirs2rapid':
            ng = (14.589)
        elif mode == 'nrsirs2':
            ng = (72.944)
        ng = np.floor(exp / (self.nint * self.nexp * ng + 58.35 * self.nint))
        if ng > 100:
            self.nint += 1
            ng = np.floor(exp / (self.nint * self.nexp * ng + 58.35 * self.nint))
        if ng <1:
            warnings.warn('exposure is too short for even one group!'+
                        'setting ng = 1') 
            ng =1 
        self.ngroups = ng
        return
    
    def NIRSpec_groups_2_exp(self):
        """
        Calculate the total exposure time for a given configuration.
        """
        self.Check_read_mode()
        self.Check_exposure_time()
        mode = self.read_mode
        if mode == 'nrsrapid':
            exp = 10.737 
        elif mode == 'nrsrapidd6':
            exp = 75.159 
        elif mode == 'nrs':
            exp = 42.947 
        elif mode == 'nrsirs2rapid':
            exp = 14.589
        elif mode == 'nrsirs2':
            exp = 72.944 
        exp = exp * self.nint * self.nexp * self.ngroups * 14.6

        self.exp_time = exp
        return
    
    def Check_all_fields(self):
        """
        Check if all parameters have been filled in correctly to enter into Pandeia
        """
        self.Check_grating()
        self.Check_filter()
        self.Check_read_mode()
        self.Check_mode()
        message = ''
        if self.ngroups is None:
            message += 'Number of groups must be defined\n'
        elif self.ngroups < 4:
            m = ('Number of groups is only {}, a value >4 is prefered,'.format(self.ngroups) +
                ' so increasing exposure time is advised.')
            warnings.warn(m)
        if self.nint is None:
            message += 'Number of integrations must be defined\n'
        if self.nexp is None:
            messae += 'Number of exposures must be defined\n'
        if self.read_mode is None:
            message += 'Read mode must be defined\n'
        if self.spec is None:
            message += 'No spectrum provided\n'
        if self.ref_wave is None:
            message += 'No reference wavelength specified (ref_wave = None)\n'
        if message != '':
            raise(ValueError(message))
        return
        
    def Get_background(self):
        """
        Add the background model to the pandeia table.
        The fits file was provided by Klausvia Pantoppindan via Nora Luetzgendorf
        """
        file = os.path.join(package_directory,'data/minzodi12_12052016.fits')
        hdu = fits.open(file)
        table = hdu[1].data
        back_lam = table.field('wavelength')
        back_all =  (table.field('background') + table.field('thermal') + 
                     table.field('straylight') + table.field('infield'))
        self.pandeiacfg['background'] = [list(back_lam),list(back_all)]
        return
    
    def Normalise_spec(self):
        imgr_data = self.imgr_data
        imgr_data['scene'][0]['spectrum']["normalization"] = {}
        imgr_data['scene'][0]['spectrum']["normalization"]["bandpass"]= "nircam,lw_imaging," + self.ref_filt
        imgr_data['scene'][0]['spectrum']["normalization"]["norm_flux"] = self.mag
        imgr_data['scene'][0]['spectrum']["normalization"]["norm_fluxunit"] =  "abmag"
        imgr_data['scene'][0]['spectrum']["normalization"]["type"] = "jwst"
        self.imgr_data = imgr_data

    def Make_config(self):
        """
        Put the configuration data into the format expected by Pandeia, 
        using the default template provided by Nora Luetzgendorf.
        """
        if self.ngroups is None:
            self.NIRSpec_exp_2_groups()

        self.Check_all_fields()
        if self.pandeiacfg is None:
            self.pandeiacfg = build_default_calc('jwst','nirspec','fixed_slit')
        

        self.pandeiacfg['configuration']['detector']['ngroup'] = self.ngroups
        self.pandeiacfg['configuration']['detector']['nint'] = self.nint
        self.pandeiacfg['configuration']['detector']['nexp'] = self.nexp
        self.pandeiacfg['configuration']['instrument']['mode'] = self.mode
        self.pandeiacfg['configuration']['instrument']['filter'] = self.filter
        self.pandeiacfg['configuration']['instrument']['disperser'] = self.grating
        self.pandeiacfg['configuration']['detector']['readout_pattern'] = self.read_mode
        self.pandeiacfg['configuration']['detector']['subarray'] = 'full'
        
        self.spec.convert('micron')
        self.spec.convert('mJy')
        self.pandeiacfg['scene'][0]['spectrum']['sed']['sed_type'] = 'input'
        self.pandeiacfg['scene'][0]['spectrum']['sed']['spectrum'] = [self.spec.wave,self.spec.flux]
        self.pandeiacfg['scene'][0]['spectrum']['sed']['unit'] = 'flam'

        if self.mag is not None:
            self.pandeiacfg['scene'][0]['spectrum']['normalization']['type'] = 'jwst'
            self.pandeiacfg['scene'][0]['spectrum']['normalization']['bandpass'] = 'nircam,sw_imaging,' + self.ref_filter
            self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_flux'] = self.mag
            self.pandeiacfg['scene'][0]['spectrum']['normalization']['norm_fluxunit'] = 'abmag'
        else:
            self.pandeiacfg['scene'][0]['spectrum']['normalization']['type'] = 'none'

                            
        #imgr_data['scene'][0]['spectrum']['normalization']['norm_flux'] = self.mag
        self.pandeiacfg['background'] = self.lambkg4ETC
        return
    
    def Pattern_select(self):
        """
        Arbitrary read mode switch based on exposure time.
        This needs to be refined.
        """
        if self.exp_time < 2000:
            self.read_mode = 'nrsirs2rapid'
        else:
            self.read_mode = 'nrsirs2'
        return
    
    def SN_for_wavelength(self,wave=None):
        if wave is not None:
            self.ref_wave = wave
        ind = np.argmin(abs(self.calc_wave - self.ref_wave))
        if abs(self.calc_wave - self.ref_wave) > .1:
            raise(ValueError('reference wavelength not in wavelength range'))
        return self.calc_sn[ind]


    def Max_SN(self):
        return np.nanmax(self.calc_sn)


    def Calculate_SNR(self,exptime=None,ngroups= None,target=None,bkgpercentile=50):
        """
        Calculate the signal to noise of the provided spectrum under the 
        given instrument configuration.

        Inputs
        ------
            exptime : float 
                total exposure time, degenerate with ngroup
            ngroups : int
                total number of groups, degenerate with exptime

        Outputs
        -------
            self.ref_sn : float
                Calculated signal to noise 
            real_t : float
                total exposure time as calculated by Pandeia 
                for the given configuration

        """
        if (exptime is not None) & (ngroups is None):
            self.exp_time = exptime
            self.ngroups = None
        elif (exptime is None) & (ngroups is not None):
            self.ngroups = ngroups
        
        self.set_background4jwst(bkgpercentile,target='EmptyERS')
        self.Make_config()
        #rint('calculate sn ng ',self.ngroups)
        results = perform_calculation(self.pandeiacfg)   # results is a dictionary
        self.calc_flux = results['1d']['extracted_flux'][1]
        self.calc_wave = results['1d']['extracted_flux'][0]
        self.calc_err = results['1d']['extracted_noise'][1]
        self.calc_sn = results['1d']['sn'][1] # S/N for full spectrum
        self.ref_sn = results['scalar']['sn'] # S/N for reference wavelength
        self.calc_exp = results['scalar']['total_exposure_time']
        return 
    
    
    def SN_minimise(self,groups,target_sn):
        """
        Not working
        """
        #print('groups ',groups)
        self.ngroups = groups
        self.Calculate_SNR()
        SN = self.SN_for_wavelength()
        res = abs(SN-target_sn)
        if res <= 0.2:
            res = 0

        return res
    
    def Calculate_exposure_for_SN(self,SN,verbose=False):
        """
        Calculates the exposure time needed and the optimal read out mode 
        to reach a given signal to noise.

        Inputs
        ------
            SN : float
                signal to noise

        Outputs
        -------
            exp_times : float
                the required exposure time to reach the desired SN
        """
        methods = ['nrsirs2rapid','nrsirs2']
        traditional = ['nrsrapid', 'nrsrapidd6','nrs']
        exp_times = []
        self.exp_time = None
        g = []
        #bds = [(100,3000)]
        for m in methods:
            self.read_mode = m
            #res = minimize(self.SN_minimise,g0,args=(SN),options={'eps':1},
            #bounds=[(0,40)])#,method='Nelder-Mead')
            groups, t = self.groups_for_SN(target=SN)
            g += [groups]
            exp_times += [t]
        exp_times = np.array(exp_times)
        m = np.argmin(exp_times)
        self.ngroups = g[m]
        self.read_mode = methods[m]
        
        if verbose:
            print('Optimal mode: ' + methods[m])
            print('exposure time: ' + str(exp_times[m]))
        return exp_times[m]
    
    def groups_for_SN(self,target):
        """
        Clunky loop to get the appropirate exposure time for a target SN

        Inputs 
        ------
            target : float
                target signal to noise

        Outputs
        -------
            g : int
                number of groups 
            t : float 
                exposure time to reach the desired SN
        """
        g = 4
        lim = 0.1
        maxiter = 34
        diff = target
        while (diff > lim) & (g < maxiter):
            self.ngroups = g
            self.Calculate_SNR()
            SN = self.SN_for_wavelength()
            diff = target - sn
            g += 1
        return g, self.calc_exp

    def Redden_spec(self,av=None):
        """
        Redden the input spectrum according to the fitzpatrick 99 relation and 
        the visible band extinction.

        Inputs
        ------
            av : float
                v band extinction that is applied via teh Fitzpatric 99 model.
        """
        if av is not None:
            self.av = av

        if self.av is not None:
            self.spec.convert('Angstrom')
            wav = self.spec.wave # convert to Angstrom for extinction
            red = apply(fitzpatrick99(wav.astype('double'),self.av,3.1),self.spec.flux)
            self.spec = S.ArraySpectrum(self.spec.wave,red,waveunits=self.spec.waveunits,
                                        fluxunits=self.spec.fluxunits)
        return

    def Distance_modulus(self,apparent,absolute):
        d = 10**((apparent - absolute + 5)/5) /1e6
        self.dist = d

    def Distance_scale(self,dist=None,z=None,apparent=None):
        """
        Scales and redshifts the spectrum to the given distance/redshift.

        Inputs
        ------
            dist : float
                distance to the source in Mpc
            z : float 
                redshift to the source
        """
        if (dist is not None):
            self.dist = dist
            self.z = None
        if z is not None:
            self.z = z
            self.dist = None

        #if apparent is not None:
        #    self.mag = apparent

        #self.Distance_modulus(self.mag,-19)

        dt = self.dist is not None
        zt = self.z is not None
        if dt & zt:
            warnings.warn('Both distance and redshift set, defaulting to redshift.')
            dt = False
        if dt & ~zt:
            d = Distance(self.dist,unit=u.Mpc)
            self.dist = d
            try:
                self.z = d.compute_z(cosmo)
            except:
                self.z = 0 


        elif zt & ~dt:
            self.dist = cosmo.luminosity_distance(self.z).to(u.Mpc)

        if self.dist > 0:
            flux = self.spec.flux * (10/(self.dist.to('pc').value))**2 # original spectral distance 10pc
            
            spec = S.ArraySpectrum(self.spec.wave, flux, 
                                    waveunits=self.spec.waveunits,fluxunits=self.spec.fluxunits)
            spec = spec.redshift(self.z)

            spec.convert('micron')
            spec.convert('flam')
            
            self.spec = spec            
        return



    def Random_realisation(self):
        flux = self.calc_flux
        wave = self.calc_wave
        error = self.calc_err

        noise = np.random.normal(0,error)

        random_noise = np.zeros((2,len(wave)))
        random_noise[0,:] = wave
        random_noise[1,:] = flux + noise

        return random_noise


    def initialize_pandeia(self,instrument,mode,ETCjsonfile=None):
        print('Initializing pandeia with %s, %s' % (instrument,mode))
        instrument = instrument.lower()
        mode = mode.lower()
        if not (instrument in self.allowed_instruments):
            raise(RuntimeError('instrument %s not in %s' % (instrument,' '.join(self.allowed_instruments))))
        self.readoutpattern=readoutpatternclass(instrument)
        if ETCjsonfile is None:
            self.pandeiacfg=build_default_calc('jwst',instrument,mode)
        else:
            with open(ETCjsonfile) as f:  # use a json file you have
                 self.pandeiacfg = json.load(f)
        return(0)

    def Simulate_spec(self,exp,dist=None,redshift=None,av=None,mag=None,norm_filt=None):
        
        self.exp_time = exp
        self.dist=dist
        self.z = redshift
        self.av = av
        self.mag = mag
        self.norm_filt = norm_filt

        if (dist is None) & (redshift is None):
            warnings.warn('No distance specified')
        else:
            self.Distance_scale()
        if av is not None:
            self.Redden_spec()

        self.Calculate_SNR()

        rand = self.Random_realisation()
        return rand

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