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

package_directory = os.path.dirname(os.path.abspath(__file__))

class NIRSpec_SN(object):
    """
    Wapper for the Pandeia engine to calculate the signal to noise and exposure times for 
    NIRSpec observations. Numerous options are available for readout mode, grating, filter,
    and observing mode.

    Important:
        ref_wave must be provided, wavelength the SN is calculated at
        a spectrum must be provided.

    Note:
        The spectrum input units should be microns for wavelength and mJy for flux.
        The normalising magnitude is as measured by NIRCam with mode = sw_imaging filter = f115w
    
    Key functions:
        Calculate_SN()
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
                 grating='prism',filt='clear',
                 read_mode='nrsirs2rapid',mode='fixed_slit',
                 spectrum = None,exptime=None):
        self.grating = grating.lower() # grating name
        self.Check_grating()
        self.filter = filt # filter name
        self.ngroups = None # number of frame groups
        self.nint = 1 # Number of integrations
        self.nexp = 4 # Number of exposures
        self.mode = mode # Spectrograph observation mode
        self.read_mode = read_mode # readout mode
        self.exp_time = exptime # Desired exposure time
        self.imgr_data = None 
        self.mag = None
        self.SN = None # signal to nosie 
        self.ref_wave = ref_wave # reference wavelength in mirons
        self.av = av # v band extinction
        self.dist = dist
        self.z = z

        if spectrum is not None:
            self.wave = spectrum[0,:]
            self.flux = spectrum[1,:]
        else:
            self.wave = None
            self.flux = None
        
    def Check_grating(self):
        """
        Check to see if the input grating value is a valid option
        """
        allowed = ['prism', 'g140h', 'g140m', 
                   'g235h', 'g235m', 'g395h', 'g395m']
        for i in range(len(allowed)):
            if self.grating == allowed[i]:
                return 
        message = ('No such grating available, please choose from:\n'
                   + 'prism\n g140h\n g140m\n g235h\n g235m\n g395h\n g395m')
        raise(ValueError, mesage)
        
    def Check_filter(self):
        """
        Check to see if the input filter value is a valid option
        """
        allowed = ['clear', 'f070lp', 'f100lp', 
                   'f110w', 'f140x', 'f170lp', 'f290lp']
        for i in range(len(allowed)):
            if self.filter == allowed[i]:
                return 
        message = ('No such filter available, please choose from:\n'
                   + 'clear\n f070lp\n f100lp\n f110w\n f140x\n f170lp\n f290lp')
        raise(ValueError, mesage)
    
    def Check_read_mode(self):
        """
        Check to see if the input read mode value is a valid option
        """
        allowed = ['nrsrapid', 'nrsrapidd6','nrs',
                   'nrsirs2rapid','nrsirs2']
        for i in range(len(allowed)):
            if self.read_mode == allowed[i]:
                return 
        message = ('No such readout mode available, please choose from:\n'
                   + 'nrsrapid\n nrsrapidd6\n nrs\n nrsirs2rapid\n nrsirs2')
        raise(ValueError, mesage)
    
    def Check_mode(self):
        """
        Check to see if the input mode value is a valid option
        """
        allowed = ['fixed_slit']
        for i in range(len(allowed)):
            if self.mode == allowed[i]:
                return 
        message = ('No such mode available, please choose from:\n'
                   + 'fixed_slit\n ')
        raise(ValueError, mesage)

    def Check_exposure_time(self):
        """
        Check to see if the input exposure time value is a valid option
        """
        if self.exp_time is None:
            message = "Exposure time can't be None"
            raise(ValueError, mesage)
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
            ng = exp / (10.737)
        elif mode == 'nrsrapidd6':
            ng = exp / (75.159)
        elif mode == 'nrs':
            ng = exp / (42.947)
        elif mode == 'nrsirs2rapid':
            ng = exp / (14.589)
        elif mode == 'nrsirs2':
            ng = exp / (72.944)
        ng = np.floor(ng / (self.nint * self.nexp))
        if ng <1:
            warnings.warn('exposure is too short for even one group!')
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
        exp = exp * self.nint * self.nexp * self.ngroups

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
        if self.wave is None:
            message += 'No spectrum provided\n'
        if self.ref_wave is None:
            message += 'No reference wavelength specified (ref_wave = None)\n'
        if message != '':
            raise Warning(message)
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
        self.imgr_data['background'] = [list(back_lam),list(back_all)]
        return
    
    def Make_config(self):
        """
        Put the configuration data into the format expected by Pandeia, 
        using the default template provided by Nora Luetzgendorf.
        """
        if self.ngroups is None:
            self.NIRSpec_exp_2_groups()

        self.Check_all_fields()
        jf = os.path.join(package_directory,'data/nirspec_fs_default.json')
        with open(jf) as f:
            imgr_data = json.load(f)   # this is a python dictionary

        imgr_data['configuration']['detector']['ngroup'] = self.ngroups
        imgr_data['configuration']['detector']['nint'] = self.nint
        imgr_data['configuration']['detector']['nexp'] = self.nexp
        imgr_data['configuration']['instrument']['mode'] = self.mode
        imgr_data['configuration']['instrument']['filter'] = self.filter
        imgr_data['configuration']['instrument']['disperser'] = self.grating
        imgr_data['configuration']['detector']['readout_pattern'] = self.read_mode
        
        imgr_data['scene'][0]['spectrum']['sed']['spectrum'] = [self.wave, self.flux]
        imgr_data['strategy']['reference_wavelength'] = self.ref_wave
        if self.mag is not None:
            imgr_data['scene'][0]['spectrum']['normalization']['norm_flux'] = self.mag
        
        self.imgr_data = imgr_data
        self.Get_background()
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
    
    def Calculate_SN(self,exptime=None,ngroup= None):
        """
        Calculate the signal to noise of the provided spectrum under the 
        given instrument configuration.

        Inputs
        ------
            exptime : float 
                total exposure time, degenerate with ngroup
            ngroup : int
                total number of groups, degenerate with exptime

        Outputs
        -------
            self.SN : float
                Calculated signal to noise 
            real_t : float
                total exposure time as calculated by Pandeia 
                for the given configuration

        """
        if (exptime is not None) & (ngroup is None):
            self.exp_time = exptime
            self.ngroup = None
        elif (exptime is None) & (ngroup is not None):
            self.ngroups = ngroup
        #else:
        #   raise ValueError('define either exptime OR ngroup')
        #if (self.exp_time is not None) :
            #self.NIRSpec_exp_2_groups
            #self.Pattern_select()
            #print(self.exptime)
            #self.NIRSpec_exp_2_groups()
        #print(self.exp_time)
        self.Make_config()
        #rint('calculate sn ng ',self.ngroups)
        results = perform_calculation(self.imgr_data)   # results is a dictionary
        self.SN = results['scalar']['sn']
        real_t = results['scalar']['total_exposure_time']
        #if real_t> 
        #print('exp',real_t)
        #print(self.SN)
        return self.SN, real_t
    
    
    def SN_minimise(self,groups,target_sn):
        """
        Not working
        """
        #print('groups ',groups)
        self.ngroups = groups
        SN,_ = self.Calculate_SN()
        #print('sncalc ',SN)
        #print(abs(SN-target_sn))
        res = abs(SN-target_sn)
        if res <= 0.2:
            res = 0
        #print('residual ',res)
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
            sn,t = self.Calculate_SN()
            diff = target - sn
            g += 1
        return g, t

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
            wav = self.wave * 1e4 # convert to Angstrom for extinction
            red = apply(fitzpatrick99(wav.astype('double'),self.av,3.1),self.flux)
            self.flux = red
        return

    def Distance_scale(self,dist=None,z=None):
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

        dt = self.dist is not None
        zt = self.z is not None
        if dt & zt:
            warnings.warn('Both distance and redshift set, defaulting to redshift.')
            dt = False
        if dt & ~zt:
            d = Distance(self.dist,unit=u.Mpc)
            try:
                self.z = d.compute_z(cosmo)
            except:
                self.z = 0 

        elif zt & ~dt:
            self.dist = cosmo.luminosity_distance(self.z).to(u.Mpc).value

        if self.dist > 0:
            self.flux = self.flux * (10/(self.dist*1e6))**2 # original spectral distance 10pc


            spec = S.ArraySpectrum(self.wave, self.flux, waveunits='micron',fluxunits='mjy')
            spec = spec.redshift(self.z)

            spec.convert('micron')
            spec.convert('mjy')

            self.wave = spec.wave
            self.flux = spec.flux

        return





