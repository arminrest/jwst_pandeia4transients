#!/usr/bin/env python3
"""
Created on Fri Oct 16 18:44:27 2020

@author: arest

some definitions based on Dan Coe's jupyter notebooks:
    https://github.com/dancoe/pandeia-imaging/blob/master/depth/NIRCam%20Imaging%20Depth.ipynb

"""

import numpy as np
import math,sys,socket,os,re
from pdastro import pdastroclass
from astropy import coordinates,units
from jwst_backgrounds import jbt  # To calculate background

def RaInDeg(RA):
    s = re.compile('\:|\ ')
    if isinstance(RA,str) and s.search(RA):
        A = coordinates.Angle(RA,units.hour)
    else:
        A = coordinates.Angle(RA,units.degree)
    return(A.degree)
       
def DecInDeg(Dec):
    A = coordinates.Angle(Dec,units.degree)
    return(A.degree)


class background4jwstclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)        
        # lambda in microns
        self.lam = 4.5
        self.thresh = 1.1
        
        self.lam4percentile=4.5
        
        # This will contain the output of jbt=background4jwsts
        self.bkg = None
        
        self.defaulttargets = {
            "ElGordo":('01 02 55.2','-49 14 29.3'), # Low background example: El Gordo galaxy cluster ACT0102-49
            "EmptyERS":('03 32 42.397','-27 42 7.93') # well studied blank field in ERS
            }
        
        self.set_position_by_target('ElGordo')
        
        self.bgname = 'total_bg'
  
    def set_position(self,ra,dec,target):
        self.ra = RaInDeg(ra)
        self.dec = DecInDeg(dec)
        self.target = target
    
    def set_position_by_target(self,target):
        if not(target in self.defaulttargets):
            print('default positions:',self.defaulttargets)
            raise RuntimeError("%s could not be found in the default positions!" % target)
        self.set_position(self.defaulttargets[target][0],self.defaulttargets[target][1],target)
        
        
    def calc_background(self,lam=None,thresh=None):
        
        if not(lam is None):
            self.lam=lam
        if not(thresh is None):
            self.thresh=thresh
    
        self.bkg = jbt.background(self.ra, self.dec, wavelength=self.lam, thresh=self.thresh)
    
    def index4bkg_precentile(self,percentile,lam4percentile=None):
        if not(lam4percentile is None):
            self.lam4percentile=lam4percentile
        
        index4percentile = self.ix_equal('lam',self.lam4percentile)
        if len(index4percentile)!=1:
            print('allowed wavelengths:\n',list(self.t['lam']))
            raise RuntimeError("The specified wavelength of %f is not int the list of wavelenghts from the background: " % (self.lam4percentile))

        return(index4percentile[0])
    
    def bkgcolname(self,percentile):
        bkgcolname = '%s_%02d' % (self.bgname,round(percentile))
        return(bkgcolname)
        
    def get_background4percentiles(self,percentilelist,lam4percentile=None):
        """
        
        This routine saves the background fluxes (type of self.bgname, e.g., 'total_bg') 
        for a list of percentiles. The flux is saved in the table with column names self.bgname+'_%d' where %d 
        is the integer value of the percentile, e.g., 'total_bg_50' for 50 percentile 
        of total background flux
        
        Parameters
        ----------
        percentilelist : list or tuple
            list or tuple of percentiles.
        lam4percentile : float, optional
            wavelength in microns for which the percentiles are calculated. 
            The default is None. If None, then self.lam4percentile is used as default

        Returns
        -------
        None.

        """
        # make sure it is a list
        if isinstance(percentilelist,float) or isinstance(percentilelist,int):percentilelist=[percentilelist]
        
        self.t['lam']=self.bkg.bkg_data['wave_array']
       
        for percentile in percentilelist:
            # float values for percentile is silly accuracy...
            percentile=int(percentile)
            print('### Calculating percentile: %d' % percentile)
            
            # get the index for the wavelength value
            index = self.index4bkg_precentile(percentile,lam4percentile=lam4percentile)

            #print(self.t.at[index,'lam'])

            #print(self.bkg.bkg_data[self.bgname][0,index])
            #print(self.bkg.bkg_data[self.bgname][1,index])
            
            # bkg_data: first index is day, second is wavelength
            bkg_vs_day = self.bkg.bkg_data[self.bgname][:,index]
            # sort bkg_vs_day by indices
            ibkg_vs_day_sorted = bkg_vs_day.argsort()

            # get index of percentile
            nbkg = len(ibkg_vs_day_sorted)
            ibkg_percentile = round(percentile / 100. * nbkg)
            ibkg_percentile = np.clip(ibkg_percentile, 0, nbkg-1)
            
            # get the index of the day of the percentile
            ibkg_day = ibkg_vs_day_sorted[ibkg_percentile]
            
            # get the background for the day
            print('day: %d' % ibkg_day)
            bkg_vs_lam = self.bkg.bkg_data[self.bgname][ibkg_day]
            
            # Save the flux in the table with the percentile in the column name
            bkgcolname = self.bkgcolname(percentile)
            self.t[bkgcolname]=bkg_vs_lam
            
        return(0)
    
    def lambkg4ETC(self,percentile,lam=None,thresh=None,lam4percentile=None, target=None):
        """
        wrapper around all the calls to get the info (wavelength, bkg flux) for the ETC. 

        Parameters
        ----------
        percentile : TYPE
            DESCRIPTION.
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
        two numpy arrays: the first one are the wavelengths, and the second the bkg fluxes

        """
        # set the RA,Dec to the target if specified
        if not(target is None):
            self.set_position_by_target(target)
        
        # get the background from pandeia
        self.calc_background(lam=lam,thresh=thresh)

        # get teh background for the given percentile
        self.get_background4percentiles([percentile],lam4percentile=lam4percentile)

        # return the wavelength array and the background array        
        bkgcolname = self.bkgcolname(percentile)
        background = (np.array(self.t['lam']),np.array(self.t[bkgcolname]))
        return(background)
        

if __name__ == '__main__':
    print('hello')
    background4jwst=background4jwstclass()
    
    # set the position for which the background will be calculated. set_position_by_target uses 
    # pre-defined positions, set_position can be used to set any position
    background4jwst.set_position_by_target('ElGordo')
    
    # Get the background from pandeia
    background4jwst.calc_background()
    #print(background4jwst.bkg.bkg_data['total_bg'][0])
    #print(background4jwst.bkg.bkg_data['wave_array'])

    # calculate the backgrounds for the list of percentiles at the given wavelength
    background4jwst.get_background4percentiles([10.0,50.0,90.0],lam4percentile=4.5)
    # write the table to the screen
    background4jwst.write()
    
    # access the panda table directly ....
    print(background4jwst.t['lam'])
    print(background4jwst.t['total_bg_50'])
    
    # ... or call everything and get wavelength and background flux as a tuple of 2 np.arrays for the ETC input
    background = background4jwst.lambkg4ETC(50,target='EmptyERS')
    print(background)