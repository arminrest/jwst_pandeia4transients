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
import astropy.io.fits as pyfits  # To load background

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
        # percentiles used by ETC web interface: Low=10, Medium=50, High=90
        self.percentile = 50.0
        
        # lambda in microns
        self.lam = 4.5
        self.thresh = 1.1
        
        self.lam4percentile=4.5
        
        # This will contain the output of jbt=background4jwsts
        self.bkg = None
        
        self.defaultpositions = {
            "ElGordo":('01 02 55.2','-49 14 29.3'), # Low background example: El Gordo galaxy cluster ACT0102-49
            "EmptyERS":('03 32 42.397','-27 42 7.93') # well studied blank field in ERS
            }
        
        self.set_position_by_name('ElGordo')
        
        self.bgname = 'total_bg'
  
    def set_position(self,ra,dec,name):
        self.ra = RaInDeg(ra)
        self.dec = DecInDeg(dec)
        self.target = name
    
    def set_position_by_name(self,name):
        if not(name in self.defaultpositions):
            print('default positions:',self.defaultpositions)
            raise RuntimeError("%s could not be found in the default positions!" % name)
        self.set_position(self.defaultpositions[name][0],self.defaultpositions[name][1],name)
        
        
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
            wavelength for which the percentiles are calculated. The default is None.

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
    
    def get_lam_bkg4ETC(self,percentile):
        bkgcolname = self.bkgcolname(percentile)
        return(np.array(background4jwst.t['lam']),np.array(background4jwst.t[bkgcolname]))

if __name__ == '__main__':
    print('hello')
    background4jwst=background4jwstclass()
    
    # set the position for which the background will be calculated. set_position_by_name uses 
    # pre-defined positions, set_position can be used to set any position
    background4jwst.set_position_by_name('ElGordo')
    
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
    
    # ... or get it as 2 arrays
    print(background4jwst.get_lam_bkg4ETC(50))