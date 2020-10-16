#!/usr/bin/env python3
"""
Created on Fri Oct 16 18:44:27 2020

@author: arest
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
        
    def get_background4percentiles(self,percentilelist,lam4percentile=None):
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
            bkgcolname = '%s_%02d' % (self.bgname,percentile)
            self.t[bkgcolname]=bkg_vs_lam
            
        return(0)


if __name__ == '__main__':
    print('hello')
    background4jwst=background4jwstclass()
    
    background4jwst.set_position_by_name('ElGordo')
    
    background4jwst.calc_background()
    #print(background4jwst.bkg.bkg_data['total_bg'][0])
    #print(background4jwst.bkg.bkg_data['wave_array'])

    background4jwst.get_background4percentiles([10.0,50.0,90.0],lam4percentile=4.5)
    background4jwst.write()