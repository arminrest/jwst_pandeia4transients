#!/usr/bin/env python
'''
wrapper around pandas with convenience functions to ease handling of tables
A. Rest
'''
import sys,os,re,types
import numpy as np
from astropy.time import Time
import astropy.io.fits as fits
import astropy
import pandas as pd

def makepath(path,raiseError=True):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=True):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)
    
#https://numpy.org/doc/stable/reference/routines.set.html
def AorB(A,B):
    if len(A) == 0:
        return(B)
    if len(B) == 0:
        return(A)
    return(np.union1d(A,B))

def AandB(A,B,assume_unique=False):
    return(np.intersect1d(A,B,assume_unique=assume_unique))

def AnotB(A,B):
    return(np.setdiff1d(A,B))

def not_AandB(A,B):
    return(np.setxor1d(A,B))

class pdastroclass:
    def __init__(self,**kwargs):
        self.t = pd.DataFrame(**kwargs)
   
        self.verbose = 0
        
        # if self.auto_convert_dtypes==True, then before .to_string() is called in self.write, self.t.convert_dtypes() is run 
        # if self.default_dtypeMapping != None, then the dtype mapping is applied to the table .to_string() is called in self.write
        # This makes sure that formatters don't throw errors if the type of a column got changed to float or object during
        # one of the table operations
        self.auto_convert_dtypes = True
        self.default_dtypeMapping = None
        # example:
        # self.default_dtypeMapping = {'counter':np.int64}
        
        # default_formatters are passed to to_string() in self.write
        self.default_formatters = None
        # example:
        # self.default_formatters = {'MJD':'{:.6f}'.format,'counter':'{:05d}'.format}
       

    def load_spacesep(self,filename,test4commentedheader=True,namesMapping=None,roundingMapping=None,
                      na_values=['None','-','--'],**kwargs):
        
        kwargs['delim_whitespace']=True

        #also test for commented header to make it compatible to old format.
        self.load(filename,na_values=na_values,test4commentedheader=test4commentedheader,
                  namesMapping=namesMapping,roundingMapping=roundingMapping,**kwargs)

        return(0)

    def load(self,filename,raiseError=True,test4commentedheader=False,namesMapping=None,roundingMapping=None,**kwargs):
        #self.t = ascii.read(filename,format='commented_header',delimiter='\s',fill_values=[('-',0),('--',0)])
        try:
            self.t = pd.read_table(filename,**kwargs)
        except Exception as e:
            print('ERROR: could not read %s!' % filename)
            if raiseError:
                raise RuntimeError(str(e))
            return(1)
        
        if test4commentedheader:
            # This is to make it compatible to my old-style commented header files!
            # commented header: make sure it doesn't count the '#' as a column!
            if self.t.columns[0]=='#':
                renamemapping = {}
                for i in range(len(self.t.columns)-1):
                    renamemapping[self.t.columns[i]]=self.t.columns[i+1]
                renamemapping[self.t.columns[-1]]='__delme'
                self.t = self.t.rename(columns=renamemapping)
                self.t = self.t.drop(columns=['__delme'])

            # make sure the # is not kept in column name!
            if self.t.columns[0][0]=='#':
                self.t = self.t.rename(columns={self.t.columns[0]:self.t.columns[0][1:]})
            
        self.formattable(namesMapping=namesMapping,roundingMapping=roundingMapping)
        
        
        return(0)

    def write(self,filename=None,indices=None,columns=None,formatters=None,raiseError=True,overwrite=True,verbose=False, 
              index=False, makepathFlag=True,**kwargs):

        # make sure indices are converted into a valid list
        indices=self.getindices(indices)

        # make the path to the file if necessary
        if not (filename is None):
            if makepathFlag:
                if self.verbose: print('Clobbering %s' % filename)
                if makepath4file(filename,raiseError=False):
                    errorstring='ERROR: could not make directory for %s' % filename
                    if raiseError:
                        raise RuntimeError(errorstring)
                    #print(errorstring)
                    return(1)

            # if overwrite, then remove the old file first
            if os.path.isfile(filename):
                if overwrite:
                    os.remove(filename)
                    if os.path.isfile(filename):
                        errorstring='ERROR: could not clobber %s' % filename
                        if raiseError:
                            raise RuntimeError(errorstring)
                        #print(errorstring)
                        return(2)
                else:
                    print('Warning: file exists, not deleting it, skipping! if you want to overwrite, use overwrite option!')
                    return(0)
        
        # Fix the dtypes if wanted
        if self.auto_convert_dtypes:
            self.t=self.t.convert_dtypes()
        if not(self.default_dtypeMapping is None):
            self.formattable(dtypeMapping=self.default_dtypeMapping)            
        
        # if both formatters and self.defaultformatters are None, then no formatting. formatters supersedes self.defaultformatters
        if formatters is None:
            formatters = self.default_formatters
        
        if indices is None:
            # no indices are passed
            if self.verbose and not(filename is None): print('Saving %d rows into %s' % (len(self.t),filename))
            if len(self.t)==0:
                # just save the header
                if filename is None:
                    print(' '.join(columns)+'\n')
                else:
                    open(filename,'w').writelines(' '.join(columns)+'\n')
            else:
                if filename is None:
                    pass
                    print(self.t.to_string(index=index, columns=columns, formatters=formatters, **kwargs))
                else:
                    self.t.to_string(filename, index=index, columns=columns, formatters=formatters, **kwargs)
        else:
            if self.verbose and not(filename is None): print('Saving %d rows into %s' % (len(indices),filename))
            if len(indices)==0:
                # just save the header
                if filename is None:
                    print(' '.join(columns)+'\n')
                else:
                    if columns is None:
                        columns = []
                    open(filename,'w').writelines(' '.join(columns)+'\n')
            else:
                if filename is None:
                    print(self.t.loc[indices].to_string(index=index, columns=columns, formatters=formatters, **kwargs))
                else:
                    self.t.loc[indices].to_string(filename, index=index, columns=columns, formatters=formatters, **kwargs)

        if not (filename is None):
            # some extra error checking...
            if not os.path.isfile(filename):
                errorstring='ERROR: could not save %s' % filename
                if raiseError:
                    raise RuntimeError(errorstring)
                #print(errorstring)
                return(3)
                
        return(0)
        
    def formattable(self,namesMapping=None,roundingMapping=None,dtypeMapping=None):
             
        if not(namesMapping is None):
            self.t = self.t.rename(columns=namesMapping)
            
        if not(roundingMapping is None):
            self.t = self.t.round(roundingMapping)

        if not(dtypeMapping is None):
            for col in dtypeMapping:
                self.t[col] = self.t[col].astype(dtypeMapping[col])
            
        return(0)
    
    def getindices(self,indices=None):
        """make indices conform (input can be None,([list],), int, str, or list). The output is a list """
        
        #If indices is tuple, return the first entry which is the array list of indicies
        if indices is None: 
            return(self.t.index.values)
        
        # If indices=([indiceslist],), then it needs to be changed to indices=indiceslist
        if isinstance(indices,tuple):
            if len(indices)==0:
                #return an empty list instead of empty tuple
                return([])
            # teh first entry is a list, then these are the relevant indices!
            if isinstance(indices[0],list) or isinstance(indices[0],np.ndarray):
                return(indices[0])
            else:
                return(list(indices))
        
        # if the passed value is an integer or str, make it a list!
        if isinstance(indices,int) or isinstance(indices,str) or isinstance(indices,float):
            return([indices])           
        
        indices=np.array(indices)
        
        #if not (isinstance(indices,list) or isinstance(indices,np.ndarray)):
        #    raise RuntimeError("Can't convert this to an indices list!",type(indices),indices)
            
        return(indices)
    
    def getcolnames(self,colnames=None):
        """Return a list of all colnames of colnames=None. If colnames=string, return a list"""
        if (colnames is None) or colnames.lower()=='all':
            colnames = self.t.columns
        else:
            if isinstance(colnames,str):
                colnames=[colnames]
        return(colnames)
            

    def ix_remove_null(self,colnames=None,indices=None):
        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        for colname in colnames:
            #print('XXX',indices)
            (notnull,) = np.where(pd.notnull(self.t.loc[indices,colname]))
            indices = indices[notnull]
            #print('YYY',notnull)
        return(indices)

    def ix_inrange(self,colnames=None,lowlim=None,uplim=None,indices=None,
                   exclude_lowlim=False,exclude_uplim=False):

        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        for colname in colnames:
            if not(lowlim is None):
                if exclude_lowlim:
                    (keep,) = np.where(self.t.loc[indices,colname].gt(lowlim))
                else:
                    (keep,) = np.where(self.t.loc[indices,colname].ge(lowlim))
                indices = indices[keep]
                #print('lowlim cut:',keep)

            if not(uplim is None):
                if exclude_uplim:
                    (keep,) = np.where(self.t.loc[indices,colname].lt(uplim))
                else:
                    (keep,) = np.where(self.t.loc[indices,colname].le(uplim))
                indices = indices[keep]
                #print('uplim cut:',keep)
        return(indices)
    
    def ix_outrange(self,colnames=None,lowlim=None,uplim=None,indices=None,
                    exclude_lowlim=False,exclude_uplim=False):

        # get the indices based on input.
        indices=self.getindices(indices)        
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        #print('BBB',indices)
        for colname in colnames:
            if not(lowlim is None):
                if exclude_lowlim:
                    (keeplow,) = np.where(self.t.loc[indices,colname].lt(lowlim))
                else:
                    (keeplow,) = np.where(self.t.loc[indices,colname].le(lowlim))
                #print('lowlim cut:',keeplow)
            else:
                keeplow=[]

            if not(uplim is None):
                if exclude_uplim:
                    (keepup,) = np.where(self.t.loc[indices,colname].gt(uplim))
                else:
                    (keepup,) = np.where(self.t.loc[indices,colname].ge(uplim))
                #print('uplim cut:',keepup)
            else:
                keepup=[]
                
            indices = indices[AorB(keeplow,keepup)]
            
        return(indices)

        
    def fitsheader2table(self,fitsfilecolname,indices=None,requiredfitskeys=None,optionalfitskey=None,raiseError=True,skipcolname=None,headercol=None):

        indices = self.getindices(indices)        

        # initialize columns if necessary
        if requiredfitskeys!=None:
            for fitskey in requiredfitskeys:
                if not (fitskey in self.t.columns):
                    self.t[fitskey]=None
        if optionalfitskey!=None:
            for fitskey in optionalfitskey:
                if not (fitskey in self.t.columns):
                    self.t[fitskey]=None

        if headercol!=None and (not (headercol in self.t.columns)):
            self.t[headercol]=None

        # loop through the images
        for index in indices:
            header = fits.getheader(self.t[fitsfilecolname][index])
            if headercol!=None:
                self.t[headercol]=header
                
            if requiredfitskeys!=None:
                for fitskey in requiredfitskeys:
                    if fitskey in header:
                        self.t[fitskey][index]=header[fitskey]
                    else:
                        if raiseError:
                            raise RuntimeError("fits key %s does not exist in file %s" % (fitskey,self.t[fitsfilecolname][index]))
                        else:
                            self.t[fitskey][index]=None
                            if skipcolname!=None:
                                 self.t[skipcolname][index]=True
                                 
            if optionalfitskey!=None:
                for fitskey in optionalfitskey:
                    if fitskey in header:
                        self.t[fitskey][index]=header[fitskey]
                    else:
                        self.t[fitskey][index]=None

    def dateobs2mjd(self,dateobscol,mjdcol,timeobscol=None):
        if not (mjdcol in self.t.columns):
            self.t[mjdcol]=None

        if timeobscol!=None:
            dateobslist = list(self.t[dateobscol]+'T'+self.t[timeobscol])
        else:
            dateobslist = list(self.t[dateobscol])

        dateobjects = Time(dateobslist,  format='isot', scale='utc')
        mjds = dateobjects.mjd

        self.t[mjdcol]=mjds