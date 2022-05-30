#!/usr/bin/env python3
# coding: utf-8

## @file
## @brief read binary output of ASWNS, compute some quantities, and generate gnuplot-friendly text file
## @details units are c=G=M☉=kB=1
## @author Giovanni Camelio
## @date first version: 2020-11-17, this version: 2022-05-17

# I repeat the documentation in Doxygen format and in python docstring format

import numpy as np

## import the stellar structure from a file in the binary format and perform some operations on it
class binstar:
    """
    import the stellar structure from a file in the binary format and perform
    some operations on it.
    
    BEWARE: all units are c=G=M☉=kB=1.
    
    EXAMPLE:
    > star= binstar("path/to/file.out")
    > star.nth # number of angular points
    """
    
    ## load binary file
    def __init__(self,filename):
        # open the binary ('b') file read only ('r')
        fileid= open(filename,'rb')
        
        ## intger, number of angular grid points
        self.nth= np.fromfile(fileid,dtype='int32',  count=1)[0]
        
        ## integer, number of radial grid points in the whole domain
        self.nr= np.fromfile(fileid,dtype='int32',  count=1)[0]
        # the domain is (nth x nr)
        
        ## integer, index of the equator in the angular direction
        self.ieqth= np.fromfile(fileid,dtype='int32',  count=1)[0] - 1
        # note that planar symmetry <==> ieqth == nth - 1
        # I subtract 1 because of the different index convection between Fortran and Python
        
        # integer (temp variable): number of radial grid points for the matter quantities
        mr= np.fromfile(fileid,dtype='int32',  count=1)[0]
        # note that mr= max(surf) < nr
        
        ## double[nr], radial grid
        self.r= np.fromfile(fileid,dtype='float64',count=self.nr)
        
        ## double[nth], angular grid
        self.th= np.fromfile(fileid,dtype='float64',count=self.nth)
        
        ## integer[nth], position of the surface
        ## @details [   th[i], r[ surf[i] ]   ] is the first point outside the star
        self.surf= np.fromfile(fileid,dtype='int32',  count=self.nth)
        
        ## double[nth,nr], rest mass density
        self.rho= np.zeros(dtype='float64',shape=(self.nth,self.nr))
        self.rho[:,:mr]= np.fromfile(fileid,dtype='float64',count=self.nth*mr).reshape(mr,self.nth).transpose()
        
        ## double[nth,nr], pressure
        self.p= np.zeros(dtype='float64',shape=(self.nth,self.nr))
        self.p[:,:mr]= np.fromfile(fileid,dtype='float64',count=self.nth*mr).reshape(mr,self.nth).transpose()
        
        ## double[nth,nr], enthalpy density
        self.hden= np.zeros(dtype='float64',shape=(self.nth,self.nr))
        self.hden[:,:mr]= np.fromfile(fileid,dtype='float64',count=self.nth*mr).reshape(mr,self.nth).transpose()
        
        ## double vphi[nth,nr]: contravariant velocity along phi
        self.vphi= np.zeros(dtype='float64',shape=(self.nth,self.nr))
        self.vphi[:,:mr]= np.fromfile(fileid,dtype='float64',count=self.nth*mr).reshape(mr,self.nth).transpose()
        
        ## double[nth,nr], spacetime drag = -contravariant shift along phi
        self.omg= np.fromfile(fileid,dtype='float64',count=self.nth*self.nr).reshape(self.nr,self.nth).transpose()
        
        ## double[nth,nr], conformal factor
        self.psi= np.fromfile(fileid,dtype='float64',count=self.nth*self.nr).reshape(self.nr,self.nth).transpose()
        
        ## double[nth,nr], lapse
        self.alp= np.fromfile(fileid,dtype='float64',count=self.nth*self.nr).reshape(self.nr,self.nth).transpose()
        
        ## double[nth,nr], entropy
        self.s= np.zeros(dtype='float64',shape=(self.nth,self.nr))
        self.s[:,:mr]= np.fromfile(fileid,dtype='float64',count=self.nth*mr).reshape(mr,self.nth).transpose()
        
        ## double[nth,nr], temperature
        self.T= np.zeros(dtype='float64',shape=(self.nth,self.nr))
        self.T[:,:mr]= np.fromfile(fileid,dtype='float64',count=self.nth*mr).reshape(mr,self.nth).transpose()

    ## nan-ify the matter outside the surface
    def nanify(self):
        """nan-ify the matter outside the surface"""
        for i in range(self.nth):
            j= self.surf[i]
            self.vphi[i,j:]= self.rho[i,j:]= self.p[i,j:]= self.hden[i,j:]= self.s[i,j:]= self.T[i,j:]= np.nan
        return
    
    ## set to zero the matter outside the surface
    def reset(self):
        """set to zero the matter outside the surface"""
        for i in range(self.nth):
            j= self.surf[i]
            self.vphi[i,j:]= self.rho[i,j:]= self.p[i,j:]= self.hden[i,j:]= self.s[i,j:]= self.T[i,j:]= 0.
        return
    
    ## set the matter quantities in the first radial point outside the surface for interpolation
    def setout(self):
        """set the matter quantities in the first radial point outside the surface for interpolation"""
        for i in range(self.nth):
            j= self.surf[i]
            # rho, s = 0 outside the surface
            self.rho[i,j:]= self.p[i,j:]= self.hden[i,j:]= self.s[i,j:]= self.T[i,j:]= 0.
            # vphi continous outside the surface
            self.vphi[i,j:]= self.vphi[i,j-1]
        return
    
    ## differentiate a quantity in the radial direction
    def ddr(self,f):
        """differentiate a quantity in the radial direction"""
        nth, nr= np.shape(f)
        
        # array f must have the same shape of binstar
        if(self.nr != nr or self.nth != nth):
            raise SystemExit("error in binstar.ddr(f): wrong f shape")
        
        # create the arrays where to store the finite difference
        df= np.empty(shape=(nth,nr),dtype='float64')
        dr= np.empty(shape=(nr + 1),dtype='float64')
    
        # compute radial grid increments
        # dr[i]= r[i] - r[i-1]
        # with r[-1]= -r[0]
        dr[1:nr]= self.r[1:nr] - self.r[0:nr-1]
        dr[0]= 2.*self.r[0]
        dr[nr]= dr[nr-1]
    
        for i in range(nr):
            # compute the coefficients
            a1= -dr[i+1]/dr[i]/(dr[i] + dr[i+1])
            a2= (dr[i+1] - dr[i])/dr[i]/dr[i+1]
            a3= dr[i]/dr[i+1]/(dr[i] + dr[i+1])
            
            # first point of the stencil
            if(i == 0):
                if(self.ieqth == self.nth - 1):
                    # planar and axial symmetry at the center
                    df[:,i]= a1*f[:,i]
                else:
                    # axial symmetry at the center
                    df[:,i]= a1*f[::-1,i]
            else:
                df[:,i]= a1*f[:,i-1]
            
            #internal point of the stencil
            df[:,i]+= a2*f[:,i]
            
            # outer point of the stencil
            if(i < nr-1):
                df[:,i]+= a3*f[:,i+1]
            else:
                # smoothness at the outer boundary
                df[:,i]+= a3*(f[:,i] + dr[i+1]/dr[i]*(f[:,i] - f[:,i-1]))
    
        return df
    
    ## differentiate a quantity in the angular direction
    def ddth(self,f):
        """differentiate a quantity in the angular direction"""
        # shape of f
        nth, nr= np.shape(f)
        
        # array f must have the same shape of binstar
        if(self.nr != nr or self.nth != nth):
            raise SystemExit("error in binstar.ddth(f): wrong f shape")
        
        # create the arrays where to store the finite difference
        df=  np.empty(shape=(nth,nr),dtype='float64')
        dth= np.empty(shape=nth + 1, dtype='float64')
        
        # compute angular grid increments
        # dth[i]= th[i] - th[i-1]
        # with th[-1]= -th[0]
        # and th[nth]= 2*pi - th[nth-1]
        dth[1:nth]= self.th[1:nth] - self.th[0:nth-1]
        dth[0]= 2.*self.th[0]
        if(self.ieqth == self.nth - 1):
            dth[nth]= dth[nth - 1]
        else:
            dth[nth]= 2.*(np.pi - self.th[nth-1])
    
        for i in range(nth):
            # compute the coefficients
            b1= -dth[i+1]/dth[i]/(dth[i] + dth[i+1])
            b2= (dth[i+1] - dth[i])/dth[i]/dth[i+1]
            b3= dth[i]/dth[i+1]/(dth[i] + dth[i+1])
            # first point of the stencil
            if(i == 0):
                # even at the north pole
                df[i,:]= b1*f[i,:]
            else:
                df[i,:]= b1*f[i-1,:]
            #internal point of the stencil
            df[i,:]+= b2*f[i,:]
            # third point of the stencil
            if(i == nth-1):
                if(self.ieqth == self.nth - 1):
                    # planar symmetry
                    df[i,:]+= b3*f[i-1,:]
                else:
                    # even at the south pole
                    df[i,:]+= b3*f[i,:]
            else:
                df[i,:]+= b3*f[i+1,:]
    
        return df


# ############################################################################# #
# ################## other useful quantities and functions #################### #
# ################## units: G = c = kB = Msol = 1 ############################# #
# ############################################################################# #

## saturation density
rhon= 4.339e-4
## golden ratio
golden_ratio= 0.5*(np.sqrt(5.) + 1.)
## conversion from [Msol] to [ms]
msol2ms= 1.4766275/2.99792458e2
## conversion from [Msol] to [km]
msol2km= 1.4766275
## mass of the neutron
mn= 8.42335e-58
## conversion from [Msol] to [MeV]
msol2mev= 939.565/mn
## conversion from angular velocity in code units to frequency in Hz
omg2khz= 1./(2.*np.pi*msol2ms)

## internal energy density
def ei(a):
    """internal energy density"""
    return a.hden - a.p - a.rho

## cylindrical radius
def R(a):
    """cylindrical radius"""
    tmp_th, tmp_r= np.meshgrid(a.th, a.r, indexing='ij')
    return tmp_r*np.sin(tmp_th)*a.psi**2

## Lorentz factor
def W(a):
    """Lorentz factor"""
    return 1./np.sqrt(1. - (R(a)*a.vphi)**2)

## angular velocity
def Omg(a):
    """angular speed"""
    return a.alp*a.vphi + a.omg

## F function
def F(a):
    """F function"""
    tmp= R(a)**2*a.vphi
    return tmp/(a.alp - tmp*a.alp*a.vphi)

## angular momentum per unit energy
def l(a):
    """angular momentum per unit energy"""
    tmp= R(a)**2*a.vphi
    return tmp/(a.alp + tmp*a.omg)

## potential of the Euler equation
def Q(a):
    """potential of Euler equation"""
    return -0.5*np.log(a.alp**2*(1. - (R(a)*a.vphi)**2))

## residual of the Euler equation (i.e. force balance equation) along the radius
def res_r(a):
    """residual of the Euler equation along r"""
    return a.ddr(Q(a)) - a.ddr(a.p)/a.hden - F(a)*a.ddr(Omg(a))

## residual of the Euler equation (i.e. force balance equation) along the angle
def res_th(a):
    """residual of the Euler equation along th"""
    return a.ddth(Q(a)) - a.ddth(a.p)/a.hden - F(a)*a.ddth(Omg(a))

## angular momentum per unit mass
def j(a):
    """angular momentum per unit mass"""
    return a.hden / a.rho * a.vphi * W(a) * R(a)**2 # Gourgoulhon[2011] Eq:3.85

## compute the stellar global quantities
def global_quantities(a):
    """
    compute and returns the stellar global quantities (c=G=M☉=kB=1):
    * Mg: gravitational mass
    * M0: baryon mass
    * Mp: proper mass
    * J:  angular momentum
    * Tr: rotational energy
    * S:  total entropy
    
    BEWARE:
    the values in the log files are more precise because XNS uses
    gaussian quadrature in the angular direction while this routine
    uses finite differences.
    The difference with 50 angular points and planar symmetry is <0.01%
    for Mg, M0, Mp, and S and much smaller for J and Tr
    """
    
    # if or not the output assumes planar symmetry
    is_planar= (a.ieqth == a.nth-1)
    
    # useful quantities
    WL= W(a)   # Lorentz factor
    Rc= R(a)   # cylindrical radius
    Og= Omg(a) # angular velocity
    
    # initialize the quantities
    Mg= 0. # gravitational mass
    M0= 0. # baryon mass
    Mp= 0. # proper mass
    J= 0.  # angular momentum
    Tr= 0. # rotational energy (== kinetic energy for stationarity)
    S= 0.  # total entropy
    
    # for each angular point
    for i in range(a.nth):
        # compute dth
        if i == 0:
            dth= 0.5*(a.th[i] + a.th[i+1])
        elif i == a.nth-1:
            if is_planar:
                dth= 0.5*(0.5*np.pi - a.th[i-1])
                # th[i] == pi/2
            else:
                dth= np.pi - 0.5*(a.th[i] + a.th[i-1])
        else:
            dth= 0.5*(a.th[i+1] - a.th[i-1])
        if is_planar:
            dth= 2.*dth
        
        # for each radial point
        for j in range(a.surf[i]):
            # compute dr
            if j == 0:
                dr= 0.5*(a.r[j] + a.r[j+1])
            else:
                dr= 0.5*(a.r[j+1] - a.r[j-1])
            # the surface should not reach the end of the radial grid
            
            # useful quantities
            det= 2.*np.pi*a.psi[i,j]**6*a.r[j]**2*np.sin(a.th[i])*dth*dr # determinant
            tmp= a.hden[i,j]*WL[i,j]**2
            v= a.vphi[i,j]*Rc[i,j]
            
            # integrate quantities
            # Gourgoulhon[2011]= arXiv:1003.5015v2
            # gravitational mass from Gourgoulhon[2011] Eq. (4.18)
            Mg= Mg + (a.alp[i,j]*(2.*a.p[i,j] + tmp*(1 + v*v)) \
                + 2.*a.omg[i,j]*tmp*v*Rc[i,j])*det
            # baryon mass from Gourgoulhon[2011] Eq: (4.6)
            M0= M0 + a.rho[i,j]*WL[i,j]*det
            # proper mass from Friedman, Ipser & Parker [1986, ApJ] Eq. (22)
            Mp= Mp + (a.hden[i,j] - a.p[i,j])*WL[i,j]*det
            # rotational energy from Friedman, Ipser & Parker [1986, ApJ] Eq. (20)
            Tr= Tr + 0.5*tmp*v*Og[i,j]*Rc[i,j]*det
            # angular momentum from Gourgoulhon[2011] Eq. (4.39)
            J= J + tmp*v*Rc[i,j]*det
            S= S + a.rho[i,j]/mn*a.s[i,j]*WL[i,j]*det
    
    return (Mg,M0,Mp,J,Tr,S)
            
## angular momentum and energy per unit mass of the ISCO of a Kerr BH
def kerr_isco(m, j):
    """angular momentum and energy per unit mass of the ISCO of a Kerr BH"""
    a= j / m # angular momentum per unit mass of the BH
    # direct (corotating orbits)
    # BPT72= Bardeen, Press & Teukolsky [1972, ApJ]
    # BPT72:eq:2.21
    z1= 1. + (1. - (a / m)**2)**(1. / 3.) * ((1. + a / m)**(1. / 3.) + (1. - a / m)**(1. / 3.))
    z2= np.sqrt(3. * (a / m)**2 + z1**2)
    r= m * (3. + z2 - np.sqrt((3. - z1) * (3. + z1 + 2. * z2))) # ISCO areal radius
    # BPT72:eq:2.12
    denom= r**0.75 * np.sqrt(r**1.5 - 3. * m * np.sqrt(r) + 2. * a * np.sqrt(m))
    e= (r**1.5 - 2. * m * np.sqrt(r) + a * np.sqrt(m)) / denom # ISCO energy per unit mass
    # BPT72:eq:2.13
    l= np.sqrt(m) * (r * r - 2. * a * np.sqrt(m * r) + a * a) / denom # ISCO angular momentum per unit mass
    return (l, e) # (L/m, E/m)

## generate a gnuplot-friendly text file from a binstar object
def bin2dat(a,filename):
    """
    minimal example generating a (gnuplot-friendly) .dat file
    
    EXAMPLE:
    > bin2dat(binstar('star.out'),'star.dat')
    """
    
    # open the output file
    file= open(filename,'w')
    
    # other important quantities
    aei= ei(a)
    
    # define the format string
    frmt= "%E %E %E %E %E %E %E %E %E %E %E\n"
    # remember that in the format string "%d" = integer, "%E" = double
    
    # print the header
    file.write("# 1:th 2:r 3:rest_mass_density 4:pressure 5:internal_energy_density 6:v^phi 7:conformal_factor 8:lapse 9:ZAMO_omega 10:entropy_per_baryon 11:temperature\n")
    file.write("# units: c = G = Msol = kB = 1\n")
    file.write("# rectangular grid nth x nr = %d x %d with " % (a.nth, a.nr))
    if(a.ieqth == a.nth - 1):
        file.write("planar and axial symmetry\n")
    else:
        file.write("axial symmetry\n")
    
    # actually print the output file
    for i in range(a.nth):
        for j in range(a.nr):
            file.write(frmt % \
                       (a.th[i], a.r[j], a.rho[i,j], a.p[i,j], aei[i,j], \
                        a.vphi[i,j], a.psi[i,j], a.alp[i,j], a.omg[i,j], \
                        a.s[i,j], a.T[i,j]))
        file.write("\n")
            
    # close the file
    file.close()
    
    return
