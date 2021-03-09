"""
basic functions for introduction to aerodynamics for engineering students
"""
#%% IMPORTS
from numpy import pi

#%% Classes
class naca:
    """
    Receive the naca number and calculate the desired coefficients.

    Parameters
    ----------
    Number : String
        Must be just the four NACA numbers.
    Variables
    ---------
    You can redeem any of the following variables.
    
    Number  : String
        The four NACA numbers.
    thetaP  : float
        .
    A0      : function
        Returns the A0 coefficient for an alpha angle.
    A1      : float
        Returns the A1 coefficient.
    An      : function
        Returns an An coefficient for a number n.
    Cmc4    : float
        Returns the pitching moment in 25% of the chord behind the leading edge of the airfoil.
    CL      : function
        Returns the lift coefficient cl for an alpha angle.
    Returns
    -------
    None.
    """
    def __init__(self, Number):
       
        #Imports
        import numpy as np
        
        #Errors 
        if type(Number) == int:
            Number = str(Number)
            
        assert len(Number)  <= 8   ,"The NACA number must be less than or equal to 8 digits!"
        assert len(Number)  >= 4   ,"The NACA number must be greater than or equal to 4 digits!"
        
        #The airfoil variables
        self.digits      = len(Number)
        
        if self.digits == 4 :
            self.Number      = Number
            self._M          = int( Number[0] )
            self._m          = self._M/100
            self._P          = int( Number[1] )
            self._p          = self._P/10
            self._T          = int( Number[2:4] )
            self._t          = self._T/100
            
            self._a          = (0.2969, -0.1260, -0.3516, 0.2843, -0.1015, -0.1036)
        else:
            assert False, "Unrecognized NACA!!"
            
        #Useful functions
        
        self.__linspace  = lambda a,b,c: np.linspace(a,b,c)
        self.__arange    = lambda a,b,c: np.arange(a,b,c)
        self.__zeros     = lambda a: np.zeros(a)
        self.__sin       = lambda a: np.sin(a)
        self.__cos       = lambda a: np.cos(a)
        self.__atan      = lambda a: np.arctan(a)
        self.__sqrt      = lambda a: np.sqrt(a)
        #For airfoil Propeties
        if self.digits == 4:
            self.thetaP     = np.arccos( 1 - 2 * self._p )
            m               = self._m
            p               = self._p
            thetaP          = self.thetaP
            def An(n, alpha):
                if   n ==0:
                    cA0     = ( m*(2*p -1)/pi ) * ( thetaP/(p*p) + ( pi - thetaP )/( (1-p)*(1-p) ) ) + ( m * np.sin( thetaP )/pi ) * ( 1/(p*p) - 1/( (1-p)*(1-p) ) )
                    A0      = lambda alpha: alpha - cA0
                    return A0(alpha)
                elif n ==1:
                    A1      = (1/pi) * ( ( m/(p*p) ) * ( 2*(2*p-1)*np.sin( thetaP ) + thetaP + np.sin(thetaP)*np.cos(thetaP) ) + ( m/(((1-p)*(1-p))) )*( pi - thetaP -2*(2*p-1)*np.sin(thetaP) - np.sin(thetaP)*np.cos(thetaP) ) ) 
                    return A1
                else:
                    An      = lambda n: 0.5*( 2*(2*p-1)*np.sin(n*thetaP)/n + np.sin((n-1)*thetaP)/(n-1) + np.sin((n+1)*thetaP)/(n+1) ) * ( 2*m/(pi*p*p) - 2*m/(pi*(1-p)*(1-p)) )
                    return An(n)
            self.An         = lambda n, alpha=0: An(n,alpha)
            self.Cmc4       = (self.An(2) - self.An(1))*pi/4
            self.cl         = lambda alpha: 2 * pi * (An(0,alpha) + self.An(1)/2)
    
    def profile(self, name = "",chord = True,numberOfPoints=500, export = False):
        """
        plots the airfoil profile (red) and the chord (green).
        Parameters
        ----------
        name    : String, optional
            The title of the graphic. The default is "".
        chord   : Boolean, optional
            To draw the chord line. The default is True.
        numberOfPoints: Integer, optional
            Number of points in the plot. The default is 500.
        export  : Boolean, optional
            To export the graphic. The default is False.
        
        Returns
        -------
        None.

        """
        m       = self._m
        p       = self._p
        t       = self._t
        x       = self.__linspace(0,1,numberOfPoints)
        a       = self._a
        
        if self.digits == 4:
            yc      = self.__zeros(numberOfPoints)
            yt      = self.__zeros(numberOfPoints)
            yu      = self.__zeros(numberOfPoints)
            yl      = self.__zeros(numberOfPoints)
            xu      = self.__zeros(numberOfPoints)
            xl      = self.__zeros(numberOfPoints)
            dyc_dx  = self.__zeros(numberOfPoints)
            theta   = self.__zeros(numberOfPoints)
            
            
            for i in range(numberOfPoints):
                #Camber and Gradient
                if x[i] >= 0 and x[i] <= p:
                    yc[i]       = (m/(p*p)) * (2*p*x[i] - x[i]*x[i])
                    dyc_dx[i]   = 2*(m/(p*p)) * (p - x[i])
                    
                elif x[i] > p and x[i] <= 1:
                    yc[i]       = (m/((1-p)**2)) * (1 - 2*p + 2*p*x[i] - (x[i]*x[i]))
                    dyc_dx[i]   = 2*(m/((1-p)**2)) * (p-x[i])
                
                theta[i]    = self.__atan(dyc_dx[i])
        
                #Thickness distribution
                yt[i]       = 5 * t * ( a[0]*self.__sqrt(x[i]) + a[1]*x[i] + a[2] * (x[i]**2) + a[3] * (x[i]**3) + a[4] * (x[i]**4) )
                
                #Upper surface points
                xu[i]       = x [i] - yt[i] * self.__sin( theta[i] )
                yu[i]       = yc[i] + yt[i] * self.__cos( theta[i] )
                
                #Lower surface points
                xl[i]       = x [i] + yt[i] * self.__sin( theta[i] )
                yl[i]       = yc[i] - yt[i] * self.__cos( theta[i] )
                
        
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=[10,10])
        if chord == True:
            plt.plot(x ,yc,'forestgreen',linewidth=1)
        plt.plot(xu,yu,'r',linewidth=2)
        plt.plot(xl,yl,'r',linewidth=2)
        plt.ylabel("Zc")
        plt.axis('equal')
        plt.xlabel("Wing Length(%)")
        if name == "":
            name =  ("NACA "+self.Number+ " Airfoil Upper Profile")
        plt.title(name)
        
        if (export == True):
            plt.savefig(name)
        plt.show()
        
    def plotCl(self, alphaMax = pi, dalpha = pi/180, unit = "deg", export = False, name = ""):
        """
        plots the approximation of the lift coefficient (CL) as function of the angle of attack.
        Parameters
        ----------
        alphaMax: float, optional
            Maximum angle of attack. The default is pi.
        dalpha  : float, optional
            The range of the angles in radian. The default is pi/180.
        unit    : string, optional
            Cold be rad for radian or deg for degrees. The default is "deg".
        export  : Boolean, optional
            To export the graphic. The default is False.
        name    : String, optional
            The title of the graphic. The default is "".

        Returns
        -------
        None.

        """
        assert type(unit) == str, "String must be passed!"
        
        unit = unit.lower()
        assert unit == "deg" or unit == "rad", "Unit must be radian(rad) or degrees(deg)!"
        
        rangeAl = self.__arange(0,alphaMax+dalpha,dalpha)
        i       = 0
        cl      = []
        
        for alpha in rangeAl:
            cl.append( 2 * pi * (self.An(0,alpha) + self.An(1)/2) )
            i += 1
        
        import matplotlib.pyplot as plt
        
        if name == "":
            name =  ("NACA "+self.Number+ " Airfoil Upper Profile")
        plt.title(name)
            
        if unit == "deg":
            graus = rangeAl * 180 /pi
            plt.plot(graus,cl)
            plt.xlabel("α in deg")
            plt.ylabel("cl")
        
        if unit == "rad":
            plt.plot(rangeAl,cl)
            plt.xlabel("α in rad")
            plt.ylabel("cl")
        
        if (export == True):
            plt.savefig(name)
            
        plt.show()

#%% Application
       
a = naca(2412)

a.profile()

