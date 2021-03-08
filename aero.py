
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
            
        assert type(Number) == str , "Integer must be passed!"
        assert  len(Number) ==   4 , "The string must contain 4 digits!"
        
        #Variables
        self.Number      = Number
        self._m          = int( Number[0] )/100
        self._p          = int( Number[1] )/10
        
        self.__linspace  = lambda a,b,c: np.linspace(a,b,c)
        self.__arange    = lambda a,b,c: np.arange(a,b,c)
            
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
    
    def uProfile(self, export = False, name = ""):
        """
        plots the upper airfoil profile.
        Parameters
        ----------
        export  : Boolean, optional
            To export the graphic. The default is False.
        name    : String, optional
            The title of the graphic. The default is "".

        Returns
        -------
        None.

        """
        m       = self._m
        p       = self._p
        
        x       = self.__linspace(0,p,100)
        y       = self.__linspace(p,1,100)
        z1      = ( m / (p*p) ) * ( 2*p*x-(x*x) )
        z2      = ( m / ((p-1)*(p-1)) ) * ( 2*p*y - y*y + (1-2*p) )
        
        import matplotlib.pyplot as plt
        
        plt.plot(x,z1, 'r')
        plt.plot(y,z2, 'b')
        plt.xlim(0, 1)
        plt.ylim(0, 0.1)
        plt.ylabel("Zc")
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




