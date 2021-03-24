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
            self.alphaL0    = - (self.An(0) + self.An(1)/2)
    
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

class joukowski:
    
    '''
        
        Parameters
        ----------
        m : float, optional
            flow center. The default is 0.1.
        epsilon : float, optional
            Airfoil thickness. The default is 0.1.
        c : float, optional
            Chord. The default is 1.0.
        U0 : float, optional
            One-dimensional flow. The default is 1.0.
        alpha : float, optional
            Angle in degrees. The default is 3.
        nPts : integer, optional
            Number of points. The default is 50.
        
        Variables
        ---------
        You can redeem any of the following variables.
        
        m : float
            flow center.
        epsilon : float
            Airfoil thickness.
        c : float
            Chord.
        U0 : float
            One-dimensional flow.
        alpha : float
            Angle in degrees.
        nPts : integer
            Number of points.
        j : complex
            complex unit.
        mu : complex
            Center of circumference.
        a : float
            Circumference radius.
        beta : float
            .
        theta : numpy array
            Distribution of circumference angles.
        f : Complex numpy array
            Circumference geometry.
        Y : Complex numpy array
            Conformal transformation of cylinder points.
        Gamma : float
            Circulation.
        k : float
            Dipole intensity.
        Wc : complex numpy array
            Cylinder speed.
        Cpc : numpy array
            Cylinder CP.
        Wa : complex numpy array
            airfoil speed.
        Cpa : numpy array
            Airfoil CP.
        Cla : float
            Airfoil Cl.
        Clnum : float
            Numeric airfoil Cl.
        ClNumericError : float
            Numerical error of the Cl value for the airfoil.
        chordNumeric : float
            Numeric chord value.
        chordError : float
            Numerical error of the chord value for the airfoil.
            
        Returns
        -------
        None.

        '''
    
    def __init__(self,m = 0.1, epsilon = 0.1,c = 1.0, U0 = 1.0, alpha = 3, nPts = 50):
        
        import numpy as np
        import matplotlib.pyplot as plt
        
        self.__np       = np
        self.__plt      = plt
        
        self.j = 1j
        self.m          = m                                                     # Centro do escoamento
        self.epsilon   = epsilon                                                # Espessura do aerofólio
        self.c          = c                                                     # Corda        
        self.U0         = U0                                                    # Velocidade do escoamento
        self.alpha      = alpha * pi/180                                        # Ângulo de ataque em Deg
        self._nPts      = nPts                                                  # Número de pontos
        
        self.mu         = c/4 * (- epsilon + m*self.j)                          # Centro da circunferência 
        self.a          = abs( c/4 - self.mu )                                  # Raio da circunferência      
        self.beta       =-np.angle(c/4 - self.mu)                               # 
        self.theta      = np.linspace(0,2*pi,nPts)                              # Distribuição dos ângulos da circunferência
        self.f          = self.mu + self.a * np.exp(self.j * self.theta)        # Geometria da circunferência
        
        self.Y          = self.f + (self.c**2)/(16*self.f)                      # Transformação conformal do cilindro
        
        self.Gamma      = 4 * pi*self.a * self.U0 * np.sin(self.alpha+self.beta)# Circulação
        self.k          = 2 * pi*(self.a**2)*U0                                 # Intensidade do dipolo
        
        #Calculo para o cilindro
        
        self.Wc         = 2 * self.U0*1j*np.exp(1j*self.theta)*( np.sin(self.alpha+self.beta)-np.sin(self.alpha-self.theta) ) #Velocidade do cilindro
        self.Cpc        = 1 - ( abs(self.Wc)**2 )/( self.U0**2 )                # Calculo do CP para o cilindro
        
        #Calculo para o Aerofólio
        
        self.Wa         = self.Wc/( 1 - (self.c**2)/(16*self.f**2) )            # Velocidade do aerofólio
        self.Cpa        = 1 - ( abs(self.Wa)**2 )/( self.U0**2 )                # Cp do aerofólio
        
        self.Cla        = 2*pi* np.sqrt((1+self.epsilon)**2+self.m**2)*np.sin(self.alpha+self.beta)#Cl do aerofólio
        #Integração numérica e erros
        self.Clnum      = np.trapz(self.Cpa,self.Y.real)/self.c                 # Cálculo do cl numérico
        self.ClNumericError = (self.Clnum - self.Cla)/self.Cla
        
        self.chordNumeric   = (self.Y.real.max() - self.Y.real.min())           # Corda calculada
        self.chordError = ( self.chordNumeric - self.c)/self.c                  #Erro da corda
        
    def plotGeometry(self,color = '#00f', g = 'all', title = "Geometry", export = False):
        '''
        
        
        Parameters
        ----------
        g : str, optional
            Type of plot. The default is 'all'.
                all -> circumference and airfoil geometry
                a   -> airfoil geometry
                c   -> circumference geometry
        color : hex string, optional
            The hex color string. The default is '#00f'.
        title : str, optional
            The title of the graphic. The default is "Geometry".
        export : boolean, optional
            To export the graphic. The default is False.

        Returns
        -------
        None.

        '''
        
        plt = self.__plt
        
        g = g.upper()
        
        if g == 'ALL' :
            fig = plt.figure()
            plt.subplot(1,2,1)
            plt.plot(self.f.real,self.f.imag,color)
            plt.axis('equal')
            plt.grid(True)
        
            ax =plt.subplot(1,2,2)
            ax.yaxis.tick_right()
            plt.plot(self.Y.real,self.Y.imag,color)
            plt.axis('equal')
            plt.grid(True)
            
            plt.suptitle(title)
            fig.add_subplot(111, frame_on=False)
            plt.tick_params(labelcolor="none", bottom=False, left=False)
            plt.xlabel("$\dfrac{x}{C}$")
            plt.ylabel("$\dfrac{y}{C}$")
            
            if (export == True):
                plt.savefig(title)
            plt.show()
        
        elif g == 'C' :
            
            plt.plot(self.f.real,self.f.imag,color)
            plt.axis('equal')
            plt.grid(True)
            
            plt.title(title)
            if (export == True):
                plt.savefig(title)
            plt.show()
                        
        elif g == 'A' :
            
            plt.plot(self.Y.real,self.Y.imag,color)
            plt.axis('equal')
            plt.grid(True)
            
            plt.title(title)
            if (export == True):
                plt.savefig(title)
            plt.show()
        
        else:
            assert False, "Unrecognized g parameter"
        
    def plotCp(self, g = 'all',color='#00f', title = "Cp value", export = False):
        '''
        
        Parameters
        ----------
        g : str, optional
            Type of plot. The default is 'c'.
                all -> circumference and airfoil geometry
                a   -> airfoil geometry
                c   -> circumference geometry
        color : hex string, optional
            The hex color string. The default is '#00f'.
        title : str, optional
            The title of the graphic. The default is "Geometry".
        export : boolean, optional
            To export the graphic. The default is False.

        Returns
        -------
        None.

        '''
        plt = self.__plt 
        g = g.upper()
        
        if g == 'ALL' :
            fig = plt.figure()
            plt.subplot(1,2,1)
            plt.plot(self.f.real,self.Cpc, color)
            plt.axis('equal')
            plt.grid(True)

            ax = plt.subplot(1,2,2)
            plt.plot(self.Y.real,self.Cpa, color)
            plt.axis('equal')
            ax.yaxis.tick_right()
            plt.grid(True)
            
            plt.suptitle(title)
            fig.add_subplot(111, frame_on=False)
            plt.tick_params(labelcolor="none", bottom=False, left=False)
            plt.xlabel("$\dfrac{x}{C}$")
            plt.ylabel("$C_p$ Number")
            
            if (export == True):
                plt.savefig(title)
            plt.show()
            
        elif g == 'C' :
            plt.plot(self.f.real,self.Cpc, color)
            plt.title(title)
            plt.xlabel('x/c')
            plt.ylabel('Cp')
            plt.axis('equal')
            
            if (export == True):
                plt.savefig(title)
            plt.show()
            
        elif g == 'A' :
            plt.plot(self.Y.real,self.Cpa, color)
            plt.title(title)
            plt.xlabel('x/c')
            plt.ylabel('Cp')
            plt.axis('equal')
            
            if (export == True):
                plt.savefig(title)
            plt.show()   
    
    def plotStreamLines(self, g = 'all', nLines = 35, JustStream = False, AirfoilColor = '#f00', streamColor=['#666','#aaa'], title = "Stream lines", export = False):
        '''
        

        Parameters
        ----------
        g : str, optional
            type of plot. The default is 'all'.
                all -> circumference and airfoil geometry
                a   -> airfoil geometry
                c   -> circumference geometry
                
                *If the word mesh is added to the end of the sentence, the generated graphic will be only the mesh.*
        
        nLines : integer, optional
            Number of lines. The default is 35.
        JustStream : boolean, optional
            If true, only the flow will be plotted (the profile will not be plotted). The default is False.
        AirfoilColor : hex color string, optional
            The hex color string. The default is '#f00'.
        streamColor : List of the hex color string , optional
            List of the hex color string. The default is ['#666','#aaa'].
        title : str, optional
            The title of the graphic. The default is "Geometry".
        export : boolean, optional
            To export the graphic. The default is False.

        Returns
        -------
        None.

        '''
        
        mesh = False
        if len(g) > 3:
            g = g.split()[0]
            mesh = True 
            
        g               = g.upper()
        
        np              = self.__np
        plt             = self.__plt
        
        theta           = self.theta
        r               = np.linspace(self.a,4*self.a,self._nPts)               # Cria pontos com raio 4*a
        
        T,R             = np.meshgrid(theta,r)                                  # Cria uma malha de pontos intermediários entre a e 4*a
        f               = self.mu + R * np.exp(1j * T)                          # Geometria das circunferências concentricas
        Y               = f + (self.c**2)/(16*f)                                # Aplica a transformação conformal para formar aerofólios concêntricos
        
        Fu              = self.U0 * f * np.exp(-1j*self.alpha)                  # Equação do potencial complexo escoamento uniforme
        Fv              = self.Gamma/(2*pi) *1j* np.log(f-self.mu)              # Equação do potencial complexo vórtice
        Fd              = self.k/(2*pi) * np.exp(1j*self.alpha)/(f-self.mu)     # Equação do potencial complexo dipolo
        
        F               = Fu + Fv + Fd                                          # Equação do potencial complexo total
        
        if mesh == True :
            pulo = int(self._nPts/nLines)
            if g == 'ALL' :
                fig = plt.figure()
                plt.subplot(1,2,1)
                plt.plot(f.real[1::pulo],f.imag[1::pulo],AirfoilColor)
                plt.plot(f.real[1::pulo].T,f.imag[1::pulo].T,AirfoilColor)
                plt.axis('equal')
            
                ax = plt.subplot(1,2,2)
                ax.yaxis.tick_right()
                plt.plot(Y.real[1::pulo],Y.imag[1::pulo],AirfoilColor)
                plt.plot(Y.real[1::pulo].T,Y.imag[1::pulo].T,AirfoilColor)
                plt.axis('equal')
                
                plt.suptitle(title)
                fig.add_subplot(111, frame_on=False)
                plt.tick_params(labelcolor="none", bottom=False, left=False)
                plt.xlabel("$\dfrac{x}{C}$")
                plt.ylabel("$\dfrac{y}{C}$")
                    
                if (export == True):
                    plt.savefig(title)
                plt.show()
            
            elif g == 'C' :
                plt.plot(f.real[1::pulo],f.imag[1::pulo],AirfoilColor)
                plt.plot(f.real[1::pulo].T,f.imag[1::pulo].T,AirfoilColor)
                plt.title(title)
                
                if (export == True):
                    plt.savefig(title)
                plt.axis('equal')   
                plt.show()
                            
            elif g == 'A' :
                plt.plot(Y.real[1::pulo],Y.imag[1::pulo],AirfoilColor)
                plt.plot(Y.real[1::pulo].T,Y.imag[1::pulo].T,AirfoilColor)
                plt.title(title)
                
                if (export == True):
                    plt.savefig(title)
                plt.axis('equal')  
                plt.show()
            
            else:
                assert False, "Unrecognized g-mesh parameter"
                
        elif g == 'ALL' :
            fig = plt.figure()
            plt.subplot(1,2,1)
            plt.contour(f.real,f.imag,F.imag,nLines,colors=streamColor)
            plt.axis('equal')
            if JustStream == False:
                plt.plot(self.f.real,self.f.imag,AirfoilColor)
        
            ax = plt.subplot(1,2,2)
            ax.yaxis.tick_right()
            plt.contour(Y.real,Y.imag,F.imag,nLines,colors=streamColor)
            plt.axis('equal')
            if JustStream == False:
                plt.plot(self.Y.real,self.Y.imag,AirfoilColor)
             
            plt.suptitle(title)
            fig.add_subplot(111, frame_on=False)
            plt.tick_params(labelcolor="none", bottom=False, left=False)
            plt.xlabel("$\dfrac{x}{C}$")
            plt.ylabel("$\dfrac{y}{C}$")
            
            if (export == True):
                plt.savefig(title)
            plt.show()
        
        elif g == 'C' :
            plt.contour(f.real, f.imag, F.imag, nLines,colors=streamColor)
            plt.title(title)
            if JustStream == False:
                plt.plot(self.f.real,self.f.imag,AirfoilColor)
            
            if (export == True):
                plt.savefig(title)
            plt.axis('equal')   
            plt.show()
                        
        elif g == 'A' :
            plt.contour(Y.real, Y.imag, F.imag, nLines,colors=streamColor)
            plt.title(title)
            if JustStream == False:
                plt.plot(self.Y.real,self.Y.imag,AirfoilColor)
            
            if (export == True):
                plt.savefig(title)
            plt.axis('equal')  
            plt.show()

        else:
            assert False, "Unrecognized g parameter"
                  
#%% Application
'''      
a = naca(2412)

a.profile()
'''

A = joukowski(nPts = 100)

A.plotGeometry()
A.plotGeometry(g='a',title='Airfoil geometry')

A.plotCp()
A.plotCp(g='a',title='Airlfoil Cp')

A.plotStreamLines(g='all', nLines=30)
A.plotStreamLines(g='a',nLines=30 ,title = 'Airfoil stream lines')

A.plotStreamLines(g='all mesh',nLines = 5, title = 'meshes')
A.plotStreamLines(g='a mesh',nLines = 5, title = 'Airfoil Mesh')

print( str(A.chordError*100)[0:4] + '% de erro na corda A' )
print( str(A.ClNumericError*100)[0:5] + '% de erro numérico A')
