Aerofólios Espessos - Joukowski

Do potencial complexo
===

Temos que $\mu = k$ e $Q_\infty = U_0$

e consequentemente temos os escoamentos:

Uniforme: $F_u(Y) = Q_\infty Y e^{i\alpha}$
Vórtice: $F_v(y) = \frac{i\Gamma}{2\pi} \ln{(Y-Y_0)}$
Dipolo: $F_d(Y) = \frac{\mu}{2\pi (Y-Y_0)}e^{i\alpha}$

Assim podemos representar o escoamento de um cilindro rotativo em um escoamento uniforme como:
$F(Y) = F_u(Y)+F_v(Y)+F_d(Y) = Q_\infty Y e^{i\alpha} + \frac{i\Gamma}{2\pi} \ln{(Y-Y_0)} + \frac{\mu}{2\pi (Y-Y_0)}e^{i\alpha}$

e o campo de velocidade pode ser representado como:

$W(Y) = \frac{dF(Y)}{dY} = \frac{d}{dY}( Q_\infty Y e^{i\alpha} + \frac{i\Gamma}{2\pi} \ln{(Y-Y_0)} + \frac{\mu}{2\pi (Y-Y_0)}e^{i\alpha} )$

$W(Y) = Q_\infty e^{-i\alpha} - \frac{Q_\infty a^2}{(Y-Y_0)^2}e^{i\alpha} + \frac{i\Gamma}{2\pi(Y-Y_0)}$

Transformação conformal
===

- A ideia da transformação conformal é mapear o escoamento ao redor de um aerofólio no domínio físico $(x,y)$ em um escoamento ao redor de um cilindro rotativo no domínio do cilindro $(g,h)$.
- Assim asolução para o escoamento do cilindro pode ser utilizado para obter a solução ao redor do perfil.
- Assim define-se uma transformação de Joukowski: $Y(f)=f+ \frac{C^2}{16f}$


![](https://www.grc.nasa.gov/www/k-12/airplane/Images/map.gif)

- ### Para impor a condição de Kutta no borto de fuga força-se que o ponto correspondente ao bordo de fuga no círculo seja uma estagnação do círculo.
- ### Ajusta-se o parâmetro $\Gamma$ para que $W(f_{te}) =0$
- ### Sendo $f_{te} = C/4$ a circulação do vórtice resulta: $\Gamma = 4\pi a Q_\infty \sin{(\alpha + \beta)}$
	- $\beta = arg(\frac{f_{te} - \mu}{a})$
- ### Assim a velocidade complexa no domínio do círculo será:
	- $W(f) = Q_\infty e^{-i\alpha} - \frac{Q_\infty a^2}{(f-\mu)^2} e^{i\alpha} + \frac{i2aQ_\infty \sin{(\alpha + \beta)}}{(f-\mu)}$

