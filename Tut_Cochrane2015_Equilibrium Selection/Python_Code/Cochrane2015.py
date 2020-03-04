# -*- coding: utf-8 -*-
"""
In this program, I reproduce the IRFs from the Cochrane working paper 2015:
    "The New Keynesian Liquidity Trap"
"""

print ("\n"*80)
import numpy as np
from numpy import linalg as LA # lineare Algebra: Berechnung der Eigenwerte

import matplotlib.pyplot as plt

# =============================================================================
# Declare Variables and Parameters
# =============================================================================

# Date when economy returns out of the ZLB
T           = 5

# define a list
# r           = [0.05, 0]
# i_n         = [0, 0.05]

# define an array
r           = np.array([0.05, 0])
i_n         = np.array([0, 0.05])
rho         = 0.05
kappa       = 1
# where the first element in the vectors display the values for periods of ZLB

'''
Eine Matrix mit Python zu erstellen ist sehr untershiedlich zu der 
Herangehensweise von Matlab. Eine Matrix per se gibt es hier nicht. 
Was dem gleich kommt ist ein multidimensionaler Array. Für diesen Zweck
ist es hilfreich das numpy Paket zu laden. So wird eine Matrix erstellt durch
das Einfügen einzelner arrays.
'''

# system matrix
A = []
A = np.array([[0, -1] , 
              [-kappa, rho]])
print(A) 
np.size(A)
np.size(A,1) # zählt die columns/Spalten
np.size(A,0) # zählt die rows/Zeilen

evalue, evec = LA.eig(A)
# evalue = LA.eigvals(A)
# Man kann auch eigvals nutzen:
#   delta, lambda = La.eigvals(A)
print(evalue)
np.size(evalue)
#np.size(evalue,1) # zählt die columns/Spalten
#np.size(evalue,0) # zählt die rows/Zeilen
#
delta       = evalue[0]

lamb        = evalue[1]



# Compute coefficient matrix linked to the general solution
B = np.array([[lamb**2, -delta**2],
             [lamb,     -delta]])
B = B/(lamb-delta)

'''
In Python kann man bei einem array auch direkt den Datentyp angeben, z.B.
np.array([1,2,3,4],dtype='float32')
'''

# Create time index
t = np.arange(2*T*100)
#t+=1: Add on period to ensure that it counts to 10
t = t/100
t = np.append(t,10)    # Füge einen Eintrag hinzu


'''
Hier führen viele Wege zum Ziel. Am einfachsten ist jedoch der folgende.
Der letzte Eintrag in np.arange() gibt die Schrittweite an
'''
# Alternativ (Buch Seite 57)
t2 = np.arange(0,2*T,0.01) # mit Schrittweite 0.01
# np.linspace() gibt es auch, wie bei Matlab


N               = np.size(t)
N_len           = len(t) # gängiger
# wo befindet sich unser T in t? Gib mir den index
T_loc_t         = np.where(t == T)
T_check         = t[T_loc_t ] # Ergebnis ist ein array
T_check         = T_check[0]  # Ruf den ersten Eintrag ab

# =============================================================================
# Simulations
# =============================================================================

# Initialisiere (wie in Matlab)
S1  =  np.zeros([2,N], dtype = float)
S2  =  np.zeros([2,N], dtype = float)
S3  =  np.zeros([2,N], dtype = float)

for i in range(N):
    if t[i]<T:
        # Deflationary solution
       ir       = i_n[0] - (-r[0])
       Const    = np.array([[rho],[1]])*ir
       C        = np.array([
           [np.exp(delta*(t[i]-T))],
           [np.exp(lamb*(t[i]-T))]
           ])*ir
       X        = Const - np.dot(B,C)
       # np.dot berechnet das Matrix-Produkt
       '''
       np.dot() entspricht matrixmultiplikation
       np.multiply() entpspricht der element-wise Multiplikation
       z.B.:
           x1 = np.arange(9.0).reshape((3, 3))
           x2 = np.arange(3.0)
           np.multiply(x1, x2)
           array([[  0.,   1.,   4.],
                  [  0.,   4.,  10.],
                  [  0.,   7.,  16.]])
           Produkt von (3X3) und (1X3) element-wise
       '''
       
       S1[:,i]    += X[:,0] # returns an array: deswegen die Einträge direkt 
                            # abrufen zum Auffüllen des S1 Arrays
       del X
       # Backward stable case
       X        = Const + delta/(lamb-delta)*np.array([[delta],[1]])*ir*np.exp(lamb*(t[i]-T))
       S2[:,i]    = X[:,0]
       '''
       Result ist ein array. Arrays lassen sich nicht ohne Weiteres in andere 
       Arrays speicherneinfügen, deswegen ruft man die items direkt auf und 
       speichert den Inhalt in einem Array
       '''
    else:
       pi_T     = 0
       S1[:,i]  += pi_T
       # Backward stable case
       X  = ((lamb)/(lamb-delta))*ir*np.array([[lamb], [1]])*np.exp(delta*(t[i]-T))
       S2[:,i]    += X[:,0]


# =============================================================================
# Plot
# =============================================================================
fig = plt.figure()
ax = fig.gca()
ax.plot(t,S1[1,:])
ax.plot(t,S2[1,:])
ax.plot([0,t[-1]],[0,0],'--k', lw=1)
ax.plot([T,T],[-10,10],'--k',  lw=1)
# index [-1] is the last item in array 
ax.set_xlim([-0.0, t[-1]])
ax.set_ylim([-0.2, 0.1])
fig.savefig('Eq_Selection.eps')
plt.show()



'''
Wie füllt man arrays mit nans?
np.full(3, np.nan)
'''















