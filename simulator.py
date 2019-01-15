#################################################################
#                      Projet Modelisation                      #
#                        Chute de dominos                       #
#################################################################
"""Description"""

__author__ = "Jean Schneider"
__version__ = "1"
__date__="15/01/2019"


###################
# PARAMETRES DE SIMULATION

nb_points = 10000 # nb de points temporels
nb_dominos = 30 # nb de dominos
duration = 10 # en s
s = 15e-3 # espace entre les dominos en m
d = 5e-3 # epaisseur d'un domino en m
h = 30e-3 # hauteur d'un domino en m
omega_0 = 4 # vitesse initiale en rad/s
m = 15e-3 # masse d'un domino
g = 9.8 # acceleration de pesanteur
display_coef = 3000 # coef d'agrandissement
screen_h = 400 # hauteur de lecran
screen_w = nb_dominos*s*display_coef #largeur de l'ecran

####################
# VARIABLES

# theta[i,t] : position angulaire du ie domino en rad à t
# omega[i,t] : vitesse angulaire du ie domino en rad/s


####################
# IMPORTS

import numpy as np
from tkinter import *
import matplotlib.pyplot as plt

###################
# FONCTIONS

def p(j,theta):
    '''fonction p telle que definie dans le document. Cette fonction est recursive'''

    # Condition d'arret
    if j==0:
        return theta

    res = p(j-1, theta)

    return res + np.arcsin( ((s+d)*np.cos(res)-d)/h )

###################
#Implémentation des frottements 

def a(n,j,i) :
    return(h*np.cos(theta[j,i]-theta[j+1,i])-(s+d)*np.sin(theta[j+1,i])-0.35*d)

def b(n,j,i) :
    return(h*np.cos(theta[j,i]-theta[j+1,i])+0.35*h*np.sin(theta[j,i]-theta[j+1,i]))

def r(n,j,i) :
    if j == 0 :
        return 1
    else :
        return(r(n,j-1,i)*a(n,n-j,i)/b(n,n-j,i))

def J(n,i):
	'''inertie'''
	res=0
	for k in range(n+1) :
		res+=omega[k,i]*r(n,k,i)
	return res


def J_moins_un(n,i) :
	'''terme d'inertie'''
	res=0
	for k in range(1,n):
		res+=omega[k,i]*r(n,k,i)
	return res 
    

###################

def w(j,theta):
    '''fonction w telle que definie dans le document. Cette fonction est recursive'''

    # Condition d'arret
    if j==0:
        return 1

    res = w(j-1, theta)

    return res*(1- ((s+d)*np.sin(p(j, theta)))/(h*np.cos(p(j, theta) - p(j-1, theta))))


###################

def K(n, omega,t):
    ''' energie cinetique des dominos'''
    s=0
    for i in range(n+1):
        s+=omega[i][t]**2

    return I/2*s

###################

def H(n, theta, i):
    '''hauteur combinee des centres de masse des dominos'''
    s=0
    for j in range(n+1):
        s+=(np.cos(theta[j,i]) - np.cos(theta_inf) + (d/h)*(np.sin(theta[j,i]) - np.sin(theta_inf)))

    return s


###################

def Inertia(n, theta, i):
    s=0
    for j in range(n+1):
        s+=w(j,theta[j,i])**2 # a voir avec j et n...
    #print('inertia', s)
    return s

###################

def eq23(n):
    '''coefficient multiplicateur de l equation 23'''
    s1=0
    s2=0
    for j in range(n+1):
        s1+=w(j,theta_c)

    for i in range(n+2):
        s2+=w(j,0)

    return s1/s2
    

###################

def eq34(n,i) :
	'''coefficient multiplicateur de l equation 34'''
	return J_moins_un(n,i)/J(n,i)

def eq35(n, theta, i):
    '''coefficient multiplicateur de l equation 35'''

    return Inertia(n-1,theta,i)/Inertia(n, theta, i)

###################

def update_frame(window, canvas, dominos_tk, t):
	'''gestion graphique'''
	global nb_dominos, end

	for i in range(nb_dominos):
		#creation du rectangle avec des lignes (car il tourné...)

		angle=theta[i][t]
		x1=d*(i+1)+i*s
		y1=0

		x2=x1-d*np.cos(angle)
		y2=d*np.sin(angle)

		x4=x1+h*np.sin(angle)
		y4=h*np.cos(angle)

		x3=x2+h*np.sin(angle)
		y3=y4+y2


		# mise a lechelle
		x1=x1*display_coef
		x2=x2*display_coef
		x3=x3*display_coef
		x4=x4*display_coef

		y1=y1*display_coef
		y2=y2*display_coef
		y3=y3*display_coef
		y4=y4*display_coef


		#rotation car le canvas est dans le mauvais sens

		y1=screen_h - y1
		y2=screen_h - y2
		y3=screen_h - y3
		y4=screen_h - y4

		#update de tk
		canvas.coords(dominos_tk[i][0], x1,y1,x2,y2)
		canvas.coords(dominos_tk[i][1], x1,y1,x4,y4)
		canvas.coords(dominos_tk[i][2], x2,y2,x3,y3)
		canvas.coords(dominos_tk[i][3], x3,y3,x4,y4)

	if t<end:
		window.after(int(dt*5000), update_frame, window, canvas, dominos_tk, t+1)


###################
# MAIN

# Creation du temps
timespace = np.linspace(0,duration, nb_points)
dt = timespace[1]-timespace[0]
I = m/3 *(h**2 + d**2)


theta=np.zeros((nb_dominos,nb_points))
omega=np.zeros((nb_dominos,nb_points))

# Conditions initiales
omega[0][0] = omega_0
n = 0 # foremost falling domino
e=1 + I/(m*g*h)*omega_0**2 # energie initiale
#print("ernegie init", e)

theta_c = np.arcsin(s/h)
theta_inf = np.arccos(d/(s+d))
print("theta_inf", theta_inf)

dominos_tk = []

vitesse_front=[]
impact_time=[]

window = Tk()

canvas = Canvas(window, width=screen_w, height=screen_h)

# creation graphique des dominos

for i in range(nb_dominos):
	result = []
	result.append(canvas.create_line((i*s)*display_coef, 0,(i*s+d)*display_coef, 0))
	result.append(canvas.create_line((i*s+d)*display_coef, 0, (i*s+d)*display_coef, h*display_coef))
	result.append(canvas.create_line((i*s)*display_coef, 0, (i*s)*display_coef, h*display_coef))
	result.append(canvas.create_line((i*s)*display_coef, h*display_coef,(i*s+d)*display_coef, h*display_coef))
	dominos_tk.append(result)

canvas.pack()



# Debut de recursion
# une boucle de temps
for i in range(1, nb_points):
    t = timespace[i]

    theta[n, i] = theta[n, i - 1] + dt * omega[n, i - 1]

    # boucle pour s'occuper de tous les dominos qui chutent
    for j in range(n):
        theta[j, i] = p(n - j, theta[n, i])

    # pour le dernier domino, il est en "chute libre":

    temp=(e - H(n, theta, i))*m*g*h / (I*Inertia(n, theta, i)) # (17)
    omega[n,i] = temp**0.5

    # on controle si le domino n+1 est touché ou pas encore
    if theta[n,i] >= theta_c:
        # alors on fait evoluer le domino suivant egalement

        if n < nb_dominos-1:
            print("==== Le domino {} a poussé le suivant ====".format(n))
            n += 1

            coef=eq35(n,theta,i)
            omega[n,i] = omega[n-1,i] * coef  # (23)
            vitesse_front.append(omega[n,i]*h)
            impact_time.append(timespace[i])

            # energie de l'ensemble
            e = H(n, theta, i) + I/(m*g*h) * Inertia(n, theta, i) * omega[n,i]**2 -min(n, 6)*0.0981

        else:
        	end=i
        	break


print("=== Please wait, loading graphics... ===")


update_frame(window, canvas, dominos_tk, 1)
window.mainloop()


plt.plot(impact_time,vitesse_front)
plt.xlabel("Temps")
plt.ylabel("Vitesse du front en m/s")
plt.show()
