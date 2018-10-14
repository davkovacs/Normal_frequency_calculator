import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d.axes3d import Axes3D


class Atoms(object):
    '''Object containing the geometry and the energy of each species'''
    def __init__(self, r, theta, E):
        self.r=r
        self.theta=theta
        self.E=E
    def get_r(self):
        return self.r
    def get_r_SI(self):
        return self.r*A_to_m
    def get_theta(self):
        return self.theta
    def get_theta_SI(self):
        return self.theta*degree_to_rad
    def get_E(self):
        return self.E
    def get_E_SI(self):
        return self.E*au_to_J
    def __str__(self):
        return "<Atoms, r="+str(self.r)+", theta="+str(self.theta)+", E="+str(self.E)+">"
    def __lt__(self,other):
        return self.E<other.E


def plot_potential_energy_surface(atoms):
    '''Creates a 3D surface plot of the potential energy surface if given a list of Atoms objects'''
    XX=[xx.get_r() for xx in atoms]
    YY=[xx.get_theta() for xx in atoms]
    ZZ=[xx.get_E() for xx in atoms]

    ax = plt.axes(projection='3d')
    ax.set_xlabel('r / A')
    ax.set_ylabel('Theta / degrees')
    ax.set_zlabel('Energy / au')
    ax.set_title('Potential Energy Surface of '+molecule)
    ax.view_init(35,290)
    ax.plot_trisurf(XX, YY, ZZ,cmap=mpl.cm.coolwarm)
    plt.savefig(molecule+'_PES.pdf')
    plt.show()
    plt.close()
    return None

def quadratic_fit_plotter(xval,yval, polynom, constant_direction,xlabel):
    '''Creates a plot to show the quality of the quadatic fit'''
    xfit=np.arange(min(xval)-5e-12,max(xval)+5e-12,(max(xval)-min(xval))/1e4)
    yfit=np.polyval(polynom,xfit)
    plt.scatter(xval,yval,color="red", marker="x")
    plt.plot(xfit, yfit)
    plt.title("PES along constant "+constant_direction)
    plt.xlabel(xlabel)
    plt.ylabel("Energy / J")
    axes = plt.gca()
    axes.set_xlim([min(xval)-(max(xval)-min(xval))*0.05,max(xval)+(max(xval)-min(xval))*0.05])
    axes.set_ylim([min(yval)-(max(yval)-min(yval))*0.05,max(yval)+(max(yval)-min(yval))*0.05])
    fig = plt.gcf()
    fig.set_size_inches(9, 4)
    plt.savefig(molecule+"_const_"+constant_direction+".pdf")
    plt.show()
    plt.close()
    return None

def get_potential_energy(R,Theta,atomlist):
    '''Given the geometry looks up the atomslist to find the corresponding energy (if contained)'''
    for atom in atomlist:
        if np.isclose(R,atom.get_r(),0.001):
            if np.isclose(Theta,atom.get_theta(),0.001):
                return atom.get_E_SI()
    return print("The specified geometry is not in the list.")

def get_evals(M):
    """Function returning the  eigenvalues of the Hermitian matrix M in ascending order"""
    evals, evecs = np.linalg.eigh(M)
    return evals

#constants:
au_to_J=4.3597e-18
m_u=1.6605e-27
A_to_m=1e-10
degree_to_rad=2*np.pi/360

#input for the molecule type and the location of files
molecule=str(input("Are you interested in H2O or H2S? "))
assert molecule=="H2O" or molecule=="H2S", 'Please type either H2O or H2S'

try:
    folder=str(input("Give me the name of the directory of your " +molecule+ " data! (eg. /dir1/dir2/ ) "))
except:
    print("Invalid user input.")
assert folder[-1]=="/" and folder[0]=="/","Please correct the name of the folder"
    
if molecule=="H2O":
    r_range=[70,191]
elif molecule=="H2S":
    r_range=[60,181]

#reading in the data into a list of Atoms objects
all_geoms=[]
for i in range(r_range[0],r_range[1],5):
    for j in range(70,161):
        f=open(folder+molecule+".r%.2f" %(i/100)+"theta%.1f" %j+".out","r")
        for line in f:
            if "SCF Done:" in line:
                l=line.split()
                all_geoms.append(Atoms(i/100,j,float(l[4])))
        f.close()

plot_potential_energy_surface(all_geoms)

#Calculation of the normal mode vibrational frequencies using quadratic fit in the minima of the PES
eqgeoms=sorted(all_geoms)  # sorts the Atomslist by the energy

print("Quadratic fitting method")

eq_theta_r_SI=[]           #constant eq_theta, r data in m
eq_theta_E_SI=[]           #constant eq_theta E data in J
eq_r_theta_SI=[]           #constant eq_r, theta data in radians
eq_r_E_SI=[]               #constant eq_r, E data in J
ind1=0                     #counts allowing for fitting only in the vicinity of the minima (on the smallest 5 geoms)
ind2=0
for atom in eqgeoms:
        if atom.get_theta()==eqgeoms[0].get_theta() and ind1<5:     #can use == as no numerical calculations done yet
            eq_theta_r_SI.append(atom.get_r()*A_to_m)
            eq_theta_E_SI.append(atom.get_E()*au_to_J)
            ind1+=1
        if atom.get_r()==eqgeoms[0].get_r() and ind2<5:
            eq_r_theta_SI.append(atom.get_theta()*degree_to_rad)
            eq_r_E_SI.append(atom.get_E()*au_to_J)
            ind2+=1

p_const_theta=np.polyfit(eq_theta_r_SI,eq_theta_E_SI,2)   #fitting the quadratic with least squared method
p_const_r=np.polyfit(eq_r_theta_SI,eq_r_E_SI,2)

quadratic_fit_plotter(eq_theta_r_SI,eq_theta_E_SI,p_const_theta,'theta','r / m')  #create plot of the fits
quadratic_fit_plotter(eq_r_theta_SI,eq_r_E_SI,p_const_r,'r','theta / radians')

freq_1=1/2/np.pi/2.998e10*(p_const_theta[0]*2/2/m_u)**0.5                              #calculate the frequencies
freq_2=1/2/np.pi/2.998e10*(p_const_r[0]*2/0.5/m_u/(eqgeoms[0].get_r()*A_to_m)**2)**0.5

print("The vibrational frequency of the symmetric stretching mode is approximately: "+str(round(freq_1,1))+" cm-1")
print("The vibrational frequency of the bending mode is approximately: "+str(round(freq_2,1))+" cm-1")

# using the Hessian to find vibrational frequencies (assuming third normal mode is orthogonal to this 2)
r0=eqgeoms[0].get_r()          #store eq geom values for better readability later
theta0=eqgeoms[0].get_theta()
H=np.empty((2,2))              #Hessian matrix elements given below explicitly by finite difference numerical differentiation at the minimum
H[0][1]=(get_potential_energy(r0+0.05,theta0+1,eqgeoms)-get_potential_energy(r0+0.05,theta0,eqgeoms)-get_potential_energy(r0,theta0+1,eqgeoms)+get_potential_energy(r0,theta0,eqgeoms))/(0.05*A_to_m*degree_to_rad*r0*A_to_m)
H[1][0]=(get_potential_energy(r0-0.05,theta0-1,eqgeoms)-get_potential_energy(r0-0.05,theta0,eqgeoms)-get_potential_energy(r0,theta0-1,eqgeoms)+get_potential_energy(r0,theta0,eqgeoms))/(0.05*A_to_m*degree_to_rad*r0*A_to_m)
H[0][0]=0.5*((get_potential_energy(r0+2*0.05,theta0,eqgeoms)-2*get_potential_energy(r0+0.05,theta0,eqgeoms)+get_potential_energy(r0,theta0,eqgeoms))/(0.05*A_to_m)**2+(get_potential_energy(r0-2*0.05,theta0,eqgeoms)-2*get_potential_energy(r0-0.05,theta0,eqgeoms)+get_potential_energy(r0,theta0,eqgeoms))/(0.05*A_to_m)**2)
H[1][1]=0.5*((get_potential_energy(r0,theta0+2*1,eqgeoms)-2*get_potential_energy(r0,theta0+1,eqgeoms)+get_potential_energy(r0,theta0,eqgeoms))/(1*degree_to_rad)**2+(get_potential_energy(r0,theta0-2*1,eqgeoms)-2*get_potential_energy(r0,theta0-1,eqgeoms)+get_potential_energy(r0,theta0,eqgeoms))/(1*degree_to_rad*r0*A_to_m)**2)
H=0.5*(H+np.transpose(H))      #symmetrise the Hessian to get better values
eigvals=get_evals(H)

freq_1_H=1/2/np.pi/2.998e10*(eigvals[0]/0.5/m_u)**0.5   #calculate the normal mode frequencies
freq_2_H=1/2/np.pi/2.998e10*(eigvals[1]/2/m_u)**0.5

print("Hessian method")
print("The vibrational frequency of the first normal mode is approximately: "+str(round(freq_2_H,1))+" cm-1")
print("The vibrational frequency of the second normal mode is approximately: "+str(round(freq_1_H,1))+" cm-1")





