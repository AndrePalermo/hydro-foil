
import numpy as np
import matplotlib.pyplot as plt

def Pz_polarization(polarization_file,vorticity=False, spin_to_polarization=2):
	'''
	Computes the z component of polarization as a function of phi from the particlizationCalc output.
	If vorticity = True plots only the contribution of the vorticity
	spin_to_polarization converts the mean spin formula to polarization. Particles are assumet to be S=1/2 by default.
	'''
	filename =  np.loadtxt(polarization_file, unpack=True)
	pT = filename[0]
	phi = filename[1]
	dndp = filename[2]
	dimP=np.size(np.unique(pT))
	dimPhi=np.size(np.unique(phi))
	if(vorticity):
        Pizu = filename[6]
	else:
        Pizu = filename[6] + filename[10]  

	Pz = Pizu.reshape((dimP,dimPhi))
	dNdP = dndp.reshape((dimP,dimPhi))
	Pt = pT.reshape((dimP,dimPhi))
	mean_spin = np.trapz(Pz*Pt,x=np.unique(pT),axis=0)
	spectra = np.trapz(dNdP*Pt,x=np.unique(pT),axis=0)

	pol = spin_to_polarization*mean_spin/spectra

	return np.unique(phi), pol

def Pj_polarization(polarization_file,vorticity=False, mass=1.115683 ,spin_to_polarization=2):
	'''
	Computes the component of polarization along J as a function of phi from the particlizationCalc output.
	If vorticity = True plots only the contribution of the vorticity
	spin_to_polarization converts the mean spin formula to polarization.
	The particle at hand is assumed to be a Lambda by default.
	'''
	filename =  np.loadtxt(polarization_file, unpack=True)
	pT = filename[0]
	phi = filename[1]
	dndp = filename[2]
	px = pT * np.cos(phi)
	py = pT * np.sin(phi)
	dimP=np.size(np.unique(pT))
	dimPhi=np.size(np.unique(phi))
	if(vorticity):
    	Pixu = filename[4]
    	Piyu = filename[5]
	else:
    	Pixu = filename[4] + filename[8]
    	Piyu = filename[5] + filename[9]  

	#boost to the Particle rest frame
	Pix_rf = np.zeros(len(Piyu))
	Piy_rf = np.zeros(len(Piyu))
	for i in range(len(Piyu)):
    	e = np.sqrt(mass*mass + pT[i]*pT[i])
    	Pix_rf[i] = Pixu[i] - (Pixu[i]px[i] + Piyu[i]*py[i])/(e(e+mass))*px[i]  
    	Piy_rf[i] = Piyu[i] - (Pixu[i]px[i] + Piyu[i]*py[i])/(e(e+mass))*py[i]
   	 
	Py = Piy_rf.reshape((dimP,dimPhi))
	dNdP = dndp.reshape((dimP,dimPhi))
	Pt = pT.reshape((dimP,dimPhi))
	mean_spin = np.trapz(Py*Pt,x=np.unique(pT),axis=0)
	spectra = np.trapz(dNdP*Pt,x=np.unique(pT),axis=0)

	pol = spin_to_polarization*mean_spin/spectra

	return np.unique(phi), -pol #The experimental y axis is opposite to the angular momentum


plt.figure("1")
plt.plot(*Pz_polarization(non_proj_ng),label=r"$\varpi+\xi$, non projected")
plt.plot(*Pz_polarization(proj_ng),label=r"$\varpi+\xi, projected$")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$P_z$")
plt.legend()
plt.xlim(0,np.pi)
plt.xticks([0,np.pi/2,np.pi],["0",r"$\pi/2$",r"$\pi$"])
plt.yticks([0])
plt.grid()

plt.figure(2)
plt.plot(*Pz_polarization(non_proj_ng,vorticity=True),label=r"$\varpi$, non projected")
plt.plot(*Pz_polarization(proj_ng,vorticity=True),label=r"$\varpi$, projected")
# plt.plot(*pl.Pz_polarization(foils,vorticity=True),label=r"foils $\varpi$")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$P_z$")
plt.legend()
plt.xlim(0,np.pi)
plt.xticks([0,np.pi/2,np.pi],["0",r"$\pi/2$",r"$\pi$"])
plt.yticks([0])
plt.grid()

plt.figure("3")
plt.plot(*Pj_polarization(non_proj_ng),label=r"$\varpi+\xi$, non projected")
plt.plot(*Pj_polarization(proj_ng),label=r"$\varpi+\xi$, projected")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$P_j$")
plt.legend()

plt.xlim(0,np.pi/2)

# plt.plot(*pl.Pj_polarization(foils),label=r"foils $\varpi+\xi$")
plt.figure("4")
plt.plot(*Pj_polarization(non_proj,vorticity=True),label=r"$\varpi$, non projeced")
plt.plot(*Pj_polarization(proj,vorticity=True),label=r"$\varpi$, projected")
# plt.plot(*pl.Pj_polarization(foils,vorticity=True),label=r"foils $\varpi$")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$P_j$")
plt.legend()

plt.xlim(0,np.pi/2)
