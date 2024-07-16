import matplotlib.pyplot as plt
# import plot_lib as pl
import numpy as np
import matplotlib
import matplotlib.patches as ptch

def Sz_rapidity_polarization(polarization_file, mass=1.115683, rap_cut=1, pt_cut=10, vorticity=False, global_integrals = False):
    '''
    Computes the z component of polarization as a function of phi from the particlizationCalc output.
    If vorticity = True plots only the contribution of the vorticity
    spin_to_polarization converts the mean spin formula to polarization. Particles are assumet to be S=1/2 by default.
    '''
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<=rap_cut) & (pT>0.5) & (pT<pt_cut)
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]

    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Pi0 = filename[4]
    Pizu = filename[7]
    if(vorticity):
        vorticity=False
    else:
        Pi0 = Pi0 + filename[8]
        Pizu = Pizu + filename[11] 

    Pi0 = Pi0[select]
    Pizu = Pizu[select]

    #BACKBOOST to Lambda RF
    mT = np.sqrt(mass*mass + pT*pT)
    Pizu = Pizu - Pi0*(mT*np.sinh(y_rap))/(mT*np.cosh(y_rap)+mass)
    
    if dimy>1:
        Pz = Pizu.reshape((dimP,dimPhi,dimy))
        dNdP = dndp.reshape((dimP,dimPhi,dimy))
        Pt = pT.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Pz,x=np.unique(pT),axis=0),x=np.unique(y_rap),axis=1)
        spectra = np.trapz(np.trapz(dNdP,x=np.unique(pT),axis=0),x=np.unique(y_rap),axis=1)
    else:
        print("Midrapidity only!")
        Pz = Pizu.reshape((dimP,dimPhi))
        dNdP = dndp.reshape((dimP,dimPhi))
        Pt = pT.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Pz,x=np.unique(pT),axis=0)
        spectra = np.trapz(dNdP,x=np.unique(pT),axis=0)
    
    if(global_integrals):
        return np.unique(phi), mean_spin, spectra

    pol = mean_spin/spectra

    return np.unique(phi), pol

def Sz_polarization_pT(polarization_file, mass=1.115683, harmonics=2, rap_cut=1, vorticity=False, global_integrals = False):
    '''
    Computes the z component of polarization as a function of pT from the hydro-foil output.
    If vorticity = True plots only the contribution of the vorticity
    '''
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<=rap_cut) 
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]

    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Pi0 = filename[4]
    Pizu = filename[7]
    if(vorticity):
        vorticity=False
    else:
        Pi0 = Pi0 + filename[8]
        Pizu = Pizu + filename[11] 

    Pi0 = Pi0[select]
    Pizu = Pizu[select]

    #BACKBOOST to Lambda RF
    mT = np.sqrt(mass*mass + pT*pT)
    Pizu = Pizu - Pi0*(mT*np.sinh(y_rap))/(mT*np.cosh(y_rap)+mass)
    
    if dimy>1:
        Pzsin = Pizu*np.sin(harmonics*phi) 
        Pz_reahsped = Pzsin.reshape((dimP,dimPhi,dimy))
        dNdP_reshaped = dndp.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Pz_reahsped,x=np.unique(phi),axis=1),x=np.unique(y_rap),axis=1)
        spectra = np.trapz(np.trapz(dNdP_reshaped,x=np.unique(phi),axis=1),x=np.unique(y_rap),axis=1)
    else:
        print("Midrapidity only!")
        Pzsin = Pizu*np.sin(harmonics*phi) 
        Pz_reahsped = Pzsin.reshape((dimP,dimPhi))
        dNdP_reshaped = dndp.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Pz_reahsped,x=np.unique(phi),axis=1)
        spectra = np.trapz(dNdP_reshaped,x=np.unique(phi),axis=1)
    
    if(global_integrals):
        return np.unique(pT), mean_spin, spectra

    pol = mean_spin/spectra

    return np.unique(pT), pol

def Sj_rapidity_polarization(polarization_file, mass=1.115683, rap_cut=1, pt_cut=10, vorticity=False, global_integrals = False):
    '''
    Computes the component of polarization along J as a function of phi from the particlizationCalc output.
    If vorticity = True plots only the contribution of the vorticity
    spin_to_polarization converts the mean spin formula to polarization. 
    The particle at hand is assumed to be a Lambda by default.
    '''
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<rap_cut) & (pT>0.5) & (pT<pt_cut)
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]
    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Pi0 = filename[4]
    Piyu = filename[6]
    if(vorticity):
        vorticity=False
    else:
        Pi0 = Pi0 + filename[8]
        Piyu = Piyu + filename[10]  

    Pi0 = Pi0[select]
    Piyu = Piyu[select]

    #BACKBOOST to Lambda RF
    mT = np.sqrt(mass*mass + pT*pT)
    Piyu = Piyu - Pi0*(pT*np.sin(phi))/(mT*np.cosh(y_rap)+mass)

    if dimy>1:
        Py = Piyu.reshape((dimP,dimPhi,dimy))
        dNdP = dndp.reshape((dimP,dimPhi,dimy))
        Pt = pT.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Py,x=np.unique(pT),axis=0),x=np.unique(y_rap),axis=1)
        spectra = np.trapz(np.trapz(dNdP,x=np.unique(pT),axis=0),x=np.unique(y_rap),axis=1)
    else:
        print("Midrapidity only!")
        Py = Piyu.reshape((dimP,dimPhi))
        dNdP = dndp.reshape((dimP,dimPhi))
        Pt = pT.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Py,x=np.unique(pT),axis=0)
        spectra = np.trapz(dNdP,x=np.unique(pT),axis=0)
    
    if(global_integrals):
        return np.unique(phi), -mean_spin, spectra

    pol = mean_spin/spectra

    return np.unique(phi), -pol

def Sj_polarization_pT(polarization_file, mass=1.115683, harmonics=2, rap_cut=1, vorticity=False, global_integrals = False):
    '''
    Computes the z component of polarization as a function of pT from the hydro-foil output.
    If vorticity = True plots only the contribution of the vorticity
    '''
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<=rap_cut) 
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]

    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Pi0 = filename[4]
    Piyu = filename[6]
    if(vorticity):
        vorticity=False
    else:
        Pi0 = Pi0 + filename[8]
        Piyu = Piyu + filename[10]  

    Pi0 = Pi0[select]
    Piyu = Piyu[select]

    #BACKBOOST to Lambda RF
    mT = np.sqrt(mass*mass + pT*pT)
    Piyu = Piyu - Pi0*(pT*np.sin(phi))/(mT*np.cosh(y_rap)+mass)

    if dimy>1:
        Py = Piyu.reshape((dimP,dimPhi,dimy))
        dNdP = dndp.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Py,x=np.unique(phi),axis=1),x=np.unique(y_rap),axis=1)
        spectra = np.trapz(np.trapz(dNdP,x=np.unique(phi),axis=1),x=np.unique(y_rap),axis=1)
    else:
        print("Midrapidity only!")
        Py = Piyu.reshape((dimP,dimPhi))
        dNdP = dndp.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Py,x=np.unique(phi),axis=1)
        spectra = np.trapz(dNdP,x=np.unique(phi),axis=1)
    
    if(global_integrals):
        return np.unique(pT), -mean_spin, spectra

    pol = mean_spin/spectra

    return np.unique(pT), -pol

def Feed_down_Sz_rapidity_polarization(polarization_file, rap_cut=1, pt_cut=10, vorticity=False, global_integrals = False):
    '''
    Computes the z component of polarization as a function of phi from the particlizationCalc output.
    If vorticity = True plots only the contribution of the vorticity
    spin_to_polarization converts the mean spin formula to polarization. Particles are assumet to be S=1/2 by default.
    '''
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<rap_cut) & (pT>0.5) & (pT<pt_cut)
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]

    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Piz = filename[6]
    if(vorticity):
        vorticity=False
    else:
        Piz = Piz + filename[9] 

    Pizu = Piz[select] ##already in the rest frame

    if dimy>1:
        Pz = Pizu.reshape((dimP,dimPhi,dimy))
        dNdP = dndp.reshape((dimP,dimPhi,dimy))
        Pt = pT.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Pz,x=np.unique(pT),axis=0),x=np.unique(y_rap),axis=1)
        spectra = np.trapz(np.trapz(dNdP,x=np.unique(pT),axis=0),x=np.unique(y_rap),axis=1)
    else:
        print("Midrapidity only")
        Pz = Pizu.reshape((dimP,dimPhi))
        dNdP = dndp.reshape((dimP,dimPhi))
        Pt = pT.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Pz,x=np.unique(pT),axis=0)
        spectra = np.trapz(dNdP,x=np.unique(pT),axis=0)

    if(global_integrals):
        return np.unique(phi), mean_spin, spectra
    
    pol = mean_spin/spectra

    return np.unique(phi), pol

def Feed_down_Sz_pT(polarization_file, rap_cut=1, vorticity=False, global_integrals = False, harmonics=2):
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<rap_cut)
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]

    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Piz = filename[6]
    if(vorticity):
        vorticity=False
    else:
        Piz = Piz + filename[9] 

    Pizu = Piz[select] ##already in the rest frame

    if dimy>1:
        Pzsin =Pizu*np.sin(harmonics*phi)
        Pz_reshaped = Pzsin.reshape((dimP,dimPhi,dimy))
        dNdP_reshaped = dndp.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Pz_reshaped,x=np.unique(phi),axis=1),x=np.unique(y_rap),axis=1)
        spectra = np.trapz(np.trapz(dNdP_reshaped,x=np.unique(phi),axis=1),x=np.unique(y_rap),axis=1)
    else:
        print("Midrapidity only")
        Pzsin =Pizu*np.sin(harmonics*phi)
        Pz_reshaped = Pzsin.reshape((dimP,dimPhi))
        dNdP_reshaped = dndp.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Pz_reshaped,x=np.unique(phi),axis=1)
        spectra = np.trapz(dNdP_reshaped,x=np.unique(phi),axis=1)

    if(global_integrals):
        return np.unique(pT), mean_spin, spectra
    
    pol = mean_spin/spectra

    return np.unique(pT), pol

def Feed_down_Sj_pT(polarization_file, rap_cut=1, vorticity=False, global_integrals = False):
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<rap_cut)
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]

    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Piy = filename[5]
    if(vorticity):
        vorticity=False
    else:
        Piy = Piy + filename[8] 

    Piyu = Piy[select] ##already in the rest fraem

    if dimy>1:
        Py = Piyu.reshape((dimP,dimPhi,dimy))
        dNdP = dndp.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Py,x=np.unique(y_rap),axis=2),x=np.unique(phi),axis=1)
        spectra = np.trapz(np.trapz(dNdP,x=np.unique(y_rap),axis=2),x=np.unique(phi),axis=1)
    else:
        print("Midrapidity only!")
        Py = Piyu.reshape((dimP,dimPhi))
        dNdP = dndp.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Py,x=np.unique(phi),axis=1)
        spectra = np.trapz(dNdP,x=np.unique(phi),axis=1)
        
    if(global_integrals):
        return np.unique(pT), -mean_spin, spectra

    pol = mean_spin/spectra

    return np.unique(pT), -pol


def Feed_down_Sj_rapidity_polarization(polarization_file, rap_cut=1, pt_cut=10, vorticity=False, global_integrals = False):
    '''
    Computes the component of polarization along J as a function of phi from the particlizationCalc output.
    If vorticity = True plots only the contribution of the vorticity
    spin_to_polarization converts the mean spin formula to polarization. 
    The particle at hand is assumed to be a Lambda by default.
    '''
    filename =  np.loadtxt(polarization_file, unpack=True)
    pT = filename[0]
    phi = filename[1]
    y_rap = filename[2]
    select = (np.abs(y_rap)<rap_cut) & (pT>0.5) & (pT<pt_cut)
    pT = pT[select]
    phi = phi[select]
    y_rap = y_rap[select]

    dndp = filename[3]
    dndp = dndp[select]
    dimP=np.size(np.unique(pT))
    dimPhi=np.size(np.unique(phi))
    dimy=np.size(np.unique(y_rap))
    Piy = filename[5]
    if(vorticity):
        vorticity=False
    else:
        Piy = Piy + filename[8]
        
    Piyu = Piy[select] ##already in the rest fraem

    if dimy>1:
        Py = Piyu.reshape((dimP,dimPhi,dimy))
        dNdP = dndp.reshape((dimP,dimPhi,dimy))
        Pt = pT.reshape((dimP,dimPhi,dimy))
        mean_spin = np.trapz(np.trapz(Py,x=np.unique(y_rap),axis=2),x=np.unique(pT),axis=0)
        spectra = np.trapz(np.trapz(dNdP,x=np.unique(y_rap),axis=2),x=np.unique(pT),axis=0)
    else:
        print("Midrapidity only!")
        Py = Piyu.reshape((dimP,dimPhi))
        dNdP = dndp.reshape((dimP,dimPhi))
        Pt = pT.reshape((dimP,dimPhi))
        mean_spin = np.trapz(Py,x=np.unique(pT),axis=0)
        spectra = np.trapz(dNdP,x=np.unique(pT),axis=0)
        
    if(global_integrals):
        return np.unique(phi), -mean_spin, spectra

    pol = mean_spin/spectra

    return np.unique(phi), -pol

def global_Lambda_spin(file,**kwargs):
    phi_z, Pz_prim, spectrz = Sz_rapidity_polarization(file, global_integrals=True, **kwargs)
    phi_j, Pj_prim, spectrj = Sj_rapidity_polarization(file, global_integrals=True, **kwargs)

    Pzsin2phi = np.trapz(Pz_prim*np.sin(2*phi_z),x=phi_z,axis=0)
    Pj = np.trapz(Pj_prim,x=phi_j,axis=0)
    dnz = np.trapz(spectrz,x=phi_j,axis=0)
    dnj = np.trapz(spectrj,x=phi_j,axis=0)

    return Pj/dnj,Pzsin2phi/dnz 

def global_feed_down_spin(file,**kwargs):
    phi_z, Pz_, spectrz = Feed_down_Sz_rapidity_polarization(file, global_integrals=True, **kwargs)
    phi_j, Pj_, spectrj = Feed_down_Sj_rapidity_polarization(file, global_integrals=True, **kwargs)

    Pzsin2phi = np.trapz(Pz_*np.sin(2*phi_z),x=phi_z,axis=0)
    Pj = np.trapz(Pj_,x=phi_j,axis=0)
    dnz = np.trapz(spectrz,x=phi_j,axis=0)
    dnj = np.trapz(spectrj,x=phi_j,axis=0)

    return Pj/dnj,Pzsin2phi/dnz

def centrality_polarization(dir_list, feed_down = True,**kwargs):
    centrality_dict = {"0-5" : 2.5,
                       "5-10" : 7.5,
                       "10-20" : 15,
                       "20-30" : 25,
                       "30-40" : 35,
                       "40-50" : 45,
                       "50-60" : 55,
                       "60-70" : 65,
                       "70-80" : 75}
    cent_ = np.array([])
    cent_Pz = np.array([])
    cent_Pj = np.array([])

    fraction_primary = 0.243
    fraction_sigma0 = 0.275*0.6
    fraction_sigmastar = 0.359
    tot = fraction_primary+fraction_sigma0+fraction_sigmastar

    for dir in dir_list:

        Pj, Pzsin2phi = global_Lambda_spin(f"{dir}/primary",**kwargs)

        for range in centrality_dict.keys():
            if range in dir:
                centrality = centrality_dict[range]
                cent_ = np.append(cent_,centrality)

        if(feed_down):
            Pjsigma0, Pzsin2phisigma0 = global_feed_down_spin(f"{dir}/FeedDown_Sigma0",**kwargs)
            Pjstar, Pzsin2phistar = global_feed_down_spin(f"{dir}/FeedDown_Sigmastar",**kwargs)
            cent_Pz = np.append(cent_Pz,(fraction_primary*Pzsin2phi+fraction_sigma0*Pzsin2phisigma0+fraction_sigmastar*Pzsin2phistar)/tot)
            cent_Pj = np.append(cent_Pj,(fraction_primary*Pj+fraction_sigma0*Pjsigma0+fraction_sigmastar*Pjstar)/tot)
        else:
            cent_Pz = np.append(cent_Pz,Pzsin2phi)
            cent_Pj = np.append(cent_Pj,Pj)
            
    sort_cent = np.argsort(cent_)    
    cent_ = np.take_along_axis(cent_,sort_cent,axis=0)
    cent_Pz = np.take_along_axis(cent_Pz,sort_cent,axis=0)
    cent_Pj = np.take_along_axis(cent_Pj,sort_cent,axis=0)

    return (cent_,cent_Pj), (cent_,cent_Pz)

def experimental_plot(x, y, yerr, yerr_syst=0, xerr=0, dx=0.05, linestyle="none", marker="o", **kwargs):

    figure = plt.errorbar(x, y, yerr, xerr=xerr, linestyle=linestyle, marker=marker, **kwargs)

    col = figure[0].get_color()
    lw = figure[0].get_linewidth()
    if np.size(yerr_syst) != 1:
        for i in range(np.size(y)):
            dy = yerr_syst[i]
            if np.size(xerr) > 1:
                dx = xerr[i]

            center = (x[i], y[i])
            anchor = (center[0] - dx, center[1] - dy)

            plt.gca().add_patch(ptch.Rectangle(anchor, 2 * dx, 2 * dy,
                                        edgecolor=col,
                                        fill=False,
                                        lw=lw))
    return plt.gca()
