#include <cmath>
#include <fstream>
#include <array>


#ifdef OPEN_MP
    #include<omp.h>
    const int NTHREADS = omp_get_max_threads();
#endif

#include "particle_data_group.h"
#include "utils.h"
#include "integrals.h"


//for a description of the functions refer to integrals.h
void polarization_midrapidity(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
    double P_vorticity[4] = {0,0,0,0};
    double P_shear[4] = {0,0,0,0};
    int t_vect[4] = {1,0,0,0};
    double Denominator = 0;

    // get particle's info
    const double mass = particle.get_mass();  //GeV
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();

    double mT = sqrt(mass * mass + pT*pT);
    array<double,4> p = {mT, pT *cos(phi), pT *sin(phi), 0};
    array<double,4> p_ = {mT, -pT *cos(phi), -pT *sin(phi), 0}; //lower indices    

    #ifdef OPEN_MP
        int threads_ = NTHREADS; 
        #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,P_vorticity,P_shear)
    #endif
    for(element cell : freeze_out_sup){ //loop over the FO hypersurface
        double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)
        array<double,4> theta_vector = {0,0,0,0};
        double theta_sq = 0.;

        for(array<int,5> indices_an_levi : non_zero_levi()){
            int mu = indices_an_levi[0];
            int nu = indices_an_levi[1];
            int rh = indices_an_levi[2];
            int sg = indices_an_levi[3]; 
            int LeviCivita = indices_an_levi[4];  //levi(mu,nu,rh,sg)
            
            theta_vector[mu] += LeviCivita
                                * p_[sg] * cell.dbeta[nu][rh]/(2*mass)*hbarC; //hbarC makes theta adimensional
		           
        }

        for (int mu = 0; mu < 4; mu++) {
            pdSigma += p[mu] * cell.dsigma[mu];
            pu += p[mu] * cell.u[mu] * gmumu[mu];
            theta_sq += theta_vector[mu]*theta_vector[mu]*gmumu[mu];
        }
        const double mutot = cell.mub*baryonNumber + cell.muq*electricCharge + cell.mus*strangeness;
        const double nf = 1 / (exp( (pu - mutot) / cell.T) + 1.0);

        Denominator += pdSigma * nf ;
        for(int mu=0; mu<4; mu++){
            // computing the vorticity induced polarization
            P_vorticity[mu] += 0.5*pdSigma * nf *  
                        (theta_vector[mu]/sqrt(-theta_sq)) * sinh(sqrt(-theta_sq)*0.5)/
                        (cosh(sqrt(-theta_sq)*0.5)+exp(-(pu - mutot)/cell.T));
        }

        for(array<int,5> indices_an_levi : non_zero_levi()){
            int mu = indices_an_levi[0];
            int nu = indices_an_levi[1];
            int rh = indices_an_levi[2];
            int sg = indices_an_levi[3]; 
            int LeviCivita = indices_an_levi[4];  //levi(mu,nu,rh,sg)

            if(nu==0)
            for(int ta=0; ta<4; ta++){
            P_shear[mu] += - pdSigma * nf * (1. - nf) * LeviCivita*t_vect[nu]* p_[sg] * p[ta] / p[0] 
                        * ( cell.dbeta[rh][ta] + cell.dbeta[ta][rh])/ (8.0 * mass);
            }
        }

    } //end surface loop

    //print to file
    fileout << "   " << pT << "   " << phi << "   " << Denominator;
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_vorticity[mu];
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_shear[mu] *hbarC; //Unit conversion to make the shear adimensional 
    fileout << endl;

}

void polarization_midrapidity_linear(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
    double P_vorticity[4] = {0,0,0,0};
    double P_shear[4] = {0,0,0,0};
    double Denominator = 0;

    // get particle's info
    const double mass = particle.get_mass();  //GeV
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();

    double mT = sqrt(mass * mass + pT*pT);
    array<double,4> p = {mT, pT *cos(phi), pT *sin(phi), 0};
    array<double,4> p_ = {mT, -pT *cos(phi), -pT *sin(phi), 0}; //lower indices    

    #ifdef OPEN_MP
        int threads_ = NTHREADS; 
        #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,P_vorticity,P_shear)
    #endif
    for(element cell : freeze_out_sup){ //loop over the FO hypersurface
        double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)

        for (int mu = 0; mu < 4; mu++) {
            pdSigma += p[mu] * cell.dsigma[mu];
            pu += p[mu] * cell.u[mu] * gmumu[mu];
        }
        const double mutot = cell.mub*baryonNumber + cell.muq*electricCharge + cell.mus*strangeness;
        const double nf = 1 / (exp( (pu - mutot) / cell.T) + 1.0);

        Denominator += pdSigma * nf ;
        for(array<int,5> indices_an_levi : non_zero_levi()){
            int mu = indices_an_levi[0];
            int nu = indices_an_levi[1];
            int rh = indices_an_levi[2];
            int sg = indices_an_levi[3]; 
            int LeviCivita = indices_an_levi[4];  //levi(mu,nu,rh,sg)      
            P_vorticity[mu] += pdSigma * nf * (1. - nf) * LeviCivita
                                    * p_[sg] * cell.dbeta[nu][rh];
                    
            if(nu==0)
            for(int ta=0; ta<4; ta++)
            P_shear[mu] += -pdSigma * nf * (1. - nf) * LeviCivita
                        * p_[sg] * (p[ta] / p[0]) *t_vector[nu]
                        * ( cell.dbeta[rh][ta] + cell.dbeta[ta][rh]);
        }
    }
    //print to file
    fileout << "   " << pT << "   " << phi << "   " << Denominator;
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_vorticity[mu]/ (8.0 * mass) *hbarC; //unit conversion to make the vorticity adimensional
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_shear[mu]/ (8.0 * mass) *hbarC; //Unit conversion to make the shear adimensional 
    fileout << endl;

}

void polarization_exact_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
    double P_vorticity[4] = {0,0,0,0};
    double P_shear[4] = {0,0,0,0};
    double Denominator = 0;

    // get particle's info
    const double mass = particle.get_mass();  //GeV
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();
    const double spin = particle.get_spin();
    const int fermi_or_bose = statistics(spin);
    const double phase_space = ((2*spin +1)/pow( 2*hbarC*PI, 3.0));

    double mT = sqrt(mass * mass + pT*pT);
    array<double,4> p = {mT *cosh(y_rap), pT *cos(phi), pT *sin(phi), mT *sinh(y_rap)};
    array<double,4> p_ = {mT *cosh(y_rap), -pT *cos(phi), -pT *sin(phi), -mT *sinh(y_rap)}; //lower indices    

    #ifdef OPEN_MP
        int threads_ = NTHREADS; 
        #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,P_vorticity,P_shear)
    #endif
    for(element cell : freeze_out_sup){ //loop over the FO hypersurface
        double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)
        array<double,4> theta_vector = {0,0,0,0};
        double theta_sq = 0.;

        for(array<int,5> indices_an_levi : non_zero_levi()){
            int mu = indices_an_levi[0];
            int nu = indices_an_levi[1];
            int rh = indices_an_levi[2];
            int sg = indices_an_levi[3]; 
            int LeviCivita = indices_an_levi[4];//levi(mu,nu,rh,sg)
		                theta_vector[mu] += LeviCivita
                                * p_[sg] * cell.dbeta[nu][rh]/(2*mass)*hbarC; //theta is adimensional
		   }

        for (int mu = 0; mu < 4; mu++){
            pdSigma += p[mu] * cell.dsigma[mu];
            pu += p[mu] * cell.u[mu] * gmumu[mu];
            theta_sq += theta_vector[mu]*theta_vector[mu]*gmumu[mu];
        }
        const double mutot = cell.mub*baryonNumber + cell.muq*electricCharge + cell.mus*strangeness;
        const double distribution = 1 / (exp( (pu - mutot) / cell.T) + fermi_or_bose);

        Denominator += phase_space*pdSigma * distribution ;
        
        for(int mu=0; mu<4; mu++){
            // computing the vorticity induced polarization
            P_vorticity[mu] += phase_space*pdSigma * distribution *  
                        (theta_vector[mu]/sqrt(-theta_sq)) * 
                        aux_exact_polarization(spin, pu, cell.T, mutot, sqrt(-theta_sq));
        }
            // computing the shear induced polarization
            for(array<int,5> indices_an_levi : non_zero_levi()){
            int mu = indices_an_levi[0];
            int nu = indices_an_levi[1];
            int sg = indices_an_levi[2];
            int ta = indices_an_levi[3];
            int LeviCivita = indices_an_levi[4];//levi(mu,nu,sg,ta)

            if(nu==0)
            for(int rh=0; rh<4; rh++){
                P_shear[mu] += - phase_space*(spin/3)*(spin+1)*pdSigma * distribution * (1. - fermi_or_bose*distribution) 
                            * LeviCivita* (p_[ta]/mass) * (p[rh]/p[0]) 
                            *hbarC*( cell.dbeta[sg][rh] + cell.dbeta[rh][sg])/ 2.0;
                            //hbarC: unit conversion to make the shear adimensional 
                    }
        }
    }

    //print to file
    fileout << "   " << pT << "   " << phi << "   " << y_rap << "   " << Denominator;
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_vorticity[mu];
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_shear[mu]; 
    fileout << endl;

}


double aux_exact_polarization(double spin, double pu, double T, double mutot, double abs_theta){
double num = 0;
double den = 1e-20;
int fermi_or_bose = statistics(spin);
for(double k=-spin; k<=spin;k++){
            num += k/(exp((pu-mutot)/T-k*abs_theta)+fermi_or_bose);
            den += 1/(exp((pu-mutot)/T-k*abs_theta)+fermi_or_bose);
        }

if(num/den != num/den){
    cout<<"NaN in aux_exact_polarization!"<<endl;
    exit(1);
}

return num/den;
}

int statistics(double spin){
    //returns 1 if Fermi, -1 if Bose statistics
    const double dim_spin = 2*spin+1;
    if(dim_spin != (int) dim_spin){
        exit(1);
    }
    if((int)dim_spin % 2 == 0){
        return 1;
    }
    else if((int)dim_spin % 2 == 1){
        return -1;
    }
    return 0;
}

void spectrum_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
    double dNd3p = 0;
    
    // get particle's info
    const double mass = particle.get_mass();  
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();
    
    const int spin = particle.get_spin(); //Dimension of the spin Hilbert space: if even -> fermions, if odd -> bosons
    int fermi_or_bose = statistics(spin); //the factor to add in the denominator of the distribution: 1 Fermi-Dirac, -1 Bose-Einstein
    
    double mT = sqrt(mass * mass + pT*pT);
    array<double,4> p = {mT *cosh(y_rap), pT *cos(phi), pT *sin(phi), mT *sinh(y_rap)};
    array<double,4> p_ = {mT *cosh(y_rap), -pT *cos(phi), -pT *sin(phi), -mT *sinh(y_rap)}; //lower indices    

    #ifdef OPEN_MP
        int threads_ = NTHREADS; 
        #pragma omp parallel for num_threads(threads_) reduction(+:dNd3p)
    #endif
    for(element cell : freeze_out_sup){ //loop over the FO hypersurface
        double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)

        for (int mu = 0; mu < 4; mu++) {
            pdSigma += p[mu] * cell.dsigma[mu];
            pu += p[mu] * cell.u[mu] * gmumu[mu];
        }
        const double mutot = cell.mub*baryonNumber + cell.muq*electricCharge + cell.mus*strangeness;
        const double distribution = (1/(pow(2*PI,3))) *1/ (exp( (pu - mutot) / cell.T) + fermi_or_bose);

        dNd3p += pdSigma * distribution ;
        
    }

    //print to file
    fileout << "   " << pT << "   " << phi << "   " << y_rap << "   " << dNd3p/(hbarC*hbarC*hbarC) << endl; //hbarC for unit conversion

}

void polarization_components(pdg_particle particle, std::array<vector<element>,5> components, std::array<string,5> fileout_list){
    int size_pt = 20;
    int size_phi = 30;
    int size_y = 20;
    std::vector<double> pT = linspace(0,6.2,size_pt);
    std::vector<double> phi =  linspace(0,2*PI,size_phi);
    std::vector<double> y_rap =  linspace(-1,1,size_y);

    ofstream fout0(fileout_list[0]),fout1(fileout_list[1]),fout2(fileout_list[2]),
	fout3(fileout_list[3]),fout4(fileout_list[4]);

	for(double ipt : pT){
		for(double iphi : phi){
			for(double iy : y_rap){
                polarization_exact_rapidity(ipt, iphi, iy, particle, components[0], fout0);
                polarization_exact_rapidity(ipt, iphi, iy, particle, components[1], fout1);
                polarization_exact_rapidity(ipt, iphi, iy, particle, components[2], fout2);
                polarization_exact_rapidity(ipt, iphi, iy, particle, components[3], fout3);
                polarization_exact_rapidity(ipt, iphi, iy, particle, components[4], fout4);
            }
        }
    }

}
 

void Lambda_polarization_FeedDown(double pT, double phi, double y_rap, pdg_particle mother, 
    interpolator &spectrum_interpolator, array<interpolator,4> &S_vorticity_interpolator, array<interpolator,4> &S_shear_interpolator, ofstream &fileout){
    pdg_particle lambda(3122);

    double P_vorticity[3] = {0,0,0};
    double P_shear[3] = {0,0,0};
    double Denominator = 1e-20;
    
	const double lambda_mass = lambda.get_mass(); 
	
	//fetch information concerning mother and second son
	const double mother_mass = mother.get_mass();  
	const double mother_baryonCharge = mother.get_b();
	const double mother_electricCharge = mother.get_q();
	const double mother_strangeness = mother.get_s();
    
    pdg_particle* second_son;

    if(mother.get_id()==3212){
        second_son = new pdg_particle(22);
    }
    else if(mother.get_id() == 3224){
        second_son = new pdg_particle(211);
    }
    else{
        cout<< "ERROR! The decay is not implemented yet!"<<endl;
        exit(1);
    }

    double second_son_mass = second_son->get_mass();

    double mT = sqrt(lambda_mass * lambda_mass + pT*pT);
    //momentum of the Lambda particle in the Lab frame
    array<double,4> p = {mT *cosh(y_rap), pT *cos(phi), pT *sin(phi), mT *sinh(y_rap)};
    array<double,4> p_ = {mT *cosh(y_rap), -pT *cos(phi), -pT *sin(phi), -mT *sinh(y_rap)}; //lower indices  

    //Energy and momentum of the Lambda particle in the mother's rest frame
    double E_Lambda_rest = (mother_mass*mother_mass + lambda_mass*lambda_mass - second_son_mass*second_son_mass)/(2*mother_mass); 
    double p_rest_abs = sqrt(E_Lambda_rest*E_Lambda_rest - lambda_mass*lambda_mass);

    int number_of_bins = 30;
    double dangle = PI/number_of_bins;
    #ifdef OPEN_MP
        int threads_ = NTHREADS;
        #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,P_vorticity,P_shear)
    #endif
    for(int itheta=0; itheta<number_of_bins; itheta++){
        double dtheta = itheta*dangle;
        double integralphi_den = 1e-20;
        double integralphi_num[3] = {0,0,0}; 
        double integralphi_numShear[3] = {0,0,0}; 
        for(double dphi=0; dphi<2*PI-1e-5; dphi+=dangle){
            double rest_frame_Pi[3] = {0,0,0}; //polarization in the rest frame of the mother
            double rest_frame_PiShear[3] = {0,0,0}; //polarization in the rest frame of the mother, shear contribution
            double S_vect[3]={0,0,0};//the particular vector depending on the decay. For details 1905.03123
            double S_vectShear[3]={0,0,0};//the particular vector depending on the decay, shear contribution
            double spectrum = 1e-20;
            double polarization_mother[4]={0,0,0,0};
            double polarization_mother_shear[4]={0,0,0,0};
    
            double p_rest[3] = {p_rest_abs*sin(dtheta)*cos(dphi),p_rest_abs*sin(dtheta)*sin(dphi),p_rest_abs*cos(dtheta)}; //spatial momentum of \Lambda in mother rest frame. p^0 is irrelevant to the calculation
            //jacobian from eq. 30 of 1905.03123
            double jacobian = pow(mother_mass,3)*pow(p[0] + E_Lambda_rest,2)*
                                (pow(p[0] + E_Lambda_rest,2)-(lambda_mass*lambda_mass+p[0]*E_Lambda_rest+p[1]*p_rest[0]+p[2]*p_rest[1]+p[3]*p_rest[2]))/
                                (E_Lambda_rest*pow(lambda_mass*lambda_mass+p[0]*E_Lambda_rest+p[1]*p_rest[0]+p[2]*p_rest[1]+p[3]*p_rest[2],3));
            //momentum of the mother
            double Energy_mother;
            double P_mother[4];
            double P_mother_[4];
            //upper indices
            P_mother[1] = (p[1]-p_rest[0])*2*mother_mass*(p[0] + E_Lambda_rest)/
                        (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[0],2)-pow(p[2]-p_rest[1],2)-pow(p[3]-p_rest[2],2));
            P_mother[2] = (p[2]-p_rest[1])*2*mother_mass*(p[0] + E_Lambda_rest)/
                        (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[0],2)-pow(p[2]-p_rest[1],2)-pow(p[3]-p_rest[2],2));
            P_mother[3] = (p[3]-p_rest[2])*2*mother_mass*(p[0] + E_Lambda_rest)/
                        (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[0],2)-pow(p[2]-p_rest[1],2)-pow(p[3]-p_rest[2],2));
            Energy_mother = sqrt(mother_mass*mother_mass+pow(P_mother[1],2)+pow(P_mother[2],2)+pow(P_mother[3],2));
            P_mother[0] = Energy_mother; 
            //lower indices
            P_mother_[1] = -P_mother[1];
            P_mother_[2] = -P_mother[2];
            P_mother_[3] = -P_mother[3];
            P_mother_[0] = Energy_mother;
            
            double pt_mom = sqrt(P_mother[1]*P_mother[1]+P_mother[2]*P_mother[2]);
            double phi_mom = atan2(P_mother[2],P_mother[1]); //azimuthal angle of the mother [0,2\pi]
            if(phi_mom<0){
                phi_mom += 2*PI;
            }
            double y_mom = atanh(P_mother[3]/P_mother[0]);

            spectrum = spectrum_interpolator.trilinear_interpol(pt_mom, phi_mom, y_mom);

            polarization_mother[0] = S_vorticity_interpolator[0].trilinear_interpol(pt_mom, phi_mom, y_mom);
            polarization_mother[1] = S_vorticity_interpolator[1].trilinear_interpol(pt_mom, phi_mom, y_mom);
            polarization_mother[2] = S_vorticity_interpolator[2].trilinear_interpol(pt_mom, phi_mom, y_mom);
            polarization_mother[3] = S_vorticity_interpolator[3].trilinear_interpol(pt_mom, phi_mom, y_mom);


            polarization_mother_shear[0] = S_shear_interpolator[0].trilinear_interpol(pt_mom, phi_mom, y_mom);
            polarization_mother_shear[1] = S_shear_interpolator[1].trilinear_interpol(pt_mom, phi_mom, y_mom);
            polarization_mother_shear[2] = S_shear_interpolator[2].trilinear_interpol(pt_mom, phi_mom, y_mom);
            polarization_mother_shear[3] = S_shear_interpolator[3].trilinear_interpol(pt_mom, phi_mom, y_mom);

            //boost to the rest frame of the mother: the inverse standard boost is given by {{e/m, -(px/m), -(py/m), -(pz/m)}, {-(px/m), 1 + px^2/(e m + m^2), (px py)/(e m + m^2), (px pz)/(e m + m^2)}, {-(py/m), (px py)/(e m + m^2), 1 + py^2/(e m + m^2), (py pz)/(e m + m^2)}, {-(pz/m), (px pz)/(e m + m^2), (py pz)/(e m + m^2), 1 + pz^2/(e m + m^2)}}
            rest_frame_Pi[0] = -((polarization_mother[0]*P_mother[1])/mother_mass) + 
                    polarization_mother[1]*(1 + pow(P_mother[1],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
                    (polarization_mother[2]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
                    (polarization_mother[3]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
                    
            rest_frame_Pi[1] = -((polarization_mother[0]*P_mother[2])/mother_mass) + 
                    polarization_mother[2]*(1 + pow(P_mother[2],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
                    (polarization_mother[1]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
                    (polarization_mother[3]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
        
            rest_frame_Pi[2] = -((polarization_mother[0]*P_mother[3])/mother_mass) + 
                    polarization_mother[3]*(1 + pow(P_mother[3],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
                    (polarization_mother[1]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
                    (polarization_mother[2]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
            
            rest_frame_PiShear[0] = -((polarization_mother_shear[0]*P_mother[1])/mother_mass) + 
                    polarization_mother_shear[1]*(1 + pow(P_mother[1],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
                    (polarization_mother_shear[2]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
                    (polarization_mother_shear[3]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
                    
            rest_frame_PiShear[1] = -((polarization_mother_shear[0]*P_mother[2])/mother_mass) + 
                    polarization_mother_shear[2]*(1 + pow(P_mother[2],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
                    (polarization_mother_shear[1]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
                    (polarization_mother_shear[3]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
        
            rest_frame_PiShear[2] = -((polarization_mother_shear[0]*P_mother[3])/mother_mass) + 
                    polarization_mother_shear[3]*(1 + pow(P_mother[3],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
                    (polarization_mother_shear[1]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
                    (polarization_mother_shear[2]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
            
			//computation of the vector S depending on the decay
            switch(mother.get_id()){
                case 3212: //Sigma0 to Lambda and photon
                for(int mu=0;mu<3;mu++){ 
                    S_vect[mu] = -(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[0]+rest_frame_Pi[1]*p_rest[1]+rest_frame_Pi[2]*p_rest[2]); 
                    S_vectShear[mu] = -(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_PiShear[0]*p_rest[0]+rest_frame_PiShear[1]*p_rest[1]+rest_frame_PiShear[2]*p_rest[2]); 
                }
                break;
                case 3224: //Sigma* to Lambda and pion
                for(int mu=0;mu<3;mu++){
                    S_vect[mu] = (2.0/5.0)*rest_frame_Pi[mu] -(1.0/5.0)*(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[0]+rest_frame_Pi[1]*p_rest[1]+rest_frame_Pi[2]*p_rest[2]); 
                    S_vectShear[mu] = (2.0/5.0)*rest_frame_PiShear[mu] -(1.0/5.0)*(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_PiShear[0]*p_rest[0]+rest_frame_PiShear[1]*p_rest[1]+rest_frame_PiShear[2]*p_rest[2]); 
                }
                break;
                default:
                cout <<"I don't know about this decay!"<< endl;
                exit(1);
            }
            
            //compute integrals in phi  
            integralphi_den += dangle*jacobian*spectrum/Energy_mother; //the n(P) in eq. 29 of 1905.03123 is no the invariant one, hence the /Energy_mother
            for(int mu=0;mu<3;mu++){
                integralphi_num[mu] += dangle*jacobian*(S_vect[mu]/spectrum) *spectrum/Energy_mother; //    S/spectrum = spin vector; spectrum/Energy_mother =  n(P) in eq 29 of 1905.03123
                integralphi_numShear[mu] += dangle*jacobian*(S_vectShear[mu]/spectrum) *spectrum/Energy_mother;
            }
        } //end loop in iph
        //compute integral in theta
        Denominator += dangle*sin(dtheta)*integralphi_den;
        for(int mu=0;mu<3;mu++){
                P_vorticity[mu] += dangle*sin(dtheta)*integralphi_num[mu];
                P_shear[mu] += dangle*sin(dtheta)*integralphi_numShear[mu];
        }
    } //end loop in ith
		

    //print to file
    fileout << "   " << pT << "   " << phi << "   " << y_rap << "   " << Denominator;
    for(int mu=0; mu<3; mu++)
        fileout << "   " << P_vorticity[mu];
    for(int mu=0; mu<3; mu++)
        fileout << "   " << P_shear[mu]; //The hbarc to make the shear adimensional is accounted for in the table
    fileout << endl;

    delete second_son;
}

// void Lambda_FeedDown_nointerpolation(double pT, double phi, double y_rap, pdg_particle mother, 
//     vector<element> &freeze_out_sup, ofstream &fileout){
//     pdg_particle lambda(3122);

//     double P_vorticity[3] = {0,0,0};
//     double P_shear[3] = {0,0,0};
//     double Denominator = 1e-20;
    
// 	const double lambda_mass = lambda.get_mass(); 
	
// 	//fetch information concerning mother and second son
// 	const double mother_mass = mother.get_mass();  
// 	const double mother_baryonCharge = mother.get_b();
// 	const double mother_electricCharge = mother.get_q();
// 	const double mother_strangeness = mother.get_s();
//     const double motherspin = mother.get_spin();
//     const int fermi_or_bosemother = statistics(motherspin);
//     const double phase_space = ((2*motherspin +1)/pow( 2*hbarC*PI, 3.0));

    
//     pdg_particle* second_son;

//     if(mother.get_id()==3212){
//         second_son = new pdg_particle(22);
//     }
//     else if(mother.get_id() == 3224){
//         second_son = new pdg_particle(211);
//     }
//     else{
//         cout<< "ERROR! The decay is not implemented yet!"<<endl;
//         exit(1);
//     }

//     double second_son_mass = second_son->get_mass();

//     double mT = sqrt(lambda_mass * lambda_mass + pT*pT);
//     //momentum of the Lambda particle in the Lab frame
//     array<double,4> p = {mT *cosh(y_rap), pT *cos(phi), pT *sin(phi), mT *sinh(y_rap)};
//     array<double,4> p_ = {mT *cosh(y_rap), -pT *cos(phi), -pT *sin(phi), -mT *sinh(y_rap)}; //lower indices  

//     //Energy and momentum of the Lambda particle in the mother's rest frame
//     double E_Lambda_rest = (mother_mass*mother_mass + lambda_mass*lambda_mass - second_son_mass*second_son_mass)/(2*mother_mass); 
//     double p_rest_abs = sqrt(E_Lambda_rest*E_Lambda_rest - lambda_mass*lambda_mass);

//     int number_of_bins = 20;
//     double dangle = PI/number_of_bins;
//     #ifdef OPEN_MP
//         int threads_ = NTHREADS;
//         #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,P_vorticity,P_shear)
//     #endif
//     for(int itheta=0; itheta<number_of_bins; itheta++){ 
//         double dtheta = itheta*dangle;
//         double integralphi_den = 1e-20;
//         double integralphi_num[3] = {0,0,0}; 
//         double integralphi_numShear[3] = {0,0,0}; 
//         for(double dphi=0; dphi<2*PI-1e-5; dphi+=dangle){
//             double rest_frame_Pi[3] = {0,0,0}; //polarization in the rest frame of the mother
//             double rest_frame_PiShear[3] = {0,0,0}; //polarization in the rest frame of the mother, shear contribution
//             double S_vect[3]={0,0,0};//the particular vector depending on the decay. For details 1905.03123
//             double S_vectShear[3]={0,0,0};//the particular vector depending on the decay, shear contribution
//             double spectrum = 1e-20;
//             double polarization_mother[4]={0,0,0,0};
//             double polarization_mother_shear[4]={0,0,0,0};
    
//             double p_rest[3] = {p_rest_abs*sin(dtheta)*cos(dphi),p_rest_abs*sin(dtheta)*sin(dphi),p_rest_abs*cos(dtheta)}; //spatial momentum of \Lambda in mother rest frame. p^0 is irrelevant to the calculation
//             //jacobian from eq. 30 of 1905.03123
//             double jacobian = pow(mother_mass,3)*pow(p[0] + E_Lambda_rest,2)*
//                                 (pow(p[0] + E_Lambda_rest,2)-(lambda_mass*lambda_mass+p[0]*E_Lambda_rest+p[1]*p_rest[0]+p[2]*p_rest[1]+p[3]*p_rest[2]))/
//                                 (E_Lambda_rest*pow(lambda_mass*lambda_mass+p[0]*E_Lambda_rest+p[1]*p_rest[0]+p[2]*p_rest[1]+p[3]*p_rest[2],3));
//             //momentum of the mother
//             double Energy_mother;
//             double P_mother[4];
//             double P_mother_[4];
//             //upper indices
//             P_mother[1] = (p[1]-p_rest[0])*2*mother_mass*(p[0] + E_Lambda_rest)/
//                         (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[0],2)-pow(p[2]-p_rest[1],2)-pow(p[3]-p_rest[2],2));
//             P_mother[2] = (p[2]-p_rest[1])*2*mother_mass*(p[0] + E_Lambda_rest)/
//                         (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[0],2)-pow(p[2]-p_rest[1],2)-pow(p[3]-p_rest[2],2));
//             P_mother[3] = (p[3]-p_rest[2])*2*mother_mass*(p[0] + E_Lambda_rest)/
//                         (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[0],2)-pow(p[2]-p_rest[1],2)-pow(p[3]-p_rest[2],2));
//             Energy_mother = sqrt(mother_mass*mother_mass+pow(P_mother[1],2)+pow(P_mother[2],2)+pow(P_mother[3],2));
//             P_mother[0] = Energy_mother; 
//             //lower indices
//             P_mother_[1] = -P_mother[1];
//             P_mother_[2] = -P_mother[2];
//             P_mother_[3] = -P_mother[3];
//             P_mother_[0] = Energy_mother;
            
//             for(element cell : freeze_out_sup){ //loop over the FO hypersurface
//                 double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)
//                 array<double,4> theta_vector = {0,0,0,0};
//                 double theta_sq = 0.;

//                 for(int mu=0; mu<4; mu++){
//                     for(int nu=0; nu<4; nu++)
//                         for(int rh=0; rh<4; rh++)
//                             for(int sg=0; sg<4; sg++){ 
//                                 theta_vector[mu] += levi(mu, nu, rh, sg)
//                                         * P_mother_[sg] * cell.dbeta[nu][rh]/(2*mother_mass); 
//                 }
//                 theta_vector[mu] *= hbarC; //now theta is adimensional 
//                 }

//                 for (int mu = 0; mu < 4; mu++){
//                     pdSigma += P_mother[mu] * cell.dsigma[mu];
//                     pu += P_mother[mu] * cell.u[mu] * gmumu[mu];
//                     theta_sq += theta_vector[mu]*theta_vector[mu]*gmumu[mu];
//                 }
//                 const double mutot = cell.mub*mother_baryonCharge + cell.muq*mother_electricCharge + cell.mus*mother_strangeness;
//                 const double distribution = 1 / (exp( (pu - mutot) / cell.T) + fermi_or_bosemother);

//                 spectrum += phase_space*pdSigma * distribution ;
                
//                 for(int mu=0; mu<4; mu++){
//                     // computing the vorticity induced polarization
//                     polarization_mother[mu] += phase_space*pdSigma * distribution *  
//                                 (theta_vector[mu]/sqrt(-theta_sq)) * 
//                                 aux_exact_polarization(motherspin, pu, cell.T, mutot, sqrt(-theta_sq));

//                     // computing the shear induced polarization
//                     for(int sg=1; sg<4; sg++)
//                         for(int ta=1; ta<4; ta++)
//                             for(int rh=0; rh<4; rh++){
//                             polarization_mother_shear[mu] += - phase_space*(motherspin/3)*(motherspin+1)*pdSigma * distribution * (1. - fermi_or_bosemother*distribution) 
//                                         * levi(mu, 0, sg, ta)* (P_mother_[ta]/mother_mass) * (P_mother[rh]/P_mother[0]) 
//                                         *hbarC*( cell.dbeta[sg][rh] + cell.dbeta[rh][sg])/ 2.0;
//                                         //hbarC: unit conversion to make the shear adimensional 
//                             }
//                 }
//             }
            
//             //boost to the rest frame of the mother: the inverse standard boost is given by {{e/m, -(px/m), -(py/m), -(pz/m)}, {-(px/m), 1 + px^2/(e m + m^2), (px py)/(e m + m^2), (px pz)/(e m + m^2)}, {-(py/m), (px py)/(e m + m^2), 1 + py^2/(e m + m^2), (py pz)/(e m + m^2)}, {-(pz/m), (px pz)/(e m + m^2), (py pz)/(e m + m^2), 1 + pz^2/(e m + m^2)}}
//             rest_frame_Pi[0] = -((polarization_mother[0]*P_mother[1])/mother_mass) + 
//                     polarization_mother[1]*(1 + pow(P_mother[1],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
//                     (polarization_mother[2]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
//                     (polarization_mother[3]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
                    
//             rest_frame_Pi[1] = -((polarization_mother[0]*P_mother[2])/mother_mass) + 
//                     polarization_mother[2]*(1 + pow(P_mother[2],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
//                     (polarization_mother[1]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
//                     (polarization_mother[3]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
        
//             rest_frame_Pi[2] = -((polarization_mother[0]*P_mother[3])/mother_mass) + 
//                     polarization_mother[3]*(1 + pow(P_mother[3],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
//                     (polarization_mother[1]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
//                     (polarization_mother[2]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
            
//             rest_frame_PiShear[0] = -((polarization_mother_shear[0]*P_mother[1])/mother_mass) + 
//                     polarization_mother_shear[1]*(1 + pow(P_mother[1],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
//                     (polarization_mother_shear[2]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
//                     (polarization_mother_shear[3]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
                    
//             rest_frame_PiShear[1] = -((polarization_mother_shear[0]*P_mother[2])/mother_mass) + 
//                     polarization_mother_shear[2]*(1 + pow(P_mother[2],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
//                     (polarization_mother_shear[1]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
//                     (polarization_mother_shear[3]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
        
//             rest_frame_PiShear[2] = -((polarization_mother_shear[0]*P_mother[3])/mother_mass) + 
//                     polarization_mother_shear[3]*(1 + pow(P_mother[3],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
//                     (polarization_mother_shear[1]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
//                     (polarization_mother_shear[2]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
            
// 			//computation of the vector S depending on the decay
//             switch(mother.get_id()){
//                 case 3212: //Sigma0 to Lambda and photon
//                 for(int mu=0;mu<3;mu++){ 
//                     S_vect[mu] =      -(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_Pi[0]     *p_rest[0]+rest_frame_Pi[1]     *p_rest[1]+rest_frame_Pi[2]     *p_rest[2]); 
//                     S_vectShear[mu] = -(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_PiShear[0]*p_rest[0]+rest_frame_PiShear[1]*p_rest[1]+rest_frame_PiShear[2]*p_rest[2]); 
//                 }
//                 break;
//                 case 3224: //Sigma* to Lambda and pion
//                 for(int mu=0;mu<3;mu++){
//                     S_vect[mu]      = (2.0/5.0)*rest_frame_Pi[mu]      -(1.0/5.0)*(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_Pi[0]     *p_rest[0]+rest_frame_Pi[1]     *p_rest[1]+rest_frame_Pi[2]     *p_rest[2]); 
//                     S_vectShear[mu] = (2.0/5.0)*rest_frame_PiShear[mu] -(1.0/5.0)*(p_rest[mu]/pow(p_rest_abs,2))*(rest_frame_PiShear[0]*p_rest[0]+rest_frame_PiShear[1]*p_rest[1]+rest_frame_PiShear[2]*p_rest[2]); 
//                 }
//                 break;
//                 default:
//                 cout <<"I don't know about this decay!"<< endl;
//                 exit(1);
//             }
            
//             //compute integrals in phi  
//             integralphi_den += dangle*spectrum*jacobian;
//             for(int mu=0;mu<3;mu++){
//                 integralphi_num[mu] += dangle*jacobian*S_vect[mu]; //the spectrum is not present because it simplify with the denominator in the polarization formula
//                 integralphi_numShear[mu] += dangle*jacobian*S_vectShear[mu];
//             }
//         } //end loop in iph
//         //compute integral in theta
//         Denominator += dangle*sin(dtheta)*integralphi_den;
//         for(int mu=0;mu<3;mu++){
//                 P_vorticity[mu] += dangle*sin(dtheta)*integralphi_num[mu];
//                 P_shear[mu] += dangle*sin(dtheta)*integralphi_numShear[mu];
//         }
//     } //end loop in ith
		

//     //print to file
//     fileout << "   " << pT << "   " << phi << "   " << y_rap << "   " << Denominator;
//     for(int mu=0; mu<3; mu++)
//         fileout << "   " << P_vorticity[mu];
//     for(int mu=0; mu<3; mu++)
//         fileout << "   " << P_shear[mu]; //The hbarc to make the shear adimensional is accounted for in the table
//     fileout << endl;

//     delete second_son;
// }
