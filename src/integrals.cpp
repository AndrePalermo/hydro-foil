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
    double Denominator = 0;

    // get particle's info
    const double mass = particle.get_mass();  //GeV
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();

    double mT = sqrt(mass * mass + pT*pT);
    std::array<double,4> p = {mT, pT *cos(phi), pT *sin(phi), 0};
    std::array<double,4> p_ = {mT, -pT *cos(phi), -pT *sin(phi), 0}; //lower indices    

    #ifdef OPEN_MP
        int threads_ = NTHREADS; 
        #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,P_vorticity,P_shear)
    #endif
    for(element cell : freeze_out_sup){ //loop over the FO hypersurface
        double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)
        array<double,4> theta_vector = {0,0,0,0};
        double theta_sq = 0.;

        for(int mu=0; mu<4; mu++){
            for(int nu=0; nu<4; nu++)
                for(int rh=0; rh<4; rh++)
                    for(int sg=0; sg<3; sg++){ //sg=3 is zero because p_[3]=0. This speeds up the program
		                theta_vector[mu] += levi(mu, nu, rh, sg)
                                * p_[sg] * cell.dbeta[nu][rh]/(2*mass); 
		   }
           theta_vector[mu] *= hbarC; //now theta is adimensional 
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

            for(int rh=0; rh<4; rh++)
                for(int sg=0; sg<3; sg++) {  //since we are at midrapidity p_[3]=0 and we can stop at sg=2
                    // computing the shear induced polarization

                    for(int ta=0; ta<4; ta++)
                    P_shear[mu] += - pdSigma * nf * (1. - nf) * levi(mu, 0, rh, sg)* p_[sg] * p[ta] / p[0] 
                                * ( cell.dbeta[rh][ta] + cell.dbeta[ta][rh])/ (8.0 * mass);
                }
        }
    }

    //print to file
    fileout << "   " << pT << "   " << phi << "   " << Denominator;
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_vorticity[mu];
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_shear[mu] *hbarC; //Unit conversion to make the shear adimensional 
    fileout << endl;

}

void polarization_projected(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
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
        double projector[4][4], projected_gradient[4][4];
        double normal_size=0; //takes into account the normalization of the normal vector

        for(int mu = 0; mu < 4; mu++){
            normal_size += cell.dsigma[mu]*cell.dsigma[mu]*gmumu[mu];
        }
        

        for(int mu = 0; mu < 4; mu++){
            for(int nu = 0; nu < 4; nu++){
                projector[mu][nu] = g(mu,nu) - cell.dsigma[mu]*cell.dsigma[nu]/(normal_size); //projector with lower indices
                projected_gradient[mu][nu] = 0;
            }
        } 

        for(int mu = 0; mu < 4; mu++){
            for(int nu = 0; nu < 4; nu++){
                for(int alpha = 0; alpha < 4; alpha++){
                    for(int sigma = 0; sigma < 4; sigma++){
                        projected_gradient[mu][nu] += projector[mu][alpha]*g(alpha,sigma)*cell.dbeta[sigma][nu]; //lower indices
                        //NOTE: cell.dbeta[sigma][nu] = \partial_\sigma \beta_\nu
                    }
                }
            }
        }    

        for (int mu = 0; mu < 4; mu++) {
            pdSigma += p[mu] * cell.dsigma[mu];
            pu += p[mu] * cell.u[mu] * gmumu[mu];
        }
        const double mutot = cell.mub*baryonNumber + cell.muq*electricCharge + cell.mus*strangeness;
        const double nf = 1 / (exp( (pu - mutot) / cell.T) + 1.0);

        Denominator += pdSigma * nf ;
        for(int mu=0; mu<4; mu++)
            for(int nu=0; nu<4; nu++)
                for(int rh=0; rh<4; rh++)
                    for(int sg=0; sg<4; sg++) {
                        
                        P_vorticity[mu] += pdSigma * nf * (1. - nf) * levi(mu, nu, rh, sg)
                                                * p_[sg] * projected_gradient[nu][rh];
                        
                        if(nu==0)
                        for(int ta=0; ta<4; ta++)
                        P_shear[mu] += -pdSigma * nf * (1. - nf) * levi(mu, nu, rh, sg)
                                    * p_[sg] * (p[ta] / p[0]) *t_vector[nu]
                                    * ( projected_gradient[rh][ta] + projected_gradient[ta][rh]);
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
    std::array<double,4> p = {mT, pT *cos(phi), pT *sin(phi), 0};
    std::array<double,4> p_ = {mT, -pT *cos(phi), -pT *sin(phi), 0}; //lower indices    

     

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
        for(int mu=0; mu<4; mu++)
            for(int nu=0; nu<4; nu++)
                for(int rh=0; rh<4; rh++)
                    for(int sg=0; sg<4; sg++) {
                        
                        P_vorticity[mu] += pdSigma * nf * (1. - nf) * levi(mu, nu, rh, sg)
                                                * p_[sg] * cell.dbeta[nu][rh];
                        
                        if(nu==0)
                        for(int ta=0; ta<4; ta++)
                        P_shear[mu] += -pdSigma * nf * (1. - nf) * levi(mu, nu, rh, sg)
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

void polarization_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
    double P_vorticity[4] = {0,0,0,0};
    double P_shear[4] = {0,0,0,0};
    double P_shear_chun[4] = {0,0,0,0};
    double Denominator = 0;

    // get particle's info
    const double mass = particle.get_mass();  //GeV
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();

    double mT = sqrt(mass * mass + pT*pT);
    array<double,4> p = {mT *cosh(y_rap), pT *cos(phi), pT *sin(phi), mT *sinh(y_rap)};
    array<double,4> p_ = {mT *cosh(y_rap), -pT *cos(phi), -pT *sin(phi), -mT *sinh(y_rap)}; //lower indices    

    #ifdef OPEN_MP
        int threads_ = NTHREADS; 
        #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,P_vorticity,P_shear,P_shear_chun)
    #endif
    for(element cell : freeze_out_sup){ //loop over the FO hypersurface
        double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)
        array<double,4> theta_vector = {0,0,0,0};
        double theta_sq = 0.;

        for(int mu=0; mu<4; mu++){
            for(int nu=0; nu<4; nu++)
                for(int rh=0; rh<4; rh++)
                    for(int sg=0; sg<4; sg++){ 
		                theta_vector[mu] += levi(mu, nu, rh, sg)
                                * p_[sg] * cell.dbeta[nu][rh]/(2*mass); 
		   }
           theta_vector[mu] *= hbarC; //now theta is adimensional 
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

            // computing the shear induced polarization
            for(int rh=0; rh<4; rh++)
                for(int sg=0; sg<4; sg++)
                    for(int ta=0; ta<4; ta++){
                    P_shear[mu] += - pdSigma * nf * (1. - nf) * levi(mu, 0, rh, sg)* p_[sg] * p[ta] / p[0] 
                                * ( cell.dbeta[rh][ta] + cell.dbeta[ta][rh])/ (8.0 * mass);
                    P_shear_chun[mu] += - pdSigma * nf * (1. - nf) * levi(mu, 0, rh, sg)* p_[sg] * p[ta] / pu 
                                * ( cell.dbeta[rh][ta] + cell.dbeta[ta][rh])/ (8.0 * mass);
                    }
        }
    }

    //print to file
    fileout << "   " << pT << "   " << phi << "   " << y_rap << "   " << Denominator;
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_vorticity[mu];
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_shear[mu] *hbarC; //Unit conversion to make the shear adimensional 
    for(int mu=0; mu<4; mu++)
        fileout << "   " << P_shear_chun[mu] *hbarC; //Unit conversion to make the shear adimensional 
    fileout << endl;

}

void spectrum_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
    double dNd3p = 0;
    
    // get particle's info
    const double mass = particle.get_mass();  
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();
    
    const int dim_spin = (2*particle.get_spin()+1); //Dimension of the spin Hilbert space: if even -> fermions, if odd -> bosons
    int fermi_or_bose = 1; //the factor to add in the denominator of the distribution: 1 Fermi-Dirac, -1 Bose-Einstein
    if(dim_spin % 2){
        fermi_or_bose *= -1;
    }

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

void Lambda_polarization_FeedDown(double pT, double phi, double y_rap, pdg_particle mother, vector<element> &freeze_out_sup, ofstream &fileout){
    pdg_particle lambda(3122);

    double P_vorticity[4] = {0,0,0,0};
    double P_shear[4] = {0,0,0,0};
    double Denominator = 0;
    
	const double mass = lambda.get_mass(); 
	
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

    std::cout<<"calculations for the decay: "<< mother.get_name()<< " to "<< lambda.get_name() << " and "<< second_son->get_name()<<endl;
    double second_son_mass = second_son->get_mass();

    double mT = sqrt(mass * mass + pT*pT);
    //momentum of the Lambda particle in the Lab frame
    array<double,4> p = {mT *cosh(y_rap), pT *cos(phi), pT *sin(phi), mT *sinh(y_rap)};
    array<double,4> p_ = {mT *cosh(y_rap), -pT *cos(phi), -pT *sin(phi), -mT *sinh(y_rap)}; //lower indices  

    //Energy and momentum of the Lambda particle in the mother's rest frame
    double E_Lambda_rest = (mother_mass*mother_mass + mass*mass - second_son_mass*second_son_mass)/(2*mother_mass); 
    double p_rest_abs = sqrt(E_Lambda_rest*E_Lambda_rest - mass*mass);

    double dangle = PI/20;
    for(double dtheta=0; dtheta<PI-1e-5; dtheta+=PI/20){
        double integralphi_den = 0;
        double integralphi_num[3] = {0,0,0}; 
        double integralphi_numShear[3] = {0,0,0}; 
        for(double dphi=0; dphi<2*PI-1e-5; dphi+=PI/20){
            double rest_frame_Pi[3] = {0,0,0}; //polarization in the rest frame of the mother
            double rest_frame_PiShear[3] = {0,0,0}; //polarization in the rest frame of the mother, shear contribution
            double S_vect[3]={0,0,0};//the particular vector depending on the decay. For details 1905.03123
            double S_vectShear[3]={0,0,0};//the particular vector depending on the decay, shear contribution
            double spectrum = 0;
            double polarization_mother[4]={0,0,0,0};
            double polarization_mother_shear[4]={0,0,0,0};
    
            double p_rest[4] = {0,p_rest_abs*sin(dtheta)*cos(dphi),p_rest_abs*sin(dtheta)*sin(dphi),p_rest_abs*cos(dtheta)};
            //jacobian from eq. 30 of 1905.03123
            double jacobian = pow(mother_mass,3)*pow(p[0] + E_Lambda_rest,2)*
                                (pow(p[0] + E_Lambda_rest,2)-(mass*mass+p[0]*E_Lambda_rest+p[1]*p_rest[1]+p[2]*p_rest[2]+p[3]*p_rest[3]))/
                                (E_Lambda_rest*pow(mass*mass+p[0]*E_Lambda_rest+p[1]*p_rest[1]+p[2]*p_rest[2]+p[3]*p_rest[3],3));
            //momentum of the mother
            double Energy_mother;
            double P_mother[4];
            double P_mother_[4];
            //upper indices
            P_mother[1] = (p[1]-p_rest[1])*2*mother_mass*(p[0] + E_Lambda_rest)/
                        (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
            P_mother[2] = (p[2]-p_rest[2])*2*mother_mass*(p[0] + E_Lambda_rest)/
                        (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
            P_mother[3] = (p[3]-p_rest[3])*2*mother_mass*(p[0] + E_Lambda_rest)/
                        (pow(p[0] + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
            Energy_mother = sqrt(mother_mass*mother_mass+pow(P_mother[1],2)+pow(P_mother[2],2)+pow(P_mother[3],2));
            P_mother[0] = Energy_mother; 
            //lower indices
            P_mother_[1] = -P_mother[1];
            P_mother_[2] = -P_mother[2];
            P_mother_[3] = -P_mother[3];
            P_mother_[0] = Energy_mother;
            

            #ifdef OPEN_MP
                int threads_ = NTHREADS; 
                #pragma omp parallel for num_threads(threads_) reduction(+:Denominator,polarization_mother,polarization_mother_shear)
            #endif
            for(element cell : freeze_out_sup){ //loop over the FO hypersurface
                double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)

                array<double,4> theta_vector = {0,0,0,0};
                double theta_sq = 0.;
                
                for(int mu=0; mu<4; mu++){
                    for(int nu=0; nu<4; nu++)
                        for(int rh=0; rh<4; rh++)
                            for(int sg=0; sg<4; sg++){ 
                                theta_vector[mu] += levi(mu, nu, rh, sg)
                                        * P_mother_[sg] * cell.dbeta[nu][rh]/(2*mother_mass); 
                            }
                    theta_vector[mu] *= hbarC; //now theta is adimensional 
                }

                for (int mu = 0; mu < 4; mu++) {
                    pdSigma += P_mother[mu] * cell.dsigma[mu];
                    pu += P_mother[mu] * cell.u[mu] * gmumu[mu];
                    theta_sq += theta_vector[mu]*theta_vector[mu]*gmumu[mu];
                }
                double mutot = cell.mub * mother_baryonCharge
                                        + cell.muq * mother_electricCharge + cell.mus * mother_strangeness;
                double nf =  1 / (exp( (pu - mutot) / cell.T) + 1.0);
                spectrum += pdSigma * nf ;
                for(int mu=0; mu<4; mu++){
                    // computing the vorticity induced polarization
                    polarization_mother[mu] += 0.5*pdSigma * nf *  
                                (theta_vector[mu]/sqrt(-theta_sq)) * sinh(sqrt(-theta_sq)*0.5)/
                                (cosh(sqrt(-theta_sq)*0.5)+exp(-(pu - mutot)/cell.T));

                    for(int rh=0; rh<4; rh++)
                        for(int sg=0; sg<4; sg++) { 
                            // computing the shear induced polarization
                            for(int ta=0; ta<4; ta++)
                            polarization_mother_shear[mu] += - pdSigma * nf * (1. - nf) * levi(mu, 0, rh, sg)* P_mother_[sg] * P_mother[ta] / P_mother[0] 
                                        * ( cell.dbeta[rh][ta] + cell.dbeta[ta][rh])/ (8.0 * mass);
                        }
                } 
            }

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
                    S_vect[mu] = -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[1]+rest_frame_Pi[1]*p_rest[2]+rest_frame_Pi[2]*p_rest[3]); 
                    S_vectShear[mu] = -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_PiShear[0]*p_rest[1]+rest_frame_PiShear[1]*p_rest[2]+rest_frame_PiShear[2]*p_rest[3]); 
                }
                break;
                case 3224: //Sigma* to Lambda and pion
                for(int mu=0;mu<3;mu++){
                    S_vect[mu] =2*rest_frame_Pi[mu] -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[1]+rest_frame_Pi[1]*p_rest[2]+rest_frame_Pi[2]*p_rest[3]); 
                    S_vectShear[mu] =2*rest_frame_PiShear[mu] -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_PiShear[0]*p_rest[1]+rest_frame_PiShear[1]*p_rest[2]+rest_frame_PiShear[2]*p_rest[3]); 
                }
                break;
                default:
                cout <<"I don't know about this decay!"<< endl;
                exit(1);
            }
            
            //compute integrals in phi  
            integralphi_den += dangle*spectrum*jacobian/P_mother[0];
            for(int mu=0;mu<3;mu++){
                integralphi_num[mu] += dangle*jacobian*S_vect[mu]/P_mother[0];
                integralphi_numShear[mu] += dangle*jacobian*S_vectShear[mu]/P_mother[0];
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
        fileout << "   " << P_shear[mu] *hbarC; //Unit conversion to make the shear adimensional 
    fileout << endl;
}

