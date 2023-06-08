#include <cmath>
#include <fstream>
#include <array>


#ifdef OPEN_MP
    #include<omp.h>
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
    array<double,4> p = {mT, pT *cos(phi), pT *sin(phi), 0};
    array<double,4> p_ = {mT, -pT *cos(phi), -pT *sin(phi), 0}; //lower indices    

    #ifdef OPEN_MP
        int threads_ = omp_get_max_threads() - 2; 
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
        int threads_ = omp_get_max_threads() - 2; 
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


void shear_correction_to_spectra(double pT, double phi, double Dt, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout){
    //Dt is the time in fermi

    double dNd3p = 0;
    double dNd3p_shear = 0;

    const double mass = particle.get_mass();  // get particle's info
    const int baryonNumber = particle.get_b();
    const int electricCharge = particle.get_q();
    const int strangeness = particle.get_s();
    
    const int dim_spin = (2*particle.get_spin()+1); //Dimension of the spin Hilbert space: if even -> fermions, if odd -> bosons
    int fermi_or_bose = 1; //the factor to add in the denominator of the distribution: 1 Fermi-Dirac, -1 Bose-Einstein
    if(dim_spin % 2){
        fermi_or_bose *= -1;
    }

    #ifdef OPEN_MP
        int threads_ = omp_get_max_threads() - 2; 
        #pragma omp parallel for num_threads(threads_) reduction(+:dNd3p,dNd3p_shear)
    #endif
    for(element cell : freeze_out_sup){ //loop over the FO hypersurface
        double mT = sqrt(mass * mass + pT*pT);
        array<double,4> p = {mT, pT *cos(phi), pT *sin(phi), 0};
        array<double,4> p_ = {mT, -pT *cos(phi), -pT *sin(phi), 0}; //lower indices
        double pdSigma = 0., pu = 0.;  //scalar products p\cdot d\Sigma and p\cdot u (u is the four velocity)
        double pp_dbeta = 0, ept_dbeta = 0; //the two corrections going in the shear formula p_rho p^sigma:(g-tt)^{rho lambda}dbeta_{lambda sigma}= p^sigma(p^lambda - e t^lambda)dbeta_{lambda sigma}

        for (int mu = 0; mu < 4; mu++) {
            pdSigma += p[mu] * cell.dsigma[mu];
            pu += p[mu] * cell.u[mu] * gmumu[mu];
            ept_dbeta += p[0]*p[mu]*cell.dbeta[0][mu]; //there should be a t^nu which select zero in the first entry
        }

        for(int mu=0; mu<4; mu++){
            for(int nu=0; nu<4; nu++){
                pp_dbeta += p[mu]*p[nu]*cell.dbeta[mu][nu];
            }
        }
        //convert dbeta to 1/(Gev*fm) units
        ept_dbeta *= hbarC;
        pp_dbeta *= hbarC;

        const double mutot = cell.mub*baryonNumber + cell.muq*electricCharge + cell.mus*strangeness;
        const double distrib = 1 / (exp( (pu - mutot) / cell.T) + fermi_or_bose);

        dNd3p += pdSigma * distrib /pow(2*PI,3);
        dNd3p_shear += (pdSigma/pow(2*PI,3)) *distrib*(1+(1-fermi_or_bose*distrib)*(Dt/p[0])*(pp_dbeta-ept_dbeta)); //The correction now is a pure number: fm/GeV *(Gev^2*(1/(GeV *fm))
    }                                                                                                               //                                      Dt/p  *( pp      *dbeta)

    dNd3p /= hbarC*hbarC*hbarC;
    dNd3p_shear /= hbarC*hbarC*hbarC;
    
    //print to file
    fileout << "   " << pT << "   " << phi << "   " << dNd3p << "   " << dNd3p_shear <<endl;
 
}
