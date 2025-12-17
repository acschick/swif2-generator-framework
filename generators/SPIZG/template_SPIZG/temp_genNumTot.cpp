//#include "FF2.h"
#include "TCanvas.h"                                                                                                                                                
#include "TROOT.h"                                                                                                                            
#include "TGraphErrors.h"                                                                                                                
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TAxis.h"                                                                                                                                                     
#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include<gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

//extern "C" {
//  void coul_(double* beamE, double* theta, double* w1, double* w2);
//}

/*
double FFscr(double ebeam, double theta){
  double mpi = .1349;

  theta = theta*TMath::Pi()/180;

  double w1, w2;
  formfactorc_(&ebeam, &theta, &w1, &w2);

  return w1*w1 + w2*w2;

*/
double FF(double q, int ztgt){
  //The following was translated from FORTRAN code written by Rory Miskimen for the electromagnetic form factor for scattering off Pb or a proton target
  //This is not currently in use, and does not take into account absorbtion. Can be turned on for proton target
  double q02, hbarc, qF, gamma, rho0, pi, norm, proton_rms;

  vector<double> R, A;

  q02 = 0.6608;
  hbarc = .197326;
  
  R.push_back(0.1);
  R.push_back(0.7);
  R.push_back(1.6);
  R.push_back(2.1);
  R.push_back(2.7);
  R.push_back(3.5);
  R.push_back(4.2);
  R.push_back(5.1);
  R.push_back(6.0);
  R.push_back(6.6);
  R.push_back(7.6);
  R.push_back(8.7);

  A.push_back(0.003845);
  A.push_back(0.009724);
  A.push_back(0.033093);
  A.push_back(0.000120);
  A.push_back(0.083107);
  A.push_back(0.080869);
  A.push_back(0.139957);
  A.push_back(0.260892);
  A.push_back(0.336013);
  A.push_back(0.033637);
  A.push_back(0.018729);
  A.push_back(0.000020);

  pi = 3.14159;
  gamma = 1.388;
  
  if(ztgt == 1){
    return 1/(1 + q*q/q02)*1/(1 + q*q/q02);
  }
  else if(ztgt == 82){
    qF = q/hbarc;
    double FF = 0;
    for (int i = 0; i < 12; i++){
      FF += A[i]*(gamma*gamma*cos(qF*R[i]) + 2*R[i]/qF*sin(qF*R[i]))/(gamma*gamma + 2*R[i]*R[i])*exp(-qF*qF*gamma*gamma/4);
      //cout << "FF: " << FF << endl;
      //cout << "2*R[i]/qF*sin(qF*R[i]))" << 2*R[i]/qF*sin(qF*R[i]) << endl;
      //cout << "qF: " << qF << endl;
      //cout << "-qF*qF*gamma*gamma/4: " << -qF*qF*gamma*gamma/4 << endl;
      //cout << "exp(-qF*qF*gamma*gamma/4): " << exp(-qF*qF*gamma*gamma/4) << endl;
      //cout << "cos(qf*R[i]): " << cos(qF*R[i]) << endl;
    }

    return FF;
    
  }
  else{
    return -999;
  }
}

// Global GSL random number generator (initialized once in main)
gsl_rng * global_rng = NULL;

double rand_num(){
  // Use the global RNG - much more efficient than creating/destroying on each call
  // Returns uniform random number between 0 and 1
  return gsl_rng_uniform(global_rng);
}
double kin2mom(double beamE, double theta){
  //This function takes in the incident photon energy and production angle of the pi0, and returns the magnitude of the pi0 momentum
  //This quadratic was gotten by working out the kinematics for exclusive pi0 photoproduction off Pb
  const double mpi0 = .134978;
  const double Mt = 193.688;
  
  double C = (4*beamE*beamE*Mt*Mt - 4*beamE*beamE*mpi0*mpi0 - 4*beamE*Mt*mpi0*mpi0 - 4*Mt*Mt*mpi0*mpi0 + mpi0*mpi0*mpi0*mpi0);
  double A = 4*beamE*beamE*cos(theta*TMath::Pi()/180)*cos(theta*TMath::Pi()/180) - 4*beamE*beamE - 8*beamE*Mt - 4*Mt*Mt;
  double B = beamE*cos(theta)*(8*beamE*Mt + 4*mpi0*mpi0);

  return (-B-sqrt(B*B - 4*A*C))/(2*A);
}
double G(double phi, double P, double phi0){
  //Normalized integral of NPP equation for phi with polarization 
  return P*sin(phi)*cos(phi-2*phi0)/(2*TMath::Pi()) + phi/(2*TMath::Pi());
}
double Gmin(double y, double P, double phi0){
  //There is no trivial solution to the inverse of the inegral of the NPP equation (to be used for production a probability density
  //function for producing phi values). So, for this generator I went with Newton's method
  double phi = y;

  //cout << "Got here too!\n";
  
  double F, Fp, dPhi;
  for (int i = 0; i < 500; i++){
    //cout << "here too!!\n";
    F = G(phi, P, phi0) - y; //Best Linear approximation of the desired
    Fp = 1/(2*TMath::Pi())*(1 + P*cos(2*(phi - phi0))); //Derivative of F
    //cout << "Fp: " << Fp;

    //Newton's method is applied below
    dPhi = -1*F/Fp; 
    phi += dPhi;
    if(abs(dPhi) < pow(10, -8)){ //once the value of phi is within 10^(-8) of the desirev value, return with that value 
      return phi*180/TMath::Pi() - 180;
    }
    
  }

  /*
  double phi2 = phi;

  double F2, Fp2, dPhi2;
  for (int i = 0; i < 100; i++){
    //cout << "here too!!\n";                                                                                                                                                                                                                            
    F2 = G(phi, P, phi0) - y;
    Fp2 = 1/(2*TMath::Pi())*(1 + P*cos(2*(phi - phi0)));
    //cout << "Fp: " << Fp;                                                                                                                                                                                                                              

    dPhi2 = F2/Fp2;
    phi2 += dPhi2;
    if(abs(dPhi2) < pow(10, -8)){
      return phi2*180/TMath::Pi() - 180;
    }

  }
  */
  return -999;
  
}
void genNumTot(gsl_rng* rng, int nevents){
   //This is the meat of the generator. It's general algorithm is as follows:
  // NOTE: RNG pointer passed in for proper random number generation
  // (1) Produce cross-section table for 1000 bins in theta (production lab angle) and 100 bins in incident photon energy using form factor tables
  // (2) Produce a random number between 0 and 1 degree for theta, and from the minimum value and maximumum value
  // of the provided beam energy distirbution for incident photon energy
  // (2) Perform accept/reject on beam energy values using provided beam energy distribution
  // (3) Perform two-parameter accept/reject on production angle and beam energy
  // (4) For each beam energy + theta combination, produce a random number between -180 and 180 degrees with the probability density function of the
  // NPP equation for phi. This is done by finding the inverse of the integral from 0 to phi (normalized by the integral from 0 to 2pi) of the NPP
  // equation, and feeding a random number between 0 and 1 into that. There is no analytic inverse to this, but newton's method provides good results
  // for polarization fractions less than .995 or so
  // (5) From generated E, theta, and phi values, determine the pi0 and recoil Pb four momenta and fill a CSV file with that
  // The end result is a CSV file that can then be converted to an hddm file (either with or without the recoil Pb, since the velocity of the Pb is very
  // small). This can then be thrown at the detectors using gluexMC.py

  //definint params:
  int SPIZG_NUMEVENTS = nevents;
  int MAXEVENTS = SPIZG_NUMEVENTS*1000;
  double SPIZG_POL_FRAC = 1.0;
  double SPIZG_POL_DEG = 135.0;
  string SPIZG_FFFILE = "/w/halld-scshelf2101/home/shannen/events/generators/FF_inter/LargeTableCoulBigBig2.csv";
  bool SPIZG_PRIM = true;
  bool SPIZG_STRO = true;
  bool SPIZG_INTER = true;
  double SPIZG_PHI = 45;
  double SPIZG_GAMMA = 7.806;
  double SPIZG_C = 500;
  bool SPIZG_ISSSIN = true;
  string SPIZP_OUT = "./";
  
  const double Mt = 193.688;
  
  const double mpi0 = .134978;

  //string SPIZP_OUT = "."
  //File containing beam energy distribution used for accept/reject, from Andrew Schick
  TFile *f1 = TFile::Open("/work/halld/home/acschick/hdsim_cobrems/BGRate_cpp/BGRate_CPP_merged.root");

  //char root_out = SPIZP_OUT + "/monteInputHistPrim.root";
  //output file
  TFile *fout = TFile::Open("vectors/monteInputHistImTest2.root", "RECREATE");
  
  if (!f1) { 
    cerr << "ERROR: Could not open beam energy file" << endl;
    return; 
  }
  if (!fout || fout->IsZombie()) { 
    cerr << "ERROR: Could not create output ROOT file" << endl;
    return; 
  }

  //Gamma is radiative width of pi0
  const double Gamma = SPIZG_GAMMA*pow(10, -9);
  //Z and alpha are the atomic number and fine structure constant
  int Z = 82;
  const double alph = .00072621641;
  
  
  //cout << rand_num() << endl;

  //Getting beam energy distribution from Andrew's root file
  TH1D * beamEdist = (TH1D*)f1->Get("dRtdkH1");

  int maxBin = beamEdist->GetMaximumBin(); //Beam energy bin with the highest counts
  double dNdEmax = beamEdist->GetBinContent(maxBin); //Getting counts from that bin
  cout << dNdEmax << endl;
  int numBin = beamEdist->GetNbinsX(); //Getting the total number of bins in the beam energy distribution 
  //cout << numBin << endl;
  double minE = beamEdist->GetBinLowEdge(1); //Get the lowest value of energy in distribution  
  //cout << minE << endl;
  double binWidth = beamEdist->GetBinCenter(1) - minE; //Bin width for beam energy bin  
  double maxE = beamEdist->GetBinCenter(numBin) + binWidth; //Use bin width to calculate the maximum energy in distribution (root doesn't have such a function
  //cout << maxE << endl;

  TH1D * beamEdistNorm = new TH1D("beamEdistNorm", "", numBin, minE, maxE);
  beamEdistNorm->Add(beamEdist, 1/beamEdist->Integral(1, numBin));
  
  //  double randE = rand_num(0)*(maxE - minE) + minE;
  //cout << maxE << endl;
  double randE2, dNdE;
  double randPhi, randTheta;
  double randX, randY;

  //Histogram declarations
  //double P = .7;
  int loops = 1000;
  TH1D * Eg = new TH1D("Eg", "Incident photon energy", numBin*2, minE, maxE);
  TH1D * thetapi0 = new TH1D("thetapi0", "pi0 angle", 100, 0, 2);
  //TH1D * randX_hist = new TH1D("randX_hist", "random X", 100, .995, 1);
  TH1D * phiPi0 = new TH1D("randPhi", "pi0 phi", 73, -182.5, 182.5);
  //TH1D * randNum_test = new TH1D("randNum_test", "randNUm test", 100, 0, 1);
  TH1D * FFcoul = new TH1D("FFcoul", "EM Form Factor", 100, 0, 1);
  TH1D * FFcoul2 = new TH1D("FFcoul2", "EM Form Factor 2", 100, 0, 1);
  //TH1D * primHist = new TH1D("prim", "primakoff cross-section", numBin*30, 0, 4);
  TH2D * primTot = new TH2D("primTot", "Primakoff dist", 4000, 0, 4, 100, minE, maxE);
  TH2D * interTot = new TH2D("interTot", "Interference cross-section dist for acc/rej", 4000, 0, 4, 100, minE, maxE);
  TH2D * stroTot = new TH2D("stroTot", "Strong dist", 4000, 0, 4, 100, minE, maxE);
  TH2D * toUTot = new TH2D("toUTot", "Dist to use", 4000, 0, 4, 100, minE, maxE);
  TH1D * interTotTheta = new TH1D("interTotTheta", "interferecne cross-section in theta", 4000, 0, 4);
  TH2D * FFvThet2 = new TH2D("FFvThet2", "Form Factor vs theta", 1000, 0, 1, 1000, 0, 1);
  TH2D * FFvThet = new TH2D("FFvThet", "Form Factor vs theta", 1000, 0, 1, 1000, 0, 1);

  //Generator output file
  std::ofstream myfile;
  myfile.open("vectors/genPrim.csv");
  //myfile.open ("genPrim.csv");

  //Below is the code for using the form factor table derived from the form factor scripts written by Ilya Larin. 
  //The scripts take in beam energy and production angle and give the form factor. I produced 10 different values for
  //the beam energy from 4.65 to 6.05 GeV, and 50 different values in production angle from 0 to 1 degree, then used scipy's interpolation package to get
  //form factor values for 100 values in beam energy and 1000 values in production angle over the same range of values
  double pMagpi0, Q, primCross; //Magnitude of pi0 momentum, transferred four momentum, form factor value, primakoff cross-section
  double D, Epi0, beta, E, theta; //Constant for primakoff cross-section, pi0 energy, pi0 velocity, photon energy, production lab angle
  int iMaxStro = 0, jMaxStro = 0;  //Maximum strong cross-section bin
  double stroMax = 0;              //Maximum strong cross-section value
  int iMaxPrim = 0, jMaxPrim = 0;  //Maximum primakoff cross-section bin
  double primMax = 0;              //Maximum prinakoff cross-section bin
  int iMaxInt = 0, jMaxInt = 0;    //Maximum inerference bin
  double intMax = 0;               //Maximum interference value
  //int loops = numEvents/10;
  //string SPIZP_FFFILE = "/w/halld-scshelf2101/home/shannen/events/generators/FF/LargeTableCoul4.csv";
  ifstream FFCSV(SPIZG_FFFILE);  //Form factor table
  string line, readVal;          //Prepare to read values line by line from csv file
  //The coulomb form factors are of the form w1 + iw2, and the strong form factors are of the form w1 - w3 + (w2 - w4)i
  //The form factor tables have interpolated over each of these terms separately
  double EVal, ThetaVal, FFCoulVal, FFStroVal; //Energy value, theta value for a given form factor and those form factors 
  double w1coul, w2coul;  //w1 and w2 for coulomb
  double w1stro, w2stro, w3, w4;  //w1 - w4 for strong form factors
  int binE3, binTheta3; //Energy and theta bins
  double phiCoul, phiStro; //Complex phases of form factors
  double crossInter, stroCross;  //Interference and strong cross-sections
  double A = 208;      //Number of nucleons in Pb208
  if(FFCSV.good()){     
    int totCount = 0;   //Total number of lines read
    while(getline(FFCSV, line)){  //Loop through each line
      //cout << line << endl;  
      stringstream ss(line);  //string from a given line
      int readCount = 0;      //Number of entries in a line
      while(getline(ss, readVal, ',')){ //Parse entries by ","
	//cout << FFval2 << endl; 
	if(readCount == 0){           //First value in a line is beam energy
	  EVal = stod(readVal);
	}
	else if(readCount == 1){      //Second is production angle
	  ThetaVal = stod(readVal);
	}
	else if(readCount == 2){      //Third is the coulomb form factor
	  FFCoulVal = stod(readVal);
	}
	else if(readCount == 3){      //Fourth is the strong form factor
          FFStroVal = stod(readVal);
        }
	else if(readCount == 4){      //Fifth is w1 for the coulomb form factor
          w1coul = stod(readVal);
        }
	else if(readCount == 5){      //Sixth is w1 for the strong form factor
          w1stro = stod(readVal);
        }
	else if(readCount == 6){      //Seventh is w2 for the strong form factor
          w2coul = stod(readVal);
        }
        else if(readCount == 7){      //Eigth is w2 for the strong form factor
          w2stro = stod(readVal);
        }
	else if(readCount == 8){      //Ninth is w3
          w3 = stod(readVal);
        }
        else if(readCount == 9){      //10th is w4
          w4 = stod(readVal);
        }
	else if(readCount == 10){     //11th is the energy bin
	  //cout << "int test: " << stod(readVal);
	  binE3 = int(stod(readVal));
	}
	else{                         //12th is the theta bin
	  binTheta3 = int(stod(readVal));
	}
	readCount += 1;
	//cout << "readCount: " << readCount << endl;
      }
      readCount = 0;

      //if(FFVal2 > 1) continue;
      //Filling test histos
      FFcoul2->Fill(FFCoulVal);
      if(binE3 == 20){
	//cout << "theta" << ThetaVal << endl;
	//cout << "FFVal2 " << FFVal2 << endl;
	FFvThet2->Fill(ThetaVal, FFCoulVal);
      }

      pMagpi0 = kin2mom(EVal, ThetaVal); //From angle and energy, get the pi0 momentum ;
      Q = pMagpi0 - EVal;               //From pi0 momentum get transferred momentum Q  

      int EdistBin = beamEdistNorm->FindBin(EVal);
      
      //D = Gamma*8*alph*Z*Z*beta*beta*beta*Epi0*Epi0*Epi0*Epi0/mpi0/mpi0/mpi0;

      Epi0 = sqrt(pMagpi0*pMagpi0 + mpi0*mpi0); //Pi0 energy from pi0 momentum
      beta = pMagpi0/Epi0;                      //pi0 velocity
      D = Gamma*8*alph*Z*Z*beta*beta*beta*Epi0*Epi0*Epi0*Epi0/mpi0/mpi0/mpi0; //constant factor out front of cross-section 

      //if(ThetaVal > 0.05 && EVal < 5.95) continue;
      
      //cout << "D: " << D << endl;
      //cout << "FFVal: " << FFVal2 << endl;
      //cout << "THeta: " << ThetaVal << endl;
      //cout << "E: " << EVal << endl;
      //cout << "Q: " << Q << endl;
      //cout << "Prim cross: " << primCross << endl;
      //if(EVal < 5.5) cout << "Got here! " << endl;
      //primCross = D*FFCoulVal*sin(ThetaVal*TMath::Pi()/180)*sin(ThetaVal*TMath::Pi()/180)/Q/Q/Q/Q;
      //stroCross = A*A*FFStroVal*sin(ThetaVal*TMath::Pi()/180)*sin(ThetaVal*TMath::Pi()/180);
      primCross = D*(w1coul*w1coul + w2coul*w2coul)*sin(ThetaVal*TMath::Pi()/180)*sin(ThetaVal*TMath::Pi()/180)/Q/Q/Q/Q*sin(ThetaVal*TMath::Pi()/180)*beamEdistNorm->GetBinContent(EdistBin);  //Primakoff cross-section value at this value of energy and production angle
      stroCross = SPIZG_C*A*A*((w1stro - w3)*(w1stro - w3) + (w2stro - w4)*(w2stro - w4))*sin(ThetaVal*TMath::Pi()/180)*sin(ThetaVal*TMath::Pi()/180)*sin(ThetaVal*TMath::Pi()/180)*beamEdistNorm->GetBinContent(EdistBin); //Strong cross-section from production angle
      phiCoul = atan(w2coul/w1coul); //Coulomb complex phase
      phiStro = atan((w2stro - w4)/(w1stro - w3)); //Strong complex phase

      crossInter = 2*sqrt(primCross)*sqrt(stroCross)*(cos(phiCoul - phiStro)*cos(SPIZG_PHI) + sin(phiCoul - phiStro)*sin(SPIZG_PHI));  //Interference cross-section from phases and other cross-section
      //crossInterSin = 2*sqrt(primCross)*sqrt(stroCross)*sin(phiCoul - phiStro); 
      if(isnan(crossInter)) continue;  //If the phaases produce a NAN cross-section, terminate loop
      //if(ThetaVal > 0.1 || E < 5.8) continue;
      //cout << "EVal: " << EVal << endl;
      //cout << "ThetaVal: " << ThetaVal << endl;
      //cout << "crossInter: " << crossInter << endl;
      //cout << "binE3:" << binE3 << endl;
      //cout << "binTHeta3: " << binTheta3 << endl;

      //Fill cross-section tables
      interTot->SetBinContent(binTheta3, binE3, crossInter);
      primTot->SetBinContent(binTheta3, binE3, primCross);
      stroTot->SetBinContent(binTheta3, binE3, stroCross);
      interTotTheta->Fill(ThetaVal, crossInter);

      //Determine the maximum values and bins for each of the different production mechnaisms
      if(intMax < crossInter){
        iMaxInt = binE3;
        jMaxInt = binTheta3;

        intMax = crossInter;
      }
      if(primMax < primCross){
        iMaxPrim = binE3;
        jMaxPrim = binTheta3;

        primMax = primCross;
      }
      if(stroMax < stroCross){
        iMaxStro = binE3;
        jMaxStro = binTheta3;

        stroMax = stroCross;
      }
    }
  }
  
  //cout << "The bin with the maximum cross section is " << iMax << ", " << jMax << endl;
  //cout << "with val " << primMax << endl;
  //cout << intMax << endl;
  //cout << stroMax << endl;

  //The below segment of code takes in whatever production mechamisms you want turned on, and finds which of them has the maximum value (only works for
  //Pb208 since the primakoff term is so dominant and the strong is much larger than the interference), and writes them to one cross-section table for 
  double toUMax = 0;
  double iMaxToU = 0;
  double jMaxToU = 0;
  if(SPIZG_PRIM){
    toUTot->Add(primTot);
    toUMax = primMax;
    iMaxToU = iMaxPrim;
    jMaxToU = jMaxPrim;
  }
  if(SPIZG_STRO){
    toUTot->Add(stroTot);
    if(!SPIZG_PRIM){
    toUMax = stroMax;
      iMaxToU = iMaxStro;
      jMaxToU = jMaxStro;
    }
    else if(stroMax > primMax){
      cout << "Got here!";
      toUMax = stroMax;
      iMaxToU = iMaxStro;
      jMaxToU = jMaxStro;
    }
  }
  if(SPIZG_INTER){
    toUTot->Add(interTot);
    if(!SPIZG_STRO && !SPIZG_PRIM){
      toUMax = intMax;
      iMaxToU =	iMaxInt;
      jMaxToU =	jMaxInt;
    }
    else if(!SPIZG_PRIM && intMax > stroMax){
      toUMax = intMax;
      iMaxToU = iMaxInt;
      jMaxToU = jMaxInt;
    }
    else if(!SPIZG_STRO && intMax > primMax){
      toUMax = intMax;
      iMaxToU = iMaxInt;
      jMaxToU = jMaxInt;
    }
    else if(intMax > stroMax && intMax > primMax){
      cout << "Got here too!";
      toUMax = intMax;
      iMaxToU = iMaxInt;
      jMaxToU = jMaxInt;
    }
    
  }

  //int maxBinTot = toUTot->GetMaximum();
  //toUMax = toUTot->GetBinContent(maxBinTot);
  //cout << "Max cross-section value was: " << toUMax << endl;
  //return;

  //And at long last we've reached accept reject
  double accVal, accValE, dNdEdTheta;
  int binE, binE2, binTheta;
  double p1, p2, p3, ET;
  int evtsSoFar = 0;
  for(int k = 0; k < MAXEVENTS; k++){ //For loop is capped at 1000x the desired number of events    
    accVal = rand_num()*toUMax;               //Accept value for primakoff cross-section 

    randE2 = rand_num()*(maxE - minE) + minE; //Random photon energy value

    //binE = beamEdist->FindBin(randE2);                                                                                                                                                                       
    //dNdE = beamEdist->GetBinContent(binE);  //Find the contents of the bin of the beam energy distribution at the random value of E

    //accValE = rand_num(k)*dNdEmax;          //Random number for accept/reject

    //if(dNdE < accValE) continue;           //Reject if the value of the beam energy distribution is below a random number between 1 and 0
                                           //times the max value of beam energy  
    randTheta = 2*rand_num();             //Produce random number in theta 
    //randTheta = acos(1 - 4*rand_num(k));
    
    //cout << "Theta" << randTheta << endl;
    //cout << "Energy: " << randE2 << endl;
    
    binE2 = toUTot->GetYaxis()->FindBin(randE2);       //Get energy bin from cross-section table  
    binTheta = toUTot->GetXaxis()->FindBin(randTheta); //Get theta bin from cross-section table 

    //cout << "BinE2: " << binE2 << endl;
    //cout << "binTheta: " << binTheta << endl;
    
    dNdEdTheta = toUTot->GetBinContent(binTheta, binE2);//*sin(randTheta*TMath::Pi()/180); //Find cross-section at that bin

    //cout << "max" << toUMax << endl;
    //cout << "dNdETheta" << endl;
    //cout << "accVal" << endl;
    
    if(dNdEdTheta < accVal) continue; //Perform accept reject with that  
    
    Eg->Fill(randE2);                  //Fill pi0 and beam energy distributions 
    thetapi0->Fill(randTheta);

    //randPhi = asin(2*TMath::Pi()*(1 - rand_num())/P);
    randPhi = Gmin(rand_num(), SPIZG_POL_FRAC, SPIZG_POL_DEG); //Produce values of phi with distribution of NPP equation 
    phiPi0->Fill(randPhi);

    //Epi0 = sqrt(pMag*pMag + mpi0*mpi0);
    double pMag = kin2mom(randE2, randTheta);  //Find pi0 momenta from kinematics
    Epi0 = sqrt(pMag*pMag + mpi0*mpi0); 
    ET = randE2 + Mt - Epi0; 

    //Find components of 3-momenta
    p1 = pMag*sin(randTheta*TMath::Pi()/180)*cos(randPhi*TMath::Pi()/180);                                                                                                                                                             
    p2 = pMag*sin(randTheta*TMath::Pi()/180)*sin(randPhi*TMath::Pi()/180);                                                                                                                                                             
    p3 = pMag*cos(randTheta*TMath::Pi()/180); 

    //Fill csv with four momenta
    myfile << Epi0 << "," <<  p1 << "," << p2 << "," << p3 << ",";
    myfile << ET << "," << -p1 << "," << -p2 << "," << randE2 - p3;                                                                                                                                                                    
    myfile << "\n";

    evtsSoFar += 1;
    if(evtsSoFar >= SPIZG_NUMEVENTS) break; //Stop once we've generated enough events
  }

  // Close CSV file
  myfile.close();

  // Write histograms to ROOT file and close
  if (fout && !fout->IsZombie()) {
    fout->cd();
    //randNum_test->Write();
    Eg->Write();
    phiPi0->Write();
    //randX_hist->Write();
    thetapi0->Write();
    //FFcoul->Write();
    //primHist->Sumw2(kFALSE);
    //primHist->Write();
    toUTot->Write();
    interTotTheta->Sumw2(kFALSE);
    interTotTheta->Write();
    stroTot->Write();
    primTot->Write();
    interTot->Write();
    FFcoul->Write();
    FFcoul2->Write();
    FFvThet->Write();
    FFvThet2->Write();
    
    fout->Close();
    delete fout;
  }
}

int main(int argc, char* argv[]){
  // Command-line arguments: <jobnum> <nevents>
  // Similar to RBHG Fortran: seed, jobnum, nevents
  // For SPIZG we use: jobnum (acts as seed base), nevents
  
  if(argc != 3){
    cerr << "Usage: " << argv[0] << " <jobnum> <nevents>" << endl;
    cerr << "  jobnum:  Job number (used for RNG seed)" << endl;
    cerr << "  nevents: Number of events to generate" << endl;
    return 1;
  }
  
  int jobnum = atoi(argv[1]);
  int nevents = atoi(argv[2]);
  
  if(jobnum < 0 || nevents <= 0){
    cerr << "Error: jobnum and nevents must be positive integers" << endl;
    return 1;
  }
  
  cout << "===============================================" << endl;
  cout << "SPIZG Generator (pi0 Primakoff Production)" << endl;
  cout << "===============================================" << endl;
  cout << "Job Number: " << jobnum << endl;
  cout << "Events:     " << nevents << endl;
  
  // Initialize GSL random number generator ONCE
  // Seed = jobnum + current time for uniqueness across jobs
  gsl_rng_env_setup();
  const gsl_rng_type * T = gsl_rng_default;
  global_rng = gsl_rng_alloc(T);
  
  struct timeval CurrentTime;
  gettimeofday(&CurrentTime, 0);
  unsigned long seed = (unsigned long)jobnum * 1000000 + CurrentTime.tv_usec;
  gsl_rng_set(global_rng, seed);
  
  cout << "RNG Seed:   " << seed << endl;
  cout << "===============================================" << endl;
  
  // Run the generator
  genNumTot(global_rng, nevents);
  
  // Clean up
  gsl_rng_free(global_rng);
  
  cout << "Generation complete!" << endl;
  return 0;
}
