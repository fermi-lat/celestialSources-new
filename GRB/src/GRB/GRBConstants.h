/*! 
 * \class GRBConstants
 * \brief Class instantiated to access general parameters and constants.
 *  
 *  This should be a temporary solution: we need to find a better way to deal
 *  with it!
 * \author Nicola Omodei nicola.omodei@pi.infn.it
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 */

#ifndef GRBCONSTANTS_HH
#define GRBCONSTANTS_HH 1

namespace cst
{
  /// Universal constants :
  const double mpc2      = 938.2;
  const double erg2MeV   = 624151.0;
  const double mpc2cm    = 3.0857e+24;

  const double G2MeV     = 815.78;
  const double st        = 6.65225e-25;
  const double mec2      = 0.510999;
  //! light speed in cm/sec.
  const double c         = 2.98e+10;
  const double c2        = c*c;
  //! Planck constant in eV*sec.
  const double hplanck   = 4.13567e-15;
  const double pi        = 3.1415926535897932385;
  const double Hubble    = 6.5e+1; 
  const double wzel      = 0.0;

  const double csi       = 1.0;
  const double alphae    = .33;
  const double alphab    = .33;
  const double p         = 2.5;
  const double viscosity = 0.;
  /// Internal Parameters
  const double enmax     = 1.0e+12;
  //! Min. photon energy detectable by GLAST: set to 1 MeV
  //  const double enph      = 1.0e+6; 
  const double enph      = 5.0e+4; //BATSE band
  //! Minimal Temporal separation between 2 photons
  const double DeadTime  = 1.0e-5; //sec
  const double enmin     = 1.0e+3;
  const double dt1       = 10000.;
  const int nstep        = 200;
  const int enstep       = 100;
  //! flag =[0,1], if ==0, No inverse compton;
  const float flagIC     = 1;
  //RandomFlag==1 -> the parameters are choosing random 
}

class GRBConstants 
{ 
  
 public:
  //! Constructor calls method readParam wich gets an external file
  GRBConstants();

  ~GRBConstants() { }
  
  //! Parameters are read from a file using facilities::Util::expandEnvVar ethod
  void ReadParam();
  void Print();
  void Save();
  double SelectRandom(double min=0.0, double max=1.0);
  //! Number of shells generated from the source
  inline int Nshell() {return nshell;}
  //  inline void setNshell(int value=10){nshell=value;}
  void setNshell(int value=10);
  
  //! redshift of the source
  inline double Redshift() {return redshift;}
  //  inline void setRedshift(double value=1.0){redshift=value;}
  void setRedshift(double value=1.0);
  
  //! Total Energy available (ergs)
  inline double Etot(){return etot;}
  //  inline void setEtot(double value){etot=value;}
  void setEtot(double value);
  
  //! Initial separation between shells (cm)
  inline double R0(){return r0;}
  //  inline void setR0(double value){r0=value;}
  void setR0(double value);
  
  //! Initial thickness of the shells (cm)
  inline double T0(){return t0;}
  inline void setT0(double value);
  
  //! Minimum Lorentz factor of the shells
  inline double Gamma0(){return g0;}
  //inline void setGamma0(double value=100.0){g0=value;}
  void setGamma0(double value=100.0);
  
  //! Maximum Lorentz factor of the shells
  inline double DGamma(){return g1;}
  //inline void setDGamma(double value=100.0){g1=value;}
  void setDGamma(double value=100.0);
 
  
  int nshell;
  double redshift;
  double etot;
  double r0;
  double t0;
  double g0,g1;

};

#endif


