/*!
 * \class RadiationProcess
 *
 * \brief This class describe a generic radiatioon process
 * 
 * Synchrotron Radiation and Inverse Comton scattering inherit from here
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
#ifndef RadiationProcess_H
#define RadiationProcess_H 1

#include "SpectObj.h"
#include"GRBConstants.h"

class RadiationProcess
{

 public:
  RadiationProcess(){;}
  RadiationProcess(SpectObj);
  virtual ~RadiationProcess(){;}
  
  virtual double processFlux(double, double, double);
  virtual double electronNumber(double, double, double, double, double,
				double, double);
  
  double timeShiftForDispersion(const double time, 
				const double E, 
				const double distance_to_source);

  double comovingTime(const double time, 
		      const double gamma, 
		      const double E, const double D);

  //! return the Spectrum object.
  inline SpectObj getSpectrumObj() {return m_spectrumObj;}

  //! Helper function. /param B is in Gauss, Umag is in erg.
  inline virtual double Umag(const double B)
    {
      return pow(B/1.0e+4,2.)/(200.0*M_PI*cst::mu0);
    }
 
 protected:
  SpectObj m_spectrumObj;

};

#endif;
