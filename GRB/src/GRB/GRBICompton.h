/*!
 * \class GRBICompton
 * \brief ICompton emission process.
 * 
 * 
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
#include <vector>
#include"GRBConstants.h"
#include "SpectObj.h"
#include "RadiationProcess.h"
#include "GRBShock.h"

#ifndef GRBICOMPTON_H
#define GRBICOMPTON_H 1

class GRBICompton : virtual public RadiationProcess
{
 public:
  GRBICompton();
  GRBICompton(SpectObj spectrumObj);
  ~GRBICompton(){;}

  /*! Given a Shock compute the spectrum at a certain /param time.
   *  /param angle is needed to compute the observed spectrum.
   */
  void load(const GRBShock *Shock,
	    const double time  = 0.0,
	    const double angle = 0.0,
	    const double distance_to_source =0.0);

  /*! General Interface for computing the IC spectrum at a certain /param time.
   */
  void load(const double time      = 0.0,
	    const double angle     = 0.0,
	    const double distance_to_source =0.0,
	    const double GAMMAF    = 500.,
	    const double B         = 1.0e+6,
	    const double N0        = 1.0e+50,
	    const double gamma_min = 1.,
	    const double gamma_max = 1.0e+10,
	    const double Vol       = 1.0e+40,
	    const double dr        = 1.0e+10); 

 private:
};
#endif


