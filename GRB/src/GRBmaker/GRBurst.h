/*!
* \class GRBurst
*
* \brief This class serves as a base class for GbmGrb and LatGrb classes.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
*
* This class implements methods to create GRBs.
*
* If the user specified burst parameters (duration, number of pulses, flux, 
* fraction, broken power law indices, and peak energy), it creates the 
* burst for specified parameters.
*
* If the class was instantiated with no parameters or only a directory name, 
* it creates n GRBs where n is specified in the GRBobsConstants.h file.  
*
* If a directory name is specified, the GRBs are written to files in that 
* directory; otherwise, they are written to files in the current working directory.
*
* The user can also specify names of directory and file containing burst data in 
* which case, it simply reads the data from the file to make it available for the
* simulation.
*/


#ifndef GRB_BURST_H
#define GRB_BURST_H

#include "GRBobsUtilities.h"
#include <vector>
#include <string>
#include <functional>  // for unary_function and plus
#include <fstream>

class HepRandomEngine;
class TimeEnergy;
class GlobalData;



class GRBurst
{
    friend class GRBobsSpectrum;  // to give it access to the private members
    friend std::ifstream &operator >> (std::ifstream &is, GRBurst &grb);
    
public:
    /*!
     * \brief Constructor.
     */
    GRBurst();

    /*!
     * \brief Constructor.
     */
    GRBurst(const std::vector<std::string> &paramVector);

    /*!
     * \brief Constructor.
     */
    GRBurst(HepRandomEngine *engine, const double duration, const int npuls, 
        const double flux, const double fraction, const double alpha, 
        const double beta, const double epeak, const double specnorm, 
        const bool flag);
      
    /*!
     * \brief Destructor.
     */
    virtual ~GRBurst();
      
    /*!
     * \brief Copy Constructor.
     */
    GRBurst(const GRBurst &right);   
      
    /*!
     * \brief Assignment Operator.
     */
    GRBurst &operator=(const GRBurst &right);
    
    

    // Accessors
    std::pair<float,float> dir() const                 { return m_grbdir; }
    double univFWHM() const                            { return m_univFWHM; }
    const std::vector<double> &specnorm() const        { return m_specnorm; }
    long nphoton() const                               { return m_nphoton; } 
    const std::vector<TimeEnergy> &photonlist() const  { return m_photonlist; }
    GlobalData *globalData() const                     { return m_globalData; }
    
    void setSpecnorm(const std::vector<double> &specnorm)   { m_specnorm = specnorm; }
    
    
    
    /*!
     * \brief Create "n" GRBs where "n" is specified in GRBobsConstants.h file.
     *
     * \param prefix    "GBM" or "LAT".
     * \param dir       Name of directory that the output files are to be written to.
     */
    void createGRB(HepRandomEngine *engine, const std::string &prefix, 
        const std::string &dir=0);
    
    /*!
     * \brief Create burst for specified input parameters.
     */
    void createGRB(HepRandomEngine *engine, const double duration, 
        const int npuls, const double flux, const double fraction, 
        const double alpha, const double beta, const double epeak, 
        const double specnorm, const bool flag);
    
protected:
    // data members
    std::pair<float,float>   m_grbdir;             //! Direction of the burst
    double                   m_univFWHM;           //! Universal width for the pulses within the current burst
    std::vector<double>      m_specnorm;           //! Spectral normalization
    long                     m_nphoton;            //! Number of photons generated by a burst
    std::vector<TimeEnergy>  m_photonlist;         //! Burst data ((time,energy) pair for each photon
    GlobalData               *m_globalData;        //! Pointer to burst specific data (duration, flux etc.)
    
    
    
    // Accessors
    std::vector<TimeEnergy> &photonlist()  { return m_photonlist; }
    
    /*!
     * \brief Calculate number of photons to be generated by the burst.
     */
    virtual long calcNphoton(HepRandomEngine *engine)   { return m_nphoton; }

    /*!
     * \brief  Make BATSE-like burst time profiles, placing GLAST photons a la cumulative BATSE intensity, but in narrower pulses.
     *
     *  (1) m_npuls = number of pulses, proportional to BATSE duration.
     *  (2) pulse peak amplitude is random (0.0=>1.0); sort amps in descending amp order.
     *  (3)	scramble amps of {1st,2nd} halves of pulses, separately (leaves profile assymmetric).
     *  (4)	center of pulse time is random within duration.  sort the times in ascending order.
     *  (5)	pulse width is drawn from BATSE width distribution for bright bursts (attributes paper), 
     *		scaled to GLAST energies using width ~E^-0.333.
     *  (6)	make m_npuls pulses with "bisigma" shapes => sum to produce time profile.
     *  (7)	form cumulative distribution of BATSE-like intensity.
     *  (8) distribute the m_nphoton photons according to cumulative intensity => GRBtimes.
     *  (9)	offset the photon times according to
     *           (a) energy dependence, width ~E^-0.333 and
     *           (b) time of peak, also proportional to E^-0.333.
     */
    void makeTimes(HepRandomEngine *engine, const double ethres);
    
private:
    /*!
     * \brief Create base string used to generate output file name
     */
    std::string baseFilename(const std::string &prefix, const std::string &dir=0) const;
    
    /*! 
     *\brief Use base name generated by baseFilename to create the name of the output file.
     */
    std::string outputFilename(const std::string &base, const long isim) const;
    
    /*!
     * Return a copy of itself
     */
    GRBurst *clone() const;
    
    /*!
     * \brief Compute the direction of the burst
     */
    std::pair<float,float> direction(HepRandomEngine *engine) const;
    
    /*!
     * /brief Use pulse data to compute the time profiles.
     */
    void getTimes(HepRandomEngine *engine, const double ethres, 
        const long nphots, const long deltbinsleft, const long iphotoff, 
        const double tmax, const std::vector<double> &pulse);
    
    /*!
     * \brief Make BATSE-like burst time profiles, placing GLAST photons a la cumulative BATSE intensity, but in narrower pulses.
     */
    virtual void makeGRB(HepRandomEngine *engine)   {}
    
    /*!
     * \brief Read burst data (photon list) generated by a previous run.
     */
    void readGRB(const std::vector<std::string> &paramVector);
    
    /*!
     * \brief Used by the assignment operator.
     */
    void swap(GRBurst &other) throw();
};



/*!
 * \brief This class provides the interface for burst specific data.
 */
class GlobalData
{
public:
    GlobalData()  {};
    
    inline void setDuration(const double duration)      { m_duration=duration; }
    inline void setFlux(const double flux)              { m_flux=flux; }
    inline void setFraction(const double fraction)      { m_fraction=fraction; }
    inline void setAlpha(const double alpha)   	        { m_alpha=alpha; }
    inline void setBeta(const double beta)			    { m_beta=beta; }
    inline void setEpeak(const double epeak)            { m_epeak=epeak; }
    inline void setSpecnorm(const double specnorm)      { m_specnorm=specnorm; }
    inline void setNpuls(const int npuls)               { m_npuls=npuls; }
    
    inline double duration()       { return m_duration; }
    inline double flux()           { return m_flux; }
    inline double fraction()       { return m_fraction; }
    inline double alpha()		   { return m_alpha; }
    inline double beta()		   { return m_beta; }
    inline double epeak()		   { return m_epeak; }
    inline double specnorm()       { return m_specnorm; }
    inline int npuls()             { return m_npuls; }
    
    GlobalData *clone() const  { return new GlobalData(*this); }
    
private:
    double m_duration;    //! Duration of the burst
    double m_flux;        //! Flux
    double m_fraction;    //! Ratio of peak flux to maximum for all generated bursts
    double m_alpha;       //! Power Law Index
    double m_beta;        //! Power Law Index
    double m_epeak;       //! Peak energy
    double m_specnorm;    //! Spectral Normalization
    int    m_npuls;       //! Number of pulses in the burst
};



/*!
 * \brief This class provides interface for the burst data.
 */
class TimeEnergy
{
public:
    TimeEnergy() { };
    
    inline void setTime(double time)   {m_time=time;}
    inline void setEnergy(double energy) {m_energy=energy;}
    
    inline double time() const   {return m_time;}
    inline double energy() const {return m_energy;}
    
private:
    double m_time;
    double m_energy;
};



/*!
 * \brief This class provides the functionality to sort data in time order.
 */
class timeCmp
{
public:
    bool operator()(const TimeEnergy &data1, const TimeEnergy &data2)
    {
        return data1.time() < data2.time();    
    }
};


//! Output operator - needed to write burst data to a file.
std::ofstream &operator<<(std::ofstream &os, const GRBurst &grb);


#endif // GRB_MAKER_H
