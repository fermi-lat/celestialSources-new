/*!
* \class GRB
*
* \brief This class serves as a base class for GbmGrb and LatGrb classes.
* \author Sandhia Bansal
*
*
* This class implements methods to create GRBs.
*
* If the user specified burst parameters (duration, number of pulses, flux, 
* fraction, broken power law indices, and peak energy), it creates the 
* specific burst and writes it to a file in the current working directory.
*
* If the user specified no input or only a directory name, it creates n GRBs 
* where n is specified in the GRBobsConstants.h file.  If a directory name is 
* specified, the GRBs are written to files in that directory; otherwise, the 
* they are written to files in the current working directory.
*
* The user can also specify names of directory and file containing GRB data in 
* which case, it simply reads the data from the file
*/


// File: GRB.h
//
// GRB simulation

// Interface:
// This class is instantiated in three ways:
// -  No input
//       In this mode, it generates nbsim number of bursts and for each burst, it creates and records a photon list (time,energy)
//		 The number nbsim is generated by the GRBobsConstants class.
// -  Input: Filename
//		 In this mode, it will read the photon list (time,energy) generated by the first option
// -  Input: duration, flux, fraction, power law index, npulse, flag
//		 In this mode, it creates a photon list for the burst specified by the input parameters and if the flag is set,
//			records it in a file
//
#ifndef GRB_H
#define GRB_H

#include "GRBobsUtilities.h"
#include <vector>
#include <string>
#include <functional>  // for unary_function and plus
#include <fstream>


class HepRandomEngine;
class TimeEnergy;
class GlobalData;


class GRB
{
    friend class GRBobsSpectrum;  // to give it access to the private members
    friend std::ifstream &operator >> (std::ifstream &is, GRB &grb);
    
public:
    // Constructors
    GRB();
    GRB(const std::vector<std::string> &paramVector);
    GRB(HepRandomEngine *engine, const std::string &prefix, const double duration, const int npuls, const double flux,
        const double fraction, const double alpha, const double beta, const double epeak, const double specnorm,
        const bool flag);
    
    
    ~GRB();
    GRB(const GRB &right);   
    GRB &operator=(const GRB &right);
    
    
    // Accessor functions
    std::pair<float,float> dir() const                 { return m_grbdir; }
    double univFWHM() const                            { return m_univFWHM; }
    const std::vector<double> &specnorm() const        { return m_specnorm; }
    long nphoton() const                               { return m_nphoton; } 
    const std::vector<TimeEnergy> &photonlist() const  { return m_photonlist; }
    GlobalData *globalData() const                     { return m_globalData; }
    
    void setSpecnorm(const std::vector<double> &specnorm)    { m_specnorm = specnorm; }
    
    
    
    // Class methods
    // Create "n" GRBs
    createGRB(HepRandomEngine *engine, const std::string &prefix, const std::string &dir=0);
    
    // Create GRB for specified input parameters
    createGRB(HepRandomEngine *engine, const std::string &prefix, const double duration, const int npuls, const double flux, 
        const double fraction, const double alpha, const double beta, const double epeak, const double specnorm, 
        const bool flag);
    
protected:
    // data members
    std::pair<float,float>   m_grbdir;
    double                   m_univFWHM;
    std::vector<double>      m_specnorm;
    long                     m_nphoton;
    std::vector<TimeEnergy>  m_photonlist;
    GlobalData               *m_globalData;
    
    
    
    // private accessor functions
    std::vector<TimeEnergy> &photonlist()  { return m_photonlist; }
    
    virtual long calcNphoton(HepRandomEngine *engine)   { return 0; }
    void makeTimes(HepRandomEngine *engine, const double ethres);
    
private:
    // Class methods
    // Create base string used to generate output file name
    std::string baseFilename(const std::string &prefix, const std::string &dir=0) const;
    
    // Use base name generated by baseFilename to create the output file name
    std::string outputFilename(const std::string &base, const long isim) const;
    
    // Return a copy of itself
    GRB *clone() const;
    
    // Compute the direction of the burst
    std::pair<float,float> direction(HepRandomEngine *engine) const;
    
    // Use pulse data to compute the time profiles
    void getTimes(HepRandomEngine *engine, const double ethres, const long nphots, const long deltbinsleft, const long iphotoff, 
        const double tmax, const std::vector<double> &pulse);
    
    // Makes BATSE-like GRB time profiles, placing GLAST photons a la cumulative BATSE intensity, but in narrower pulses:
    virtual void makeGRB(HepRandomEngine *engine)   {}
    
    // Read GRB from input file
    readGRB(const std::vector<std::string> &paramVector);
    
    void swap(GRB &other) throw();
};



class GlobalData
{
public:
    GlobalData()  {};
    
    inline void setDuration(const double duration)            { m_duration=duration; }
    inline void setFlux(const double flux)                    { m_flux=flux; }
    inline void setFraction(const double fraction)            { m_fraction=fraction; }
    inline void setAlpha(const double alpha)				   { m_alpha=alpha; }
    inline void setBeta(const double beta)					   { m_beta=beta; }
    inline void setEpeak(const double epeak)                  { m_epeak=epeak; }
    inline void setSpecnorm(const double specnorm)            { m_specnorm=specnorm; }
    inline void setNpuls(const int npuls)                     { m_npuls=npuls; }
    
    inline double duration()       { return m_duration; }
    inline double flux()           { return m_flux; }
    inline double fraction()       { return m_fraction; }
    inline double alpha()		    { return m_alpha; }
    inline double beta()		    { return m_beta; }
    inline double epeak()		    { return m_epeak; }
    inline double specnorm()       { return m_specnorm; }
    inline int npuls()             { return m_npuls; }
    
    GlobalData *clone() const  { return new GlobalData(*this); }
    
private:
    double m_duration;
    double m_flux;
    double m_fraction;
    double m_alpha;
    double m_beta;
    double m_epeak;
    double m_specnorm;
    int    m_npuls;
};



class TimeEnergy
{
public:
    TimeEnergy() { };
    
    inline void setTime(double time)   {m_time=time;}
    inline void setEnergy(double energy) {m_energy=energy;}
    
    inline double time()   {return m_time;}
    inline double energy() {return m_energy;}
    
private:
    double m_time;
    double m_energy;
};



class timeCmp
{
public:
    bool operator()(TimeEnergy &data1, TimeEnergy &data2)
    {
        return data1.time() < data2.time();    
    }
};


// Output operator
std::ofstream &operator<<(std::ofstream &os, const GRB &grb);

#endif // GRB_MAKER_H
