// FILE: GRBobsSpectrum.cxx

#include <memory>  // for auto_ptr
#include <fstream>


#include "GRBmaker.h"
#include "GRB.h"
#include "GRBobsSpectrum.h"
#include "GRBobsConstants.h"
#include "src/GPS.h"
#include "src/SpectrumFactory.h"
#include "CLHEP/Random/RandFlat.h"



static std::ofstream os("gleam.lis");


// Constructors

GRBobsSpectrum::GRBobsSpectrum(const std::string &params)
    :ISpectrum(),
     m_title("GRBobsSpectrum"),
     m_particleName("gamma")
{
    std::vector<std::string> paramVector;
    parseParamList(params, paramVector);
    
    std::auto_ptr<GRBmaker> grbMaker(new GRBmaker);
    std::auto_ptr<GRB> p(grbMaker->create(paramVector));
    m_grb = p.release();
}


GRBobsSpectrum::GRBobsSpectrum(const double duration, const int npuls, const double flux, const double fraction, 
                               const double alpha, const double beta, const double epeak, const double specnorm, const bool flag)
                               :ISpectrum(),
                               m_title("GRBobsSpectrum"),
                               m_particleName("gamma")
{
    std::auto_ptr<GRBmaker> grbMaker(new GRBmaker);
    std::auto_ptr<GRB> p(grbMaker->create(duration, npuls, flux, fraction, alpha, beta, epeak, specnorm, flag));
    m_grb = p.release();
}


// Copy Constructor
GRBobsSpectrum::GRBobsSpectrum(const GRBobsSpectrum &right)
    :ISpectrum(),
     m_title(right.m_title),
     m_particleName(right.m_particleName)
{
    if (right.m_grb)
        m_grb = right.m_grb->clone();
}


// Destructor
GRBobsSpectrum::~GRBobsSpectrum()
{
    delete m_grb;
    m_grb = 0;
}


// Helper function for the assignment operator
void GRBobsSpectrum::swap(GRBobsSpectrum &other) throw()
{
    std::swap(m_title, other.m_title);
    std::swap(m_particleName, other.m_particleName);
    std::swap(m_grb, other.m_grb);
}


// Assignment Operator
GRBobsSpectrum &GRBobsSpectrum::operator=(const GRBobsSpectrum &right)
{
    GRBobsSpectrum temp(right);   // does all the work
    swap(temp);   // this can't throw
    return *this;
}


// Parse input parameter list obtained from the xml file
void GRBobsSpectrum::parseParamList(const std::string &input, std::vector<std::string>& output) const
{ 
    std::string::size_type i = input.find_last_of("/");
    
    if (i > 0)
    {
        if (i == input.size()-1)
            output.push_back(input);
        else
        {
            std::string temp = input.substr(0,i+1);
            output.push_back(temp);
            temp = input.substr(i+1);
            output.push_back(temp);
        }
    }
    
    else
    {
        std::string temp = input;
        output.push_back(temp);
    }
}


double GRBobsSpectrum::flux(double time) const
{ 
    return m_grb->globalData()->flux();
}


float GRBobsSpectrum::fraction(float energy) 
{
    return m_grb->globalData()->fraction();
}


double GRBobsSpectrum::interval(double time)
{
    static double currentPTime = 0.0;
    double nextPTime = 0.0;
    double intv = 0.0;
    
    
    nextPTime = nextTime();
    intv = nextPTime - currentPTime;
    currentPTime = nextPTime;
    
    return intv;
}


//JCT needs const to match pure virtual method
float GRBobsSpectrum::operator () (float randomNumber) const  
{
    return (float) nextEnergy();
}


// returns next available energy
// expects to be called once per time - so iterator should never become invalid before times run out
double GRBobsSpectrum::nextEnergy() const
{
    static std::vector<TimeEnergy>::iterator it = m_grb->photonlist().begin();
    
    if (it != m_grb->photonlist().end())
        return (*it++).energy();
    else
        return -1;
}


// returns next available energy to the simulation
double GRBobsSpectrum::energySrc(HepRandomEngine *engine, double time)
{
    return nextEnergy();
}


// returns next available time
double GRBobsSpectrum::nextTime() const
{
    static std::vector<TimeEnergy>::iterator it = m_grb->photonlist().begin();
    
    if (it != m_grb->photonlist().end())
        return (*it++).time();
    else
        return 1.0e+6;
}


std::pair<float,float> GRBobsSpectrum::dir(float energy) const
{
    return m_grb->dir();
}


std::pair<double,double> GRBobsSpectrum::dir(double energy, HepRandomEngine *engine)
{
    return dir(energy);
}