#ifndef __BODY_H
#define __BODY_H

#include <string>
#include <math.h>
#include <cmath>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "tools.h"
#include "physconst.h"
#include "global.h"

//#define EarthAtm_DEBUG

using namespace std;

class Body{
    public:
        string name;
        
        Body();
        
        class Track{
            public:
                double x;
                double xini;
                double xend;
		double L,LE;
                Track(double,double);
                Track();
        };
        virtual double density(Track*);
        virtual double ye(Track*);
};

// type defining
typedef Body::Track GenericTrack;

class Vacuum: public Body {
    public:        
        Vacuum();
        
        class Track: public Body::Track {
            public :      
                Track(double,double);
                Track();
        };
        
        double density(GenericTrack*);
        double ye(GenericTrack*);
};

class ConstantDensity: public Body{
    public:
        double constant_density;
        double constant_ye;
        
        ConstantDensity(double,double);
        
        class Track: public Body::Track{
            public :              
                Track(double,double);
                Track();            
        };
        
        double density(GenericTrack*);
        double ye(GenericTrack*);
};

class VariableDensity: public Body{
    public:
        double rhomin;
        double rhomax;
        double constant_ye;
        
        VariableDensity(double,double,double);    
        
        class Track: public Body::Track{
            public :          
                Track(double,double);
                Track();            
        };
        
        double density(GenericTrack*);
        double ye(GenericTrack*);        
};

class Star: public Body{
    public:
        double radius;
        double initialradius;
        double initialdensity;
        double constant_ye;
        
        Star(double,double,double);
        
        double rdensity(double);
        
        class Track: public Body::Track{
            public :              
                Track(double,double);
                Track();            
        };
        
        double density(GenericTrack*);
        double ye(GenericTrack*);  
};

class Earth: public Body{
    public:
        double radius;
        double constant_ye;        
        Earth();
        
        double rdensity(double);
        
        class Track: public Body::Track{
            public :
                double baseline;
                Track(double,double,double);
                Track();            
        };
        
        double density(Body::Track*);
        double ye(GenericTrack*);          
};

class AGN: public Body{
    public:
        double radius;
        double initialdensity;
        double constant_ye;
        
        AGN(double,double,double);
        
        double rdensity(double);
        
        class Track: public Body::Track{
            public :            
                Track(double,double);
                Track();            
        };
        
        double density(Body::Track*);
        double ye(Body::Track*);  
};

class Sun: public Body{
    public:
        Table sun_model;
        Table sun_model_nele;
        double* sun_radius;
        double* sun_density;
        double* sun_xh;
        double* sun_nele_radius;
        double* sun_nele;
        
        int arraysize,arraysize_2;

        double radius;
        
        gsl_spline * inter_density;
        gsl_interp_accel * inter_density_accel;
        
        gsl_spline * inter_rxh;
        gsl_interp_accel * inter_rxh_accel;
        
        gsl_spline * inter_nele;
        gsl_interp_accel * inter_nele_accel;
        
        Sun();
        
        double rdensity(double);
        double rxh(double);
        
        class Track: public Body::Track{
            public :            
                Track(double,double);
                Track();            
        };
        
        double density(Body::Track*);
        double ye(Body::Track*);

};

class SunASnu: public Body{
    public:
        Table sun_model;
        double* sun_radius;
        double* sun_density;
        double* sun_xh;
        
        int arraysize;

        double radius;
        
        gsl_spline * inter_density;
        gsl_interp_accel * inter_density_accel;
        
        gsl_spline * inter_rxh;
        gsl_interp_accel * inter_rxh_accel;
        
        SunASnu();
        
        double rdensity(double);
        double rxh(double);
        
        class Track: public Body::Track{
            public :
                double b_impact;
                double radius_nu;
                Track(double,double);
                Track();            
        };
        
        double density(Body::Track*);
        double ye(Body::Track*);
};

class EarthAtm: public Body{
    public:
        double radius;
        
        EarthAtm();
        
        double rdensity(double);
        double rxh(double);
        
        class Track: public Body::Track{
            public :
                double phi;
                double cosphi;
                double radius_nu;
                double atmheight;
                
                Track(double);
                Track();            
        };
        
        double density(Body::Track*);
        double ye(Body::Track*);
};



// type defining
typedef Body::Track Track;

#endif
