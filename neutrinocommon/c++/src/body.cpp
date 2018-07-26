#include "body.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )

PhysConst param;

/*
-----------------------------------------------------------------------
         BODY CLASS DEFINITIONS
-----------------------------------------------------------------------         
*/


Body::Body()
        {
            name = "GenericBody";
        }        
        
// track constructor        
Body::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Body::Track::Track()
        {
            //pass
        }
        
double Body::density(Track *track_input){
            return 0.0;
        }
        
double Body::ye(Track *track_input){
            return 1.0;
        }

/*
----------------------------------------------------------------------
         VACUUM CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// vacuum constructor
Vacuum::Vacuum()
        {
            name = "Vacuum";
        }        
        
// track constructor        
Vacuum::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Vacuum::Track::Track()
        {
            //pass
        }        

double Vacuum::density(GenericTrack *track_input){
            return 0.0;
        }
        
double Vacuum::ye(GenericTrack *track_input){
            return 1.0;
        }
        
/*
----------------------------------------------------------------------
         ConstantDensity CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
ConstantDensity::ConstantDensity(double density_input,double ye_input)
        {
            name = "ConstantDensity";
            constant_density = density_input;
            constant_ye = ye_input;
        }
// track constructor        
ConstantDensity::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
ConstantDensity::Track::Track()
        {
            //pass
        }
  
double ConstantDensity::density(GenericTrack *track_input)
        {
            return constant_density;
        }
        
double ConstantDensity::ye(GenericTrack *track_input)
        {
            return constant_ye;
        }        

/*
----------------------------------------------------------------------
         VariableDensity CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
VariableDensity::VariableDensity(double rhomin_input,double rhomax_input,double ye_input)
        {
            name = "VariableDensity";
            rhomin = rhomin_input;
            rhomax = rhomax_input;
            constant_ye = ye_input;
        }
// track constructor        
VariableDensity::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
VariableDensity::Track::Track()
        {
            //pass
        }
        
double VariableDensity::density(GenericTrack *track_input)
        {
            // linear density
            double m;
            m = (rhomax-rhomin)/(track_input->xend-track_input->xini);
            return rhomin+m*(track_input->x-track_input->xini);
        }        
double VariableDensity::ye(GenericTrack *track_input)
        {
            return constant_ye;
        }
        
/*
----------------------------------------------------------------------
         Star CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
Star::Star(double radius_input,double iniradius_input,double inidensity_input)
        {
            name = "Star";
            constant_ye = 0.86;
            
            radius = radius_input;
            initialradius = iniradius_input;
            initialdensity = inidensity_input;
        }
// track constructor        
Star::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Star::Track::Track()
        {
            //pass
        }
double Star::rdensity(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : Starradius
            double r = radius*x;
            if (r <= radius and r> initialradius){
                return initialdensity*exp( - r/initialradius );
            } else if ( r < initialradius){
                return initialdensity;
            } else
                return 0.0;
        }
        
double Star::density(GenericTrack *track_input)
        {
            double rkm = track_input->x/param.km;
            return rdensity(rkm/radius);
        }        
double Star::ye(GenericTrack *track_input)
        {
            return constant_ye;
        }             
        
/*
----------------------------------------------------------------------
         Earth CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
Earth::Earth()
        {
            name = "Earth";
            radius = 6371.0; // [km]
            constant_ye = 0.494;
        }
// track constructor        
Earth::Track::Track(double xini_input, double xend_input,double baseline_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
            baseline = baseline_input;
        }
Earth::Track::Track()
        {
            //pass
        }
        
double Earth::rdensity(double x){
        // Calcula la densidad de la Tierra segun el PREM
        // R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        // Arxiv : 9512364 pag. 23
        // x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        double dne = 0.0;
        double r = radius*x;
        if (r <= 1221.50)
            dne = 13.08850 - 8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150 - 1.26380*x - 3.64260*SQR(x) - 5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650 - 6.47610*x + 5.52830*SQR(x) - 3.08070*SQR(x)*x;
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970 - 1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940 - 8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890 - 3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910 + 0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else if (r>= radius)
            dne=0.0;
        return dne;
        }
        
double Earth::density(Body::Track *track_input)
        {
            Earth::Track* track_earth = static_cast<Earth::Track*>(track_input);
            double xkm = track_earth->x/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-(track_earth->baseline/param.km)*xkm);
            return rdensity(r/radius);
        }        
double Earth::ye(Body::Track *track_input)
        {
            return constant_ye;
        }       
        
/*
----------------------------------------------------------------------
         AGN CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
AGN::AGN(double radius_input,double inidensity_input, double ye_input)
        {
            name = "AGN";
                        
            radius = radius_input;
            initialdensity = inidensity_input;
            constant_ye = ye_input;            
        }
// track constructor        
AGN::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
AGN::Track::Track()
        {
            //pass
        }
double AGN::rdensity(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : AGNradius
            return initialdensity*exp( -x * 10);
        }
        
double AGN::density(Body::Track *track_input)
        {   
            double rparsec = track_input->x/param.parsec;
            return rdensity(rparsec/radius);
        }        
double AGN::ye(Body::Track *track_input)
        {
            return constant_ye;
        }
        
/*
----------------------------------------------------------------------
         SUN CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
Sun::Sun()
        {
            name = "Sun";
            //radius = 694439.0;
            radius = 695980.0;
            
            // import sun model
            sun_model = quickread(SUN_MODEL_LOCATION);
            arraysize = sun_model.size();
            
            sun_radius = new double[arraysize];
            sun_density = new double[arraysize];
            sun_xh = new double[arraysize];
            
            for (int i=0; i < arraysize;i++){
                sun_radius[i] = sun_model[i][1];
                sun_density[i] = sun_model[i][3];
                sun_xh[i] = sun_model[i][6];
            }
            
            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,sun_radius,sun_density,arraysize);
            
            inter_rxh = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_rxh_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_rxh,sun_radius,sun_xh,arraysize);
            
            /*
            // import electron number density for this model
            sun_model_nele = quickread(SUN_MODEL_NELECTRON_LOCATION);
            arraysize_2 = sun_model_nele.size();
            sun_nele_radius = new double[arraysize_2];
            sun_nele = new double[arraysize_2];
            
            double Na = 6.0221415e+23;
            
            for (int i=0; i < arraysize;i++){
                sun_nele_radius[i] = sun_model_nele[i][0];
                sun_nele[i] = exp(sun_model_nele[i][1])*Na;
            }
            
            inter_nele = gsl_spline_alloc(gsl_interp_cspline,arraysize_2);
            inter_nele_accel = gsl_interp_accel_alloc();
            gsl_spline_init (inter_nele,sun_nele_radius,sun_nele,arraysize_2);
            */
            
        }
// track constructor        
Sun::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Sun::Track::Track()
        {
            //pass
        }
double Sun::rdensity(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_density[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_density,x,inter_density_accel);
            }
        }
        
double Sun::rxh(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_xh[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_rxh,x,inter_rxh_accel);
            }
        }
        
double Sun::density(GenericTrack *track_input)
        {
            double r = track_input->x/(radius*param.km);
            return rdensity(r);
        }        
double Sun::ye(GenericTrack *track_input)
        {
            double r = track_input->x/(radius*param.km);
            return 0.5*(1.0+rxh(r));
        }
        
/*
----------------------------------------------------------------------
         SUN ASNU CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
SunASnu::SunASnu()
        {
            name = "SunASnu";
            radius = 694439.0;
    
            sun_model = quickread(SUN_MODEL_LOCATION);
            arraysize = sun_model.size();
            
            sun_radius = new double[arraysize];
            sun_density = new double[arraysize];
            sun_xh = new double[arraysize];
            
            for (int i=0; i < arraysize;i++){
                sun_radius[i] = sun_model[i][1];
                sun_density[i] = sun_model[i][3];
                sun_xh[i] = sun_model[i][6];
            }

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,sun_radius,sun_density,arraysize);
            
            inter_rxh = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_rxh_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_rxh,sun_radius,sun_xh,arraysize);
        }
// track constructor        
SunASnu::Track::Track(double xini_input, double b_impact_input)
        {   
            radius_nu = 694439.0*param.km;
            x = xini_input;
            xini = xini_input;
            b_impact = b_impact_input;
            xend = 2.0*sqrt(SQR(radius_nu)+SQR(b_impact));
        }
SunASnu::Track::Track()
        {
            //pass
        }        
        
double SunASnu::rdensity(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_density[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_density,x,inter_density_accel);
            }
        }
        
double SunASnu::rxh(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_xh[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_rxh,x,inter_rxh_accel);
            }
        }
        
double SunASnu::density(Body::Track *track_input)
        {
            SunASnu::Track* track_sunasnu = static_cast<SunASnu::Track*>(track_input);
            double xkm = track_sunasnu->x/param.km;
            double bkm = track_sunasnu->b_impact/param.km;
            
            double r = sqrt(SQR(radius)+SQR(xkm)-2.0*xkm*sqrt(SQR(radius)-SQR(bkm)))/radius;
            
            return rdensity(r);
        }        
double SunASnu::ye(Body::Track *track_input)
        {
            SunASnu::Track* track_sunasnu = static_cast<SunASnu::Track*>(track_input);
            double xkm = track_sunasnu->x/param.km;
            double bkm = track_sunasnu->b_impact/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-2.0*xkm*sqrt(SQR(radius)-SQR(bkm)))/radius;
            return 0.5*(1.0+rxh(r));
        }
        
/*
----------------------------------------------------------------------
         EARTHATM CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
EarthAtm::EarthAtm()
        {
            name = "EarthAtm";
            radius = 6371.0;
        }
// track constructor        
EarthAtm::Track::Track(double phi_input)
        {   
            radius_nu = 6371.0*param.km;
            atmheight = 20.0*param.km;
            
            phi = phi_input;
            cosphi = cos(phi);
            /*
            if(cosphi<=0.0){
                L = 2.0*radius_nu*abs(cosphi);
            } else {
                L = atmheight/abs(cosphi);
            }
            */
            
            double R = radius_nu;
            double r = atmheight;
            double mm = tan(phi);
            if(cosphi<=0.0){
                L = sqrt(((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*r*R + 2.0*R*
                    (R + sqrt((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*
                    r*R + R*R)))/(1.0 + mm*mm));
		LE =(2.*R)/sqrt(1.0+mm*mm); 
            } else {
                L = sqrt(((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*r*R + 2.0*R*
                    (R - sqrt((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*
                    r*R + R*R)))/(1.0 + mm*mm));
		LE = 0.0; 
            }
            
            x = 0.0;
            xini = 0.0;
            xend = L;
            
            #ifdef EarthAtm_DEBUG
                cout << "== Init Track ==" << endl;
                cout << " phi = " << phi <<
                ", cos(phi) = " << cosphi <<
                ", L = " << L/param.km << endl;
                cout << "==" << endl;
            #endif
            
        }
EarthAtm::Track::Track()
        {
            //pass
        }
double EarthAtm::rdensity(double x){
        // Calcula la densidad de la Tierra segun el PREM
        // R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        // Arxiv : 9512364 pag. 23
        // x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        double dne = 0.0;
        double r = radius*x;
        if (r <= 1221.50)
            dne = 13.08850 - 8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150 - 1.26380*x - 3.64260*SQR(x) - 5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650 - 6.47610*x + 5.52830*SQR(x) - 3.08070*SQR(x)*x;
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970 - 1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940 - 8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890 - 3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910 + 0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368.0)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else if (r>= radius)
            dne=0.0;
        return dne;
        }
        
double EarthAtm::density(Body::Track *track_input)
        {
            EarthAtm::Track* track_earthatm = static_cast<EarthAtm::Track*>(track_input);
            double xkm = track_earthatm->x/param.km;
            
            if(track_earthatm->cosphi <= 0.0){
                double r = sqrt(SQR(radius) + SQR(xkm) - (track_earthatm->L/param.km)*xkm);
                
                #ifdef EarthAtm_DEBUG
                
                cout << "r : " << r << " L : " << (track_earthatm->L/param.km)
                     << " x : " << xkm << " R : " << radius 
                << endl;
                #endif
                
                return rdensity(r/radius);
            } else {
                double h = xkm*track_earthatm->cosphi;
                double h0 = 25.0;
                return 1.05*exp(h/h0);
            }
        }        
double EarthAtm::ye(Body::Track *track_input)
        {
            //EarthAtm::Track* track_earthatm = static_cast<EarthAtm::Track*>(track_input);
            return 0.494;
        }      
