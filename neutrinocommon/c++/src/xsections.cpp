#include "xsections.h"

NeutrinoCrossSections::NeutrinoCrossSections(double Emin_in,double Emax_in,int div_in){

       Emin = Emin_in;
       Emax = Emax_in;
       div = div_in;
       
       bool use_existing_files = true;
       bool all_files_exist;

       //string root = "/Users/carguelles/Programs/neutrinocommons/neutrinocommon/c++/data/generate/";
       string root = XSECTION_LOCATION ;
       string filename_format = "_"+toString(Emin)+"_"+toString(Emax)+"_"+toString(div)+".dat";
       
       string filename_dsde_CC = root+"dsde_CC"+filename_format;
       string filename_dsde_NC = root+"dsde_NC"+filename_format;
       string filename_sigma_CC = root+"sigma_CC"+filename_format;
       string filename_sigma_NC = root+"sigma_NC"+filename_format;
    
       #ifndef NO_STD_OUTPUT
       std::cout << "== Trying to load NuSigma Cross Sections ==" << std::endl;
       #endif
       // check if files exist for this energies and divisions
       if(
          fexists(filename_dsde_CC) and
          fexists(filename_dsde_NC) and
          fexists(filename_sigma_CC) and
          fexists(filename_sigma_NC) and
          use_existing_files
          )
       {
            all_files_exist = true;
       } else {
            all_files_exist = false;
       }
       
       if(all_files_exist) {
        // read files
        #ifndef NO_STD_OUTPUT
            std::cout << " XS found. Loading. " << std::endl;
        #endif
            dsde_CC_tbl = quickread(filename_dsde_CC);
            dsde_NC_tbl = quickread(filename_dsde_NC);
            sigma_CC_tbl = quickread(filename_sigma_CC);
            sigma_NC_tbl = quickread(filename_sigma_NC);
       } else {
        #ifndef NO_STD_OUTPUT
            std::cout << " XS not found. Generating." << std::endl;
        #endif
       // generate and save
        vector<double> E_range_GeV = logspace(Emin/1.0e9,Emax/1.0e9,div);
        int e_size = E_range_GeV.size();

        int target = 2; // isospin
        int cc = 1;
        int nc = 0;
        
        for (int e1 = 0 ; e1 < e_size ; e1 ++){
            double Enu1 = E_range_GeV[e1];
            for (int e2 = 0 ; e2 < e_size ; e2 ++){
                double Enu2 = E_range_GeV[e2];
                Row dsde_cc;
                Row dsde_nc;

                dsde_cc.push_back(Enu1);
                dsde_nc.push_back(Enu1);
                dsde_cc.push_back(Enu2);
                dsde_nc.push_back(Enu2);
                
                for (int flavor = 1 ; flavor <= 6; flavor ++){
                    dsde_cc.push_back(dsde_(&Enu1,&Enu2,&flavor,&target,&cc));
                    dsde_nc.push_back(dsde_(&Enu1,&Enu2,&flavor,&target,&nc));
                }
                dsde_CC_tbl.push_back(dsde_cc);
                dsde_NC_tbl.push_back(dsde_nc);
            }
            Row row_cc;
            Row row_nc;
            
            row_cc.push_back(Enu1);
            row_nc.push_back(Enu1);
            
            for (int flavor = 1 ; flavor <= 6; flavor ++){
                row_cc.push_back(sigma_(&Enu1,&flavor,&target,&cc));
                row_nc.push_back(sigma_(&Enu1,&flavor,&target,&nc));
            }
            sigma_CC_tbl.push_back(row_cc);
            sigma_NC_tbl.push_back(row_nc);
            
            #ifndef NO_STD_OUTPUT
            if(e1 % 10){
                cout << e1 << " out of " << e_size << " done." << endl;
            }
            #endif
        }
        
        #ifdef NCS_Init_DEBUG
            std::cout << "== SIGMA-CC-Tbl ==" << std::endl;
            for(unsigned int i = 0 ; i < sigma_CC_tbl.size(); i++){
                for(unsigned int j = 0 ; j < sigma_CC_tbl.size(); j++){
                    cout << sigma_CC_tbl[i][j] << " ";
                }
                cout << endl; 
            }
            cout << "==" << endl;
        #endif
        // saving files
        quickwrite(filename_dsde_CC,dsde_CC_tbl);
        quickwrite(filename_dsde_NC,dsde_NC_tbl);
        quickwrite(filename_sigma_CC,sigma_CC_tbl);
        quickwrite(filename_sigma_NC,sigma_NC_tbl);
       }
       #ifndef NO_STD_OUTPUT
       std::cout << "XS Generated/Loaded." << std::endl;
       #endif
};

// differential cross sections

double NeutrinoCrossSections::dsde_CC(int i_enu, int i_ele, int flv, int neutype){
    if (i_ele > i_enu) {
        return 0.0;
    } else {
        int ii = i_enu*(div+1) + (i_ele);        
        return dsde_CC_tbl[ii][2+2*flv+neutype];
    }
};

double NeutrinoCrossSections::dsde_NC(int i_enu, int i_ele, int flv, int neutype){
    if (i_ele > i_enu) {
        return 0.0;
    } else {
        int ii = i_enu*(div+1) + i_ele;
        return (double) dsde_NC_tbl[ii][2+neutype];
    }
};

//total cross sections

double NeutrinoCrossSections::sigma_CC(int i_enu, int flv, int neutype){    
    return sigma_CC_tbl[i_enu][1+2*flv+neutype];
};

double NeutrinoCrossSections::sigma_NC(int i_enu, int flv, int neutype){
    return sigma_NC_tbl[i_enu][1+2*flv+neutype];
};