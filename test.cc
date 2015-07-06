//g++ test.cc -std=c++0x -I.. -lfftw3 -L/opt/local/lib -I/opt/local/include -I/opt/local/include/netpbm  -L.. -lnPhysImageF && ./a.out

#include "nPhysImageF.h"
#include "nPhysFormats.h"
#include "nPhysMaths.h"
#include <iomanip>

#include <glob.h>
#include <vector>
#include <string>


using namespace std;


inline vector<string> glob(const string& pat){
    using namespace std;
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

int main(int argc, char **argv) {		
	
    //string dirname="/Volumes/groups/Groupes_recherches/PHYHDEL/Experiments/14-NS-F06_decaying_shock/SOP/calib";
    string dirname="/Volumes/duetera/calib";
//     string dirname="/Volumes/netapp/Groups/Groupes_recherches/PHYHDEL/Experiments/14-NS-F06_decaying_shock/SOP/calib";

    vector<pair<string,double> > data;
    
    data.push_back(make_pair("calib_shutter_closed_",0));
    data.push_back(make_pair("calib_000001_ns",1.179));
    data.push_back(make_pair("calib_000002_ns",1.995));
    data.push_back(make_pair("calib_000005_ns",5.918));
    data.push_back(make_pair("calib_000010_ns",10.947));
    data.push_back(make_pair("calib_000010_ns2_",10.947));
    data.push_back(make_pair("calib_000050_ns",49.100));
    data.push_back(make_pair("calib_000100_ns",98.526));
    data.push_back(make_pair("calib_000200_ns",191.025));
    data.push_back(make_pair("calib_000500_ns",507.852));
    data.push_back(make_pair("calib_001000_ns",1003.027));
    data.push_back(make_pair("calib_002000_ns",1990.008));
    data.push_back(make_pair("calib_005000_ns",5076.382));
    data.push_back(make_pair("calib_010000_ns",9919.118));
    data.push_back(make_pair("calib_020000_ns",19816.223));
    data.push_back(make_pair("calib_050000_ns",48585.712));
    data.push_back(make_pair("calib_100000_ns",97750.030));
    data.push_back(make_pair("calib_200000_ns",197578.876));
    data.push_back(make_pair("calib_500000_ns",499607.410));
    data.push_back(make_pair("calib_1000000_ns",1001235.518));

//     data.push_back("calib_0000001_ns_450_40nm");
//     data.push_back("calib_0000002_ns_450_40nm");
//     data.push_back("calib_0000005_ns_450_40nm");
//     data.push_back("calib_0000010_ns_450_40nm");
//     data.push_back("calib_0000020_ns_450_40nm");
//     data.push_back("calib_0000050_ns_450_40nm");
//     data.push_back("calib_0000100_ns_450_40nm");
//     data.push_back("calib_0000500_ns_450_40nm");
//     data.push_back("calib_0002000_ns_450_40nm");
//     data.push_back("calib_0010000_ns_450_40nm");
//     data.push_back("calib_0050000_ns_450_40nm");
//     data.push_back("calib_0200000_ns_450_40nm");
//     data.push_back("calib_0000001_ns_rose");
//     data.push_back("calib_0000010_ns_rose");
//     data.push_back("calib_0010000_ns_rose");
//     data.push_back("calib_0100000_ns_rose");
//     data.push_back("calib_0000001_ns_rose_gain30_");
//     data.push_back("calib_0000010_ns_rose_gain30_");
//     data.push_back("calib_0000100_ns_rose_gain30_");
//     data.push_back("calib_0001000_ns_rose_gain30_");
//     data.push_back("calib_0010000_ns_rose_gain30_");
//     data.push_back("calib_0100000_ns_rose_gain30_");
//     data.push_back("calib_1000000_ns_rose_gain30_");
//     data.push_back("calib_0100000_ns_30gain");
//     data.push_back("calib_1000000_ns_30gain");
//     data.push_back("calib_0000010_ns_600nm");
//     data.push_back("calib_0100000_ns_600nm");
//     data.push_back("calib_1000000_ns_2_5V");
    
    for (int n=0;n<data.size();n++) {
        double sum=0.0;
        double sumstdev=0.0;
        string globstring=data[n].first+"???.img";
        
        vector<string> fnames=glob(dirname+"/"+globstring);
        for (int i=0;i<fnames.size();i++) {
        
	        nPhysD myimage=physDouble_img(fnames[i]);
	        
    
            int dx=myimage.getW();
            int dy=myimage.getH();
    
            nPhysD una=myimage.sub(0.25*dx,0,0.5*dx,dy);
            
            double surf=una.getSurf();
            
            double mean=phys_sum_points(una)/surf;
            phys_subtract(una,mean);
            double stdev = sqrt(phys_sum_square_points(una))/surf;
            
            sum+=mean;
            sumstdev+=stdev;
            
//             cout << "# " << setw(4)<< n << setw(4)<< i+1 << setw(14) << mean;
//             cout << setw(14) << stdev << setw(14) << data[n].second;
// 
//             cout << " " << fnames[i].substr(dirname.size()+1) <<  endl; 

            
        }
        cout << n << " " << fnames.size() << " " << sum/fnames.size()  << " " << sumstdev/fnames.size() << " " << data[n].second << " " << globstring << endl; 
    }
    
}