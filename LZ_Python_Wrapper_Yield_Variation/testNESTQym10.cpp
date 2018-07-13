/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   testNEST.cpp
 * Author: brodsky3
 *
 * Created on August 1, 2017, 1:03 PM
 */

/*
 * Original testNEST.cpp has been modified with python wrapper so that after compiling the NEST can be imported into python space as a module
 * Author for python wrpper: D.Q. Huang
 * Date: June 8, 2018
 * version 1.0
 * version 2.0 - 20180615 - add anther class (nest_py_wrapper_mass) for vertization, the nest_py_wrapper class has now been renamed as nest_py_wrapper_std
 * version 2.1 - 20180617
 *         - nest_py_wrapper_mass is updated such that when either position or field are not inputed, use default random posiiton or random field.
 *         - no longer need to input a seed list. now either no seed input or just one seed value
 *         - added reset_input function
 * version 2.3 - 20180702
 *         - add posiion x,y,z of interaction into result struct
 *         - add field into result struct
 * Parts has been added or modifed based on the original testNEST.cpp are indicated by
 * "// this part to be added for python wrapper --->" and "// this part to be added for python wrapper <---"
 * where "--->" indicates the beginning of added part and "<----" indicates the end of the added part
 * The wrapper contains five parts:
 *  1. std_vector_to_py_list function, which is used to convert c++ vector to python-style list
 *  2. nestOP stuct, which is the output type when run NEST in python
 *  3. nest_py_wrapper,the wrapper class
 *  4. testNEST function, which is a function of nest_py_wrapper class, and the original main function
 *  5. BOOST_PYTHON_MODULE(wrapper), at the very end of file, which wrappes up class and function defined to be callable by python
 * Follow "readme.txt" to learn how to set up boost python and compile the code and run the module in python
 * In case there is a new NEST output(e.g. nph, ne, etc.) that needed to be added, follow the following steps to modify the code (use nph as example)
 *  1. add "boost::python::list nph;" into "nestOP" struct
 *  2. add "vector<int> nph;" into "nestOP nest_py_wrapper_std::testNEST(boost::python::list command_input)" function
 *     add "vector<int> nph;" into "nestOP nest_py_wrapper_mass::testNEST(boost::python::list command_input)" function
 *  3. add "nest_op.nph = (std_vector_to_py_list<int>(nph));" into "nestOP nest_py_wrapper_std::testNEST(boost::python::list command_input)" function
 *     add "nest_op.nph = (std_vector_to_py_list<int>(nph));" into "nestOP nest_py_wrapper_mass::testNEST(int numEvt)" function
 *  4. add "nph" argument into "testNEST(argc, argv, nph, ne, s1_raw, s1_zCorr, s1c_spike, ne_extract, s2_raw, s2_zCorr);" which is in "nestOP nest_py_wrapper::testNEST(boost::python::list command_input)" function
 *  5. add "vector<int>& nph" argument into "int testNEST(int argc, char** argv, vector<int>& nph, vector<int>& ne, vector<double>& s1_raw, vector<double>& s1_zCorr, vector<double>& s1c_spike, vector<int>& ne_extract, vector<double>& s2_raw, vector<double>& s2_zCorr);"
 *  6. add "nph.push_back(quanta.photons);" inside the body of "int testNEST(int argc, char** argv, vector<int>& nph, vector<int>& ne, vector<double>& s1_raw, vector<double>& s1_zCorr, vector<double>& s1c_spike, vector<int>& ne_extract, vector<double>& s2_raw, vector<double>& s2_zCorr);"
 *  7. add ".def_readonly("nph",&nestOP::nph)" into "BOOST_PYTHON_MODULE(wrapper)"
 * In the output
 boost::python::list nph;
 boost::python::list ne;
 boost::python::list s1_n_hits;
 boost::python::list s1_n_phe;
 boost::python::list s1_raw_area_phe;
 boost::python::list s1c_area_phe;
 boost::python::list s1_raw_area_phd;
 boost::python::list s1c_area_phd;
 boost::python::list s1_raw_spike;
 boost::python::list s1c_spike;
 boost::python::list s2_ne_extract;
 boost::python::list s2_n_ph;
 boost::python::list s2_n_hits;
 boost::python::list s2_n_phe;
 boost::python::list s2_raw_area_phe;
 boost::python::list s2c_area_phe;
 boost::python::list s2_raw_area_phd;
 boost::python::list s2c_area_phd;
 boost::python::list deposite_energy;
 boost::python::list pos_x_mm
 boost::python::list pos_y_mm
 boost::python::list pos_z_mm
 boost::python::list field_V_per_cm
 */

#include "NEST.hh"
#include "TestSpectra.hh"
#include "analysis.hh"

#include "DetectorExample_XENON10.hh"
#include "LZ_Detector.hh"

// this part to be added for python wrapper --->
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
// this part to be added for python wrapper <---

using namespace std;
using namespace NEST;

/*
 * 
 */

double band[NUMBINS_MAX][6], energies[3];
vector<vector<double>> GetBand ( vector<double> S1s, vector<double> S2s, bool resol );
void GetEnergyRes ( vector<double> Es );

// this part to be added for python wrapper --->
template<class T> boost::python::list std_vector_to_py_list(const std::vector<T>& v)
{
    boost::python::object get_iter = boost::python::iterator<std::vector<T> >();
    boost::python::object iter = get_iter(v);
    boost::python::list l(iter);
    return l;
}

struct nestOP{
    boost::python::list nph;
    boost::python::list ne;
    boost::python::list s1_n_hits;
    boost::python::list s1_n_phe;
    boost::python::list s1_raw_area_phe;
    boost::python::list s1c_area_phe;
    boost::python::list s1_raw_area_phd;
    boost::python::list s1c_area_phd;
    boost::python::list s1_raw_spike;
    boost::python::list s1c_spike;
    boost::python::list s2_ne_extract;
    boost::python::list s2_n_ph;
    boost::python::list s2_n_hits;
    boost::python::list s2_n_phe;
    boost::python::list s2_raw_area_phe;
    boost::python::list s2c_area_phe;
    boost::python::list s2_raw_area_phd;
    boost::python::list s2c_area_phd;
    boost::python::list deposite_energy;
    boost::python::list pos_x_mm;
    boost::python::list pos_y_mm;
    boost::python::list pos_z_mm;
    boost::python::list field_V_per_cm;
};

class nest_py_wrapper{
protected:
    int argc;
    char** argv;
    nestOP nest_op;
    int print_result = 1;
    int testNEST(int argc, char** argv, int print_result,
                 vector<int>& nph, vector<int>& ne, vector<int>& s1_n_hits, vector<int>& s1_n_phe,
                 vector<double>& s1_raw_area_phe, vector<double>& s1c_area_phe, vector<double>& s1_raw_area_phd, vector<double>& s1c_area_phd, vector<double>& s1_raw_spike, vector<double>& s1c_spike,
                 vector<int>& s2_ne_extract, vector<double>& s2_n_ph, vector<double>& s2_n_hits, vector<double>& s2_n_phe,
                 vector<double>& s2_raw_area_phe, vector<double>& s2c_area_phe, vector<double>& s2_raw_area_phd, vector<double>& s2c_area_phd,
                 vector<double>& deposite_energy, vector<double>& pos_x_mm, vector<double>& pos_y_mm, vector<double>& pos_z_mm, vector<double>& field_V_per_cm);
public:
    void print_result_or_not(int yn){print_result = yn;}
};

class nest_py_wrapper_mass: public nest_py_wrapper{
private:
    string type;
    boost::python::list energy;
    boost::python::list x;
    boost::python::list y;
    boost::python::list z;
    boost::python::list field;
    int seed;
    int seed_flag = 0;
public:
    void reset_input(){
        type = ""; seed_flag = 0; energy = boost::python::list(); x = boost::python::list(); y = boost::python::list(); z = boost::python::list(); field = boost::python::list();nest_op = nestOP();
    }
    void input_type_interaction(string type_){type = type_;}
    void input_energy_list(boost::python::list energy_){energy = energy_;}
    void input_x_list(boost::python::list x_){x = x_;}
    void input_y_list(boost::python::list y_){y = y_;}
    void input_z_list(boost::python::list z_){z = z_;}
    void input_field_list(boost::python::list field_){field = field_;}
    bool set_seed(int seed_){seed = seed_;seed_flag = 1;return true;}
    nestOP testNEST(int numEvt);
};

nestOP nest_py_wrapper_mass::testNEST(int numEvt){
    vector<int> nph;
    vector<int> ne;
    vector<int> s1_n_hits;
    vector<int> s1_n_phe;
    vector<double> s1_raw_area_phe;
    vector<double> s1c_area_phe;
    vector<double> s1_raw_area_phd;
    vector<double> s1c_area_phd;
    vector<double> s1_raw_spike;
    vector<double> s1c_spike;
    vector<int> s2_ne_extract;
    vector<double> s2_n_ph;
    vector<double> s2_n_hits;
    vector<double> s2_n_phe;
    vector<double> s2_raw_area_phe;
    vector<double> s2c_area_phe;
    vector<double> s2_raw_area_phd;
    vector<double> s2c_area_phd;
    vector<double> deposite_energy;
    vector<double> pos_x_mm;
    vector<double> pos_y_mm;
    vector<double> pos_z_mm;
    vector<double> field_V_per_cm;
    
    //vector<boost::python::list> result;
    
    float energy_;
    float x_;
    float y_;
    float z_;
    float field_;
    
    int l_energy = len(energy);
    int l_x = len(x);
    int l_y = len(y);
    int l_z = len(z);
    int l_field = len(field);
    if(l_energy!=1 && l_energy!=numEvt){cout<<"Error: make sure the energy list length is either 1 or num of evts"<<endl;return nest_op;}
    if(l_x!=0 && l_x!=1 && l_x!=numEvt){cout<<"Error: make sure the x list length is either 0, 1 or num of evts"<<endl;return nest_op;}
    if(l_y!=0 && l_y!=1 && l_y!=numEvt){cout<<"Error: make sure the y list length is either 0, 1 or num of evts"<<endl;return nest_op;}
    if(l_z!=0 && l_z!=1 && l_z!=numEvt){cout<<"Error: make sure the z list length is either 0, 1 or num of evts"<<endl;return nest_op;}
    if(l_field!=0 && l_field!=1 && l_field!=numEvt){cout<<"Error: make sure the field list length is either 0, 1 or num of evts"<<endl;return nest_op;}
    //  boost::python::object print =
    //  boost::python::import("__main__").attr("__builtins__").attr("print");
    if(seed_flag == 0){srand(time(NULL));}
    else {srand(seed);}
    argc = 8;
    argv = (char**)std::malloc((argc)*sizeof(char *));
    std:: string str_[argc];
    str_[0] = "ignore";
    argv[0] = const_cast<char *>(str_[0].c_str());
    str_[1] = "1";
    argv[1] = const_cast<char *>(str_[1].c_str());
    str_[2] = type;
    argv[2] = const_cast<char *>(str_[2].c_str());
    for(int i = 0; i < numEvt; i++){
        if(l_energy == 1)
            energy_ = boost::python::extract<float>(energy[0]);
        else
            energy_ = boost::python::extract<float>(energy[i]);
        if(l_x == 0 || l_y == 0 || l_z == 0){}
        else{
            if(l_x == 1)
                x_ = boost::python::extract<float>(x[0]);
            else
                x_ = boost::python::extract<float>(x[i]);
            if(l_y == 1)
                y_ = boost::python::extract<float>(y[0]);
            else
                y_ = boost::python::extract<float>(y[i]);
            if(l_z == 1)
                z_ = boost::python::extract<float>(z[0]);
            else
                z_ = boost::python::extract<float>(z[i]);
        }
        if(l_field == 0){}
        else{
            if(l_field == 1)
                field_ = boost::python::extract<float>(field[0]);
            else
                field_ = boost::python::extract<float>(field[i]);
        }
        str_[3] = std::to_string(energy_);
        argv[3] = const_cast<char *>(str_[3].c_str());
        str_[4] = std::to_string(energy_);
        argv[4] = const_cast<char *>(str_[3].c_str());
        if(l_field == 0)
        {str_[5] = "-1";}
        else
        {str_[5] = std::to_string(field_);}
        argv[5] = const_cast<char *>(str_[5].c_str());
        if(l_x == 0 || l_y == 0 || l_z == 0)
        {str_[6] = "-1";}
        else
        {str_[6] = std::to_string(x_) + string(",") + std::to_string(y_) + string(",") + std::to_string(z_);}
        argv[6] = const_cast<char *>(str_[6].c_str());
        str_[7] = std::to_string(rand()).c_str();
        argv[7] = const_cast<char *>(str_[7].c_str());
        nest_py_wrapper::testNEST(argc, argv, print_result, nph, ne, s1_n_hits, s1_n_phe,
                                  s1_raw_area_phe, s1c_area_phe, s1_raw_area_phd, s1c_area_phd, s1_raw_spike, s1c_spike,
                                  s2_ne_extract, s2_n_ph, s2_n_hits, s2_n_phe,
                                  s2_raw_area_phe, s2c_area_phe, s2_raw_area_phd, s2c_area_phd,
                                  deposite_energy, pos_x_mm, pos_y_mm, pos_z_mm, field_V_per_cm);
        
    }
    free(argv);
    nest_op.nph = (std_vector_to_py_list<int>(nph));
    nest_op.ne = (std_vector_to_py_list<int>(ne));
    nest_op.s1_n_hits = (std_vector_to_py_list<int>(s1_n_hits));
    nest_op.s1_n_phe = (std_vector_to_py_list<int>(s1_n_phe));
    nest_op.s1_raw_area_phe = (std_vector_to_py_list<double>(s1_raw_area_phe));
    nest_op.s1c_area_phe = (std_vector_to_py_list<double>(s1c_area_phe));
    nest_op.s1_raw_area_phd = (std_vector_to_py_list<double>(s1_raw_area_phd));
    nest_op.s1c_area_phd = (std_vector_to_py_list<double>(s1c_area_phd));
    nest_op.s1_raw_spike = (std_vector_to_py_list<double>(s1_raw_spike));
    nest_op.s1c_spike = (std_vector_to_py_list<double>(s1c_spike));
    nest_op.s2_ne_extract = (std_vector_to_py_list<int>(s2_ne_extract));
    nest_op.s2_n_ph = (std_vector_to_py_list<double>(s2_n_ph));
    nest_op.s2_n_hits = (std_vector_to_py_list<double>(s2_n_hits));
    nest_op.s2_n_phe = (std_vector_to_py_list<double>(s2_n_phe));
    nest_op.s2_raw_area_phe = (std_vector_to_py_list<double>(s2_raw_area_phe));
    nest_op.s2c_area_phe = (std_vector_to_py_list<double>(s2c_area_phe));
    nest_op.s2_raw_area_phd = (std_vector_to_py_list<double>(s2_raw_area_phd));
    nest_op.s2c_area_phd = (std_vector_to_py_list<double>(s2c_area_phd));
    nest_op.deposite_energy = (std_vector_to_py_list<double>(deposite_energy));
    nest_op.pos_x_mm = (std_vector_to_py_list<double>(pos_x_mm));
    nest_op.pos_y_mm = (std_vector_to_py_list<double>(pos_y_mm));
    nest_op.pos_z_mm = (std_vector_to_py_list<double>(pos_z_mm));
    nest_op.field_V_per_cm = (std_vector_to_py_list<double>(field_V_per_cm));
    return nest_op;
}

class nest_py_wrapper_std: public nest_py_wrapper{
private:
public:
    nestOP testNEST(boost::python::list command_input);
};

nestOP nest_py_wrapper_std::testNEST(boost::python::list command_input){
    //int argc;
    //char** argv;
    vector<int> nph;
    vector<int> ne;
    vector<int> s1_n_hits;
    vector<int> s1_n_phe;
    vector<double> s1_raw_area_phe;
    vector<double> s1c_area_phe;
    vector<double> s1_raw_area_phd;
    vector<double> s1c_area_phd;
    vector<double> s1_raw_spike;
    vector<double> s1c_spike;
    vector<int> s2_ne_extract;
    vector<double> s2_n_ph;
    vector<double> s2_n_hits;
    vector<double> s2_n_phe;
    vector<double> s2_raw_area_phe;
    vector<double> s2c_area_phe;
    vector<double> s2_raw_area_phd;
    vector<double> s2c_area_phd;
    vector<double> deposite_energy;
    vector<double> pos_x_mm;
    vector<double> pos_y_mm;
    vector<double> pos_z_mm;
    vector<double> field_V_per_cm;
    
    //vector<boost::python::list> result;
    int l = len(command_input);
    argc = l+1;
    argv = (char**)std::malloc((argc)*sizeof(char *));
    std::string str_[argc];
    for(int i = 0; i < argc; i++){
        if(i==0){argv[i] = const_cast<char *>("ignore");}
        else{
            str_[i] = boost::python::extract<std::string>(command_input[i-1]);
            argv[i] = const_cast<char *>(str_[i].c_str());
        }
    }
    nest_py_wrapper::testNEST(argc, argv, print_result, nph, ne, s1_n_hits, s1_n_phe,
                              s1_raw_area_phe, s1c_area_phe, s1_raw_area_phd, s1c_area_phd, s1_raw_spike, s1c_spike,
                              s2_ne_extract, s2_n_ph, s2_n_hits, s2_n_phe,
                              s2_raw_area_phe, s2c_area_phe, s2_raw_area_phd, s2c_area_phd,
                              deposite_energy, pos_x_mm, pos_y_mm, pos_z_mm, field_V_per_cm);
    
    
    free(argv);
    
    // define the print function
    boost::python::object print =
    boost::python::import("__main__").attr("__builtins__").attr("print");
    
    nest_op.nph = (std_vector_to_py_list<int>(nph));
    nest_op.ne = (std_vector_to_py_list<int>(ne));
    nest_op.s1_n_hits = (std_vector_to_py_list<int>(s1_n_hits));
    nest_op.s1_n_phe = (std_vector_to_py_list<int>(s1_n_phe));
    nest_op.s1_raw_area_phe = (std_vector_to_py_list<double>(s1_raw_area_phe));
    nest_op.s1c_area_phe = (std_vector_to_py_list<double>(s1c_area_phe));
    nest_op.s1_raw_area_phd = (std_vector_to_py_list<double>(s1_raw_area_phd));
    nest_op.s1c_area_phd = (std_vector_to_py_list<double>(s1c_area_phd));
    nest_op.s1_raw_spike = (std_vector_to_py_list<double>(s1_raw_spike));
    nest_op.s1c_spike = (std_vector_to_py_list<double>(s1c_spike));
    nest_op.s2_ne_extract = (std_vector_to_py_list<int>(s2_ne_extract));
    nest_op.s2_n_ph = (std_vector_to_py_list<double>(s2_n_ph));
    nest_op.s2_n_hits = (std_vector_to_py_list<double>(s2_n_hits));
    nest_op.s2_n_phe = (std_vector_to_py_list<double>(s2_n_phe));
    nest_op.s2_raw_area_phe = (std_vector_to_py_list<double>(s2_raw_area_phe));
    nest_op.s2c_area_phe = (std_vector_to_py_list<double>(s2c_area_phe));
    nest_op.s2_raw_area_phd = (std_vector_to_py_list<double>(s2_raw_area_phd));
    nest_op.s2c_area_phd = (std_vector_to_py_list<double>(s2c_area_phd));
    nest_op.deposite_energy = (std_vector_to_py_list<double>(deposite_energy));
    nest_op.pos_x_mm = (std_vector_to_py_list<double>(pos_x_mm));
    nest_op.pos_y_mm = (std_vector_to_py_list<double>(pos_y_mm));
    nest_op.pos_z_mm = (std_vector_to_py_list<double>(pos_z_mm));
    nest_op.field_V_per_cm = (std_vector_to_py_list<double>(field_V_per_cm));
    
    
    // test print
    if(0)
    {
        for(int i = 0; i < len(nest_op.nph);i++){
            print(boost::python::object(nest_op.nph[i]));
        }
    }
    return nest_op;
}
// this part to be added for python wrapper <-------

//this method declaration used to be that of the main method. Simply change the name (and arguments if needed)
// this part to be modified for python wrapper --->
int nest_py_wrapper::testNEST(int argc, char** argv, int print_result,
                              vector<int>& nph, vector<int>& ne, vector<int>& s1_n_hits, vector<int>& s1_n_phe,
                              vector<double>& s1_raw_area_phe, vector<double>& s1c_area_phe, vector<double>& s1_raw_area_phd, vector<double>& s1c_area_phd, vector<double>& s1_raw_spike, vector<double>& s1c_spike,
                              vector<int>& s2_ne_extract, vector<double>& s2_n_ph, vector<double>& s2_n_hits, vector<double>& s2_n_phe,
                              vector<double>& s2_raw_area_phe, vector<double>& s2c_area_phe, vector<double>& s2_raw_area_phd, vector<double>& s2c_area_phd,
                              vector<double>& deposite_energy, vector<double>& pos_x_mm, vector<double>& pos_y_mm, vector<double>& pos_z_mm, vector<double>& field_V_per_cm){
    // this part to be modified for python wrapper <---
  // Instantiate your own VDetector class here, then load into NEST class constructor
	//DetectorExample_XENON10* detector = new DetectorExample_XENON10();
	LZ_Detector* detector = new LZ_Detector();
	
	// Custom parameter modification functions
	//detector->ExampleFunction();

	// Construct NEST class using detector object
	NESTcalc n(detector);

	vector<double> signal1,signal2,signalE, vTable, NuisParam={1.,.9}; //scaling factors, for now just for NR Ly & Qy. But must initialize!
	string position, delimiter, token; size_t loc; int index;
	double g2,pos_x,pos_y,pos_z,r,phi,driftTime, field, vD,vD_middle, atomNum=0,massNum=0, keVee=0.0, origX,origY;
  YieldResult yieldsMax;
  
  if (argc < 7)
    {
      cout << "This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector:" << endl;
      cout << "\t./testNEST numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}" << endl << endl;
      cout << "For 8B, numEvts is kg-days of exposure with everything else same. For WIMPs:" << endl;
      cout << "\t./testNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}" << endl << endl;
      cout << "For cosmic-ray muons or other similar particles with elongated track lengths:" << endl;
      cout << "\t./testNEST numEvts {MIP} LET[MeV*cm^2/gram] x,y-position[mm](Initial) field_drift[V/cm] x,y,z-position[mm](Final) {optional:seed}" << endl << endl;
      return 0;
    }
  
  unsigned long int numEvts = atoi(argv[1]);
  string type = argv[2];
  double eMin = atof(argv[3]);
  if ( eMin == -1. ) eMin = 0.;
  double eMax = atof(argv[4]);
  if ( eMax == -1. && eMin == 0. )
    eMax = 1e2; //the default energy max is 100 keV
  if(eMax==0.){
    cerr << "ERROR: The maximum energy cannot be 0 keV!" << endl;
    return 0; }
  double inField = atof(argv[5]);
  position = argv[6];
  double fPos = atof(argv[6]);
  
  if ( argc == 8 ) {
    if ( atoi(argv[7]) == -1 )
      RandomGen::rndm()->SetSeed( time (NULL) );
    else
      RandomGen::rndm()->SetSeed(atoi(argv[7]));
  }
  
  INTERACTION_TYPE type_num;
	TestSpectra spec;
	if ( type == "NR" || type == "neutron" || type == "-1" ) type_num = NR; //-1: default particle type is also NR
  else if (type == "WIMP")
    {
      if (atof(argv[3])<0.44) { cerr << "WIMP mass too low, you're crazy!" << endl; return 0; }
      type_num = WIMP;
      spec.wimp_spectrum_prep = spec.WIMP_prep_spectrum(atof(argv[3]), E_step);
      numEvts = RandomGen::rndm()->poisson_draw(spec.wimp_spectrum_prep.integral * 1.0 * atof(argv[1]) * atof(argv[4]) / 1e-36);
    }
  else if ( type == "B8" || type == "Boron8" || type == "8Boron" || type == "8B" || type == "Boron-8" )
    {
      type_num = B8;
      numEvts = RandomGen::rndm()->poisson_draw(0.0026 * atof(argv[1]));
    } else if ( type == "DD" || type == "D-D" ) type_num = DD;
  else if ( type == "AmBe" ) type_num = AmBe;
  else if ( type == "Cf" || type == "Cf252" || type == "252Cf" || type == "Cf-252" ) type_num = Cf;
  else if ( type == "ion" || type == "nucleus" || type == "alpha" ) {
    type_num = ion;
    if ( type == "alpha" ) {
      atomNum = 2; massNum = 4;
    }
    else {
      cerr << "Atomic Number: "; cin >> atomNum;
      cerr << "Mass Number: "; cin >> massNum;
    } if ( atomNum == ATOM_NUM ) type_num = NR;
  }
  else if ( type == "gamma" || type == "gammaRay" ||
	    type == "x-ray" || type == "xray" || type == "xRay" || type == "X-ray" || type == "Xray" || type == "XRay" )
    type_num = gammaRay; //includes photo-absorption and electron capture
  else if ( type == "Kr83m" || type == "83mKr" || type == "Kr83" ) type_num = Kr83m;
  else if ( type == "CH3T" || type == "tritium" ) type_num = CH3T;
  else if ( type == "C14" || type == "Carbon14" || type == "14C" ) type_num = C14;
  else if ( type == "beta" || type == "ER" || type == "Compton" || type == "compton" || type == "electron" || type == "e-" ||
	    type == "muon" || type == "MIP" || type == "LIP" || type == "mu" || type == "mu-" )
    type_num = beta; //default electron recoil model
  else {
    cerr << "UNRECOGNIZED PARTICLE TYPE!! VALID OPTIONS ARE:" << endl;
    cerr << "NR or neutron," << endl;
    cerr << "WIMP," << endl;
    cerr << "B8 or Boron8 or 8Boron or 8B or Boron-8," << endl;
    cerr << "DD or D-D," << endl;
    cerr << "AmBe," << endl;
    cerr << "Cf or Cf252 or 252Cf or Cf-252," << endl;
    cerr << "ion or nucleus," << endl;
    cerr << "alpha," << endl;
    cerr << "gamma or gammaRay," << endl;
    cerr << "x-ray or xray or xRay or X-ray or Xray or XRay," << endl;
    cerr << "Kr83m or 83mKr or Kr83," << endl;
    cerr << "CH3T or tritium," << endl;
    cerr << "Carbon14 or 14C or C14," << endl;
    cerr << "beta or ER or Compton or compton or electron or e-, and" << endl;
    cerr << "muon or MIP or LIP or mu or mu-" << endl;
    return 0;
  }
  
  if ( type_num == Kr83m ) {
    if ( eMin == 9.4 && eMax == 9.4 ) {}
    else if ( eMin == 32.1 &&
	      eMax == 32.1 ) {}
    else
      { cerr << "ERROR: For Kr83m, put both energies as 9.4 or both as 32.1 keV please." << endl;
	return 0; }
  }
  
  if ( (eMin < 10. || eMax < 10.) && type_num == gammaRay ) {
    cerr << "WARNING: Typically beta model works better for ER BG at low energies as in a WS." << endl;
    cerr << "ER data is often best matched by a weighted average of the beta & gamma models." << endl;
  }
  
  double rho = n.SetDensity ( detector->get_T_Kelvin(), detector->get_p_bar() ); //cout.precision(12);
  if ( rho < 1. ) detector->set_inGas(true);
  
	// Calculate and print g1, g2 parameters (once per detector)
	vector<double> g2_params = n.CalculateG2(print_result);
	g2 = fabs(g2_params.back()); double g1 = detector->get_g1();
  
	if ( inField == -1. ) {
    vTable = n.SetDriftVelocity_NonUniform(rho, z_step);
    vD_middle = vTable[int(floor(.5*detector->get_TopDrift()/z_step))];
  }
  else vD_middle = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, inField);
  
	if (print_result) {
		cout << "Density = " << rho << " g/mL" << "\t";
		cout << "central vDrift = " << vD_middle << " mm/us\n";
		cout << "\t\t\t\t\t\t\t\t\t\tNegative numbers are flagging things below threshold!   phe=(1+P_dphe)*phd & phd=phe/(1+P_dphe)\n";

		if ( type_num == Kr83m && eMin == 9.4 && eMax == 9.4 )
			fprintf(stdout, "t [ns]\t\tE [keV]\t\tfield [V/cm]\ttDrift [us]\tX,Y,Z [mm]\tNph\tNe-\tS1 [PE or phe]\tS1_3Dcor [phd]\tspikeC(NON-INT)\tNe-Extr\tS2_rawArea [PE]\tS2_3Dcorr [phd]\n");
		else
			fprintf(stdout, "E [keV]\t\tfield [V/cm]\ttDrift [us]\tX,Y,Z [mm]\tNph\tNe-\tS1 [PE or phe]\tS1_3Dcor [phd]\tspikeC(NON-INT)\tNe-Extr\tS2_rawArea [PE]\tS2_3Dcorr [phd]\n");
 	}

  if ( type_num == WIMP ) {
    yieldsMax = n.GetYields(      NR, 25.0, rho, detector->FitEF(0., 0., detector->get_TopDrift()/2.),
				  double(massNum), double(atomNum), NuisParam);
  }
  else if ( type_num == B8 ) {
    yieldsMax = n.GetYields(      NR, 4.00, rho, detector->FitEF(0., 0., detector->get_TopDrift()/2.),
				  double(massNum), double(atomNum), NuisParam);
  }
  else {
    double energyMaximum;
    if ( eMax < 0. ) energyMaximum = 1. / fabs(eMax);
    if ( type_num == Kr83m )
      yieldsMax = n.GetYields(  beta  , energyMaximum, rho, detector->FitEF(0., 0., detector->get_TopDrift()/2.),
			      double(massNum), double(atomNum), NuisParam); //the reason for this: don't do the special Kr stuff when just checking max
    else
      yieldsMax = n.GetYields(type_num, energyMaximum, rho, detector->FitEF(0., 0., detector->get_TopDrift()/2.),
			      double(massNum), double(atomNum), NuisParam);
  }
  if ( (g1*yieldsMax.PhotonYield) > (2.*maxS1) && eMin != eMax )
    cerr << "\nWARNING: Your energy maximum may be too high given your maxS1.\n";
  
  double keV = -999.;
  for (unsigned long int j = 0; j < numEvts; j++) {
    
    if (eMin == eMax && eMin >= 0. && eMax > 0. ) {
      keV = eMin;
    } else {
      switch (type_num) {
      case CH3T:
	keV = spec.CH3T_spectrum(eMin, eMax);
	break;
      case C14:
	keV = spec.C14_spectrum(eMin, eMax);
	break;
      case B8: //normalize this to ~3500 / 10-ton / year, for E-threshold of 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
	keV = spec.B8_spectrum(eMin, eMax);
	break;
      case AmBe: //for ZEPLIN-III FSR from HA (Pal '98)
	keV = spec.AmBe_spectrum(eMin, eMax);
	break;
      case Cf:
	keV = spec.Cf_spectrum(eMin, eMax);
	break;
      case DD:
	keV = spec.DD_spectrum(eMin, eMax);
	break;
      case WIMP:
	{
          keV = spec.WIMP_spectrum(spec.wimp_spectrum_prep, atof(argv[3]));
	}
	break;
      default:
	if ( eMin < 0. ) return 0;
	if ( eMax > 0. )
	  keV = eMin + (eMax - eMin) * RandomGen::rndm()->rand_uniform();
	else { //negative eMax signals to NEST that you want to use an exponential energy spectrum profile
	  if ( eMin == 0. ) return 0;
	  keV=1e100; //eMin will be used in place of eMax as the maximum energy in exponential scenario
	  while ( keV > eMin )
	    keV = eMax * log ( RandomGen::rndm()->rand_uniform() );
	}
	break;
      }
    }
    
    if ( type_num != WIMP && type_num != B8 && eMax > 0. ) {
      if ( keV > eMax ) keV = eMax;
      if ( keV < eMin ) keV = eMin;
    }
    
  Z_NEW:
    if ( fPos == -1. ) { // -1 means default, random location mode
      pos_z = 0. + ( detector->get_TopDrift() - 0. ) * RandomGen::rndm()->rand_uniform(); // initial guess
      r = detector->get_radius() * sqrt ( RandomGen::rndm()->rand_uniform() );
      phi = 2.*M_PI*RandomGen::rndm()->rand_uniform();
      pos_x = r * cos(phi); pos_y = r * sin(phi); origX = pos_x; origY = pos_y;
    }
    else {
      delimiter = ",";
      loc = 0; int i = 0;
      while ( (loc = position.find(delimiter)) != string::npos ) {
	token = position.substr(0,loc);
	if ( i == 0 ) pos_x = stof(token);
	else pos_y = stof(token);
	position.erase(0,loc+delimiter.length());
	i++;
      }
      pos_z = stof(position);
      if ( stof(position) == -1. )
	pos_z = 0. + ( detector->get_TopDrift() - 0. ) * RandomGen::rndm()->rand_uniform();
      if ( stof(token) == -999. ) {
	r = detector->get_radius() * sqrt ( RandomGen::rndm()->rand_uniform() );
	phi = 2.*M_PI*RandomGen::rndm()->rand_uniform();
	pos_x = r * cos(phi); pos_y = r * sin(phi);
      }
      if ( j == 0 ) { origX = pos_x; origY = pos_y; }
    }
    
    if ( inField == -1. ) { // -1 means use poly position dependence
			field = detector->FitEF(pos_x, pos_y, pos_z);
    }
    else field = inField; //no fringing
    
    if ( field <= 0. )
      cerr << "\nWARNING: A LITERAL ZERO FIELD MAY YIELD WEIRD RESULTS. USE A SMALL VALUE INSTEAD.\n";
    if ( field > 12.e3 )
      cerr << "\nWARNING: Your field is >12,000 V/cm. No data out here. Are you sure about this?\n";
    
    if ( inField == -1. ) {
      //for ( int jj = 0; jj < vTable.size(); jj++ ) //DEBUG
      //cerr << double(jj)*z_step << "\t" << vTable[jj] << endl;
      index = int(floor(pos_z/z_step));
      vD = vTable[index];
    }
    else
      vD = n.SetDriftVelocity(detector->get_T_Kelvin(),rho,field);
    driftTime = ( detector->get_TopDrift() - pos_z ) / vD; // (mm - mm) / (mm / us) = us
    if ( inField != -1. && detector->get_dt_min() > ( detector->get_TopDrift() - 0. ) / vD && field >= FIELD_MIN )
      { cerr << "ERROR: dt_min is too restrictive (too large)" << endl; return 0; }
    if ( (driftTime > detector->get_dt_max() || driftTime < detector->get_dt_min()) && (fPos == -1. || stof(position) == -1.) && field >= FIELD_MIN )
      goto Z_NEW;
    if ( detector->get_dt_max() > (detector->get_TopDrift()-0.)/vD && !j && field >= FIELD_MIN )
      { cerr << "WARNING: dt_max is greater than max possible" << endl; }
    
    // The following should never happen: this is simply a just-in-case code-block dealing with user error
    if ( pos_z <= 0. ) {
      cerr << "ERROR: unphysically low Z coordinate (vertical axis of detector) of " << pos_z << " mm" << endl;
      return 0;
    }
    if ( (pos_z > (detector->get_TopDrift()+z_step) || driftTime < 0.0) && field >= FIELD_MIN ) {
      cerr << "ERROR: unphysically big Z coordinate (vertical axis of detector) of " << pos_z << " mm" << endl;
      return 0;
    }
    
    YieldResult yields; QuantaResult quanta;
    if ( type == "muon" || type == "MIP" || type == "LIP" || type == "mu" || type == "mu-" ) {
      double xi = -999., yi = -999.;
      if ( atof(argv[4]) == -1. ) {
	r = detector->get_radius()
	  * sqrt ( RandomGen::rndm()->rand_uniform() );
	phi = 2.*M_PI*RandomGen::rndm()->rand_uniform();
	xi = r * cos(phi); yi = r * sin(phi);
      }
      else {
	position = argv[4];
	delimiter = ",";
	loc = 0; int ii = 0;
	while ( (loc = position.find(delimiter)) != string::npos ) {
	  token = position.substr(0,loc);
	  if ( ii == 0 ) xi = stof(token);
	  else yi = stof(token);
	  position.erase(0,loc+delimiter.length());
	  ii++;
	}
	yi = stof(position);
      }
      double dEOdx = eMin, eStep = dEOdx * rho * z_step * 1e2, refEnergy = 1e6; keV = 0.;
      int Nph = 0, Ne = 0;
      double xx = xi, yy = yi, zz = detector->get_TopDrift();
      double xf = pos_x; double yf = pos_y;
      double distance = sqrt(pow(xf-xi,2.)+pow(yf-yi,2.)+pow(detector->get_TopDrift(),2.));
      double norm[3];
      norm[0] = ( xf - xi ) / distance;
      norm[1] = ( yf - yi ) / distance;
      norm[2] = -detector->get_TopDrift() / distance; //have not yet tested muons which leave before hitting Z=0, would have to modify code here
      while ( zz > 0. && sqrt(pow(xx,2.)+pow(yy,2.)) < detector->get_radius() ) { //stop making S1 and S2 if particle exits Xe vol
	yields = n.GetYields ( beta, refEnergy, rho, detector->FitEF ( xx, yy, zz ), double(massNum), double(atomNum), NuisParam );
	quanta = n.GetQuanta ( yields, rho );
	Nph+= quanta.photons * (eStep/refEnergy);
	index = int(floor(zz/z_step));
	if ( index >= vTable.size() )
	  index = vTable.size() - 1;
	vD = vTable[index];
	driftTime = ( detector->get_TopDrift() - zz ) / vD;
	if ( pos_z >= detector->get_cathode() )
	  Ne += quanta.electrons*(eStep/refEnergy)*exp(-driftTime/detector->get_eLife_us());
	keV+= eStep;
	xx += norm[0] * z_step;
	yy += norm[1] * z_step;
	zz += norm[2] * z_step; //cout << xx << " " << yy << " " << zz << endl;
      }
      quanta.photons = Nph; quanta.electrons = Ne;
      pos_z = detector->get_TopDrift() / 2.; driftTime = 0.00;
      vD = vD_middle; //approximate things not already done right in loop as middle of detector since muon traverses whole length
      pos_x = .5*(xi+xf); pos_y = .5*(yi+yf);
      field = detector->FitEF(pos_x,pos_y,pos_z);
    }
    else {
      yields = n.GetYields(type_num,keV,rho,field,
			   double(massNum),double(atomNum),NuisParam);
      quanta = n.GetQuanta(yields,rho);
    }
    
    double Nphd_S2 = g2 * quanta.electrons * exp(-driftTime/detector->get_eLife_us());
    if ( !MCtruthPos && Nphd_S2 > PHE_MIN ) {
      vector<double> xySmeared(2); xySmeared = n.xyResolution ( origX, origY, Nphd_S2 );
      pos_x = xySmeared[0];
      pos_y = xySmeared[1];
    }

    vector<long int> wf_time;
    vector<double>   wf_amp;
    vector<double> scint = n.GetS1(quanta,pos_x,pos_y,pos_z,
				   vD,vD_middle,type_num,j,field,keV,useTiming, print_result,
				   wf_time, wf_amp);
    if ( usePD == 0 && fabs(scint[3]) > minS1 && scint[3] < maxS1 )
      signal1.push_back(scint[3]);
    else if ( usePD == 1 && fabs(scint[5]) > minS1 && scint[5] < maxS1 )
      signal1.push_back(scint[5]);
    else if ( usePD >= 2 && fabs(scint[7]) > minS1 && scint[7] < maxS1 )
      signal1.push_back(scint[7]);
    else signal1.push_back(-999.);
    
    if ( pos_z < detector->get_cathode() ) quanta.electrons = 0;
    vector<double> scint2= n.GetS2(quanta.electrons, pos_x, pos_y, driftTime, vD, j, field, useTiming, print_result, wf_time, wf_amp, g2_params);
    if ( usePD == 0 && fabs(scint2[5]) > minS2 && scint2[5] < maxS2 )
      signal2.push_back(scint2[5]);
    else if ( usePD >= 1 && fabs(scint2[7]) > minS2 && scint2[7] < maxS2 )
      signal2.push_back(scint2[7]); //no spike option for S2
    else signal2.push_back(-999.);

    if ( !MCtruthE ) {
      double Nph, Ne, Wq_eV = 1.9896 + ( 20.8 - 1.9896 ) / ( 1. + pow ( rho / 4.0434, 1.4407 ) ); // out-of-sync danger: copied from NEST.cpp
      if ( usePD == 0 )
	Nph= fabs(scint[3]) / (g1*(1.+detector->get_P_dphe()));
      else if ( usePD == 1 ) Nph = fabs(scint[5]) / g1;
      else Nph = fabs(scint[7]) / g1;
      if ( usePD == 0 )
	Ne = fabs(scint2[5]) / (g2*(1.+detector->get_P_dphe()));
      else Ne = fabs(scint2[7]) / g2;
      if ( signal1.back() <= 0. )
	Nph= 0.;
      if ( signal2.back() <= 0. )
	Ne = 0.;
      if ( yields.Lindhard > DBL_MIN && Nph > 0. && Ne > 0. ) {
	keV = ( Nph + Ne ) * Wq_eV * 1e-3 / yields.Lindhard;
	keVee+=( Nph + Ne ) * Wq_eV * 1e-3; //as alternative, use W_DEFAULT in both places, but won't account for density dependence
      }
      else
	keV = 0.;
    }
    if ( (signal1.back() <= 0. || signal2.back() <= 0.) && field >= FIELD_MIN )
      signalE.push_back(0.);
    else
      signalE.push_back(keV);
    
		// Possible outputs from "scint" vector
		// scint[0] = nHits; // MC-true integer hits in same OR different PMTs, NO double phe effect
		// scint[1] = Nphe; // MC-true integer hits WITH double phe effect (Nphe > nHits)
		// scint[2] = pulseArea; // floating real# smeared DAQ pulse areas in phe, NO XYZ correction
		// scint[3] = pulseAreaC; // smeared DAQ pulse areas in phe, WITH XYZ correction
		// scint[4] = Nphd; // same as pulse area, adjusted/corrected *downward* for 2-PE effect (LUX phd units)
		// scint[5] = NphdC; // same as Nphd, but XYZ-corrected
		// scint[6] = spike; // floating real# spike count, NO XYZ correction
		// scint[7] = spikeC; // floating real# spike count, WITH XYZ correction
		// scint[8] = fdetector->get_g1(); // g1 (light collection efficiency in liquid)

		// Possible outputs from "scint2" vector
		// scint2[0] = Nee; // integer number of electrons unabsorbed in liquid then getting extracted
		// scint2[1] = Nph; // raw number of photons produced in the gas gap
		// scint2[2] = nHits; // MC-true integer hits in same OR different PMTs, NO double phe effect
		// scint2[3] = Nphe; // MC-true integer hits WITH double phe effect (Nphe > nHits). S2 has more steps than S1 (e's 1st)
		// 
		// If S2 threshold is set to positive (normal mode)
		// scint2[4] = pulseArea; // floating real# smeared DAQ pulse areas in phe, NO XYZ correction
		// scint2[5] = pulseAreaC; // smeared DAQ pulse areas in phe, WITH XYZ correction
		// scint2[6] = Nphd; // same as pulse area, adjusted/corrected *downward* for 2-PE effect (LUX phd units)
		// scint2[7] = NphdC; // same as Nphd, but XYZ-corrected
		//
		// If S2 threshold is set to negative (switches from S2 -> S2 bottom, NOT literally negative)
		// scint2[4] = S2b; // floating real# smeared pulse areas in phe ONLY including bottom PMTs, NO XYZ correction
		// scint2[5] = S2bc; // floating real# smeared pulse areas in phe ONLY including bottom PMTs, WITH XYZ correction
		// scint2[6] = S2b / (1.+fdetector->get_P_dphe()); // same as S2b, but adjusted for 2-PE effect (LUX phd units)
		// scint2[7] = S2bc / (1.+fdetector->get_P_dphe()); // same as S2bc, but adjusted for 2-PE effect (LUX phd units)
		// scint2[8] = g2; // g2 = ExtEff * SE, light collection efficiency of EL in gas gap (from CalculateG2)
	  
      // this part to be  modified for python wrapper --->  (print_result was set to be 0 in original NEST code)
      if ( print_result ) { //fabs(scint[7]) > DBL_MIN && fabs(scint2[7]) > DBL_MIN ) { //if you want to skip specific below-threshold events, then please comment in this if statement
          // this part to be modified for python wrapper <----

      if ( type_num == Kr83m && eMin == 9.4 && eMax == 9.4 ) printf ( "%.6f\t", yields.DeltaT_Scint );
			printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%d\t%d\t",keV,field,driftTime,pos_x,pos_y,pos_z,quanta.photons,quanta.electrons); //comment this out when below line in
      //printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%lf\t%lf\t",keV,field,driftTime,pos_x,pos_y,pos_z,yields.PhotonYield,yields.ElectronYield); //for when you want means
      if ( pos_z < detector->get_cathode() && print_result ) printf("g-X ");
      if ( keV > 1000. || scint[5] > maxS1 || scint2[7] > maxS2 ||
	   //switch to exponential notation to make output more readable, if energy is too high (>1 MeV)
	   type == "muon" || type == "MIP" || type == "LIP" || type == "mu" || type == "mu-" ) {
	printf("%e\t%e\t%e\t", scint[2], scint[5], scint[7]);
	printf("%li\t%e\t%e\n", (long)scint2[0], scint2[4], scint2[7]);
      }
      else {
	printf("%.6f\t%.6f\t%.6f\t", scint[2], scint[5], scint[7]); //see GetS1 inside of NEST.cpp for full explanation of all 8 scint return vector elements. Sample 3 most common
	printf("%i\t%.6f\t%.6f\n", (int)scint2[0], scint2[4], scint2[7]); //see GetS2 inside of NEST.cpp for full explanation of all 8 scint2 vector elements. Change as you desire
      }
    } //always execute statement, if(1) above, because if is just place-holder in case you want to drop all sub-threshold data
      // this part to be added for python wrapper --->
      nph.push_back(quanta.photons);
      ne.push_back(quanta.electrons);
      s1_n_hits.push_back(scint[0]);
      s1_n_phe.push_back(scint[1]);
      s1_raw_area_phe.push_back(scint[2]);
      s1c_area_phe.push_back(scint[3]);
      s1_raw_area_phd.push_back(scint[4]);
      s1c_area_phd.push_back(scint[5]);
      s1_raw_spike.push_back(scint[6]);
      s1c_spike.push_back(scint[7]);
      s2_ne_extract.push_back(scint2[0]);
      s2_n_ph.push_back(scint2[1]);
      s2_n_hits.push_back(scint2[2]);
      s2_n_phe.push_back(scint2[3]);
      s2_raw_area_phe.push_back(scint2[4]);
      s2c_area_phe.push_back(scint2[5]);
      s2_raw_area_phd.push_back(scint2[6]);
      s2c_area_phd.push_back(scint2[7]);
      deposite_energy.push_back(keV);
      pos_x_mm.push_back(pos_x);
      pos_y_mm.push_back(pos_y);
      pos_z_mm.push_back(pos_z);
      field_V_per_cm.push_back(field);
      // this part to be added for python wrapper <----
      
  }
  
	if (print_result) {
		if ( eMin != eMax ) {
			if ( useS2 == 2 )
				GetBand ( signal2, signal1, false );
			else
				GetBand ( signal1, signal2, false );
			fprintf(stderr,"Bin Center\tBin Actual\tHist Mean\tMean Error\tHist Sigma\t\tEff[%%>thr]\n");
			for ( int j = 0; j < numBins; j++ ) {
				fprintf(stderr,"%lf\t%lf\t%lf\t%lf\t%lf\t\t%lf\n",band[j][0],band[j][1],band[j][2],band[j][4],band[j][3],band[j][5]*100.);
				if ( band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 || band[j][3] <= 0.0 || band[j][4] <= 0.0 || band[j][5] <= 0.0 ||
						 std::isnan(band[j][0]) || std::isnan(band[j][1]) || std::isnan(band[j][2]) || std::isnan(band[j][3]) || std::isnan(band[j][4]) || std::isnan(band[j][5]) )
		{ if ( eMax != -999. ) {
				if( ( (g1*yieldsMax.PhotonYield) < maxS1 || (g2*yieldsMax.ElectronYield) < maxS2 ) && j != 0)
					cerr << "WARNING: Insufficient number of high-energy events to populate highest bins is likely.\n";
				else
					cerr << "WARNING: Insufficient number of low-energy events to populate lowest bins is likely. Increase minS1 and/or minS2.\n";
			}
			eMax = -999.; }
			}
		}
		else {
			GetBand ( signal1, signal2, true );
			GetEnergyRes ( signalE );
			if ( type_num == NR ) {
				fprintf(stderr,"S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc [keVnr]\tEc Res[%%]\tEff[%%>thr]\tEc [keVee]\n");
				keVee /= numEvts;
			}
			else
				fprintf(stderr,"S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc Mean\t\tEc Res[%%]\tEff[%%>thr]\n"); //the C here refers to the combined (S1+S2) energy scale
			for ( int j = 0; j < numBins; j++ ) {
				fprintf(stderr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",band[j][0],band[j][1]/band[j][0]*100.,
					band[j][2],band[j][3]/band[j][2]*100.,energies[0],energies[1]/energies[0]*100.,energies[2]*100.);
				if ( type_num == NR ) fprintf(stderr,"%lf\n",keVee/energies[2]); else fprintf(stderr,"\n");
				if ( (band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 || band[j][3] <= 0.0 ||
				      std::isnan(band[j][0]) || std::isnan(band[j][1]) || std::isnan(band[j][2]) || std::isnan(band[j][3])) && field >= FIELD_MIN ) {
				  if ( numEvts > 1 )
				    cerr << "CAUTION: YOUR S1 and/or S2 MIN and/or MAX may be set to be too restrictive, please check.\n";
				  else cerr << "CAUTION: Poor stats. You must have at least 2 events to calculate S1 and S2 and E resolutions.\n";
				}
				else if ( (energies[0] == eMin || energies[0] == eMax || energies[1] <= 0.0) && field >= FIELD_MIN )
		cerr << "If your energy resolution is 0% then you probably still have MC truth energy on." << endl;
				else ;
			}
		}
	}
  
  return 1;
  
}

vector<vector<double>> GetBand ( vector<double> S1s,
				 vector<double> S2s, bool resol ) {
  
  vector<vector<double>> signals;
  signals.resize(numBins,vector<double>(1,-999.));
  double binWidth, border;
  if ( useS2 == 2 ) {
    binWidth = ( maxS2 - minS2 ) / double(numBins);
    border = minS2;
  }
  else {
    binWidth = ( maxS1 - minS1 ) / double(numBins);
    border = minS1;
  }
  int i = 0, j = 0; double s1c, numPts;
  unsigned long reject[NUMBINS_MAX] = {0};
  
  if ( resol ) {
    numBins = 1;
    binWidth = DBL_MAX;
  }
  
  for ( i = 0; i < S1s.size(); i++ ) {
    for ( j = 0; j < numBins; j++ ) {
      s1c = border + binWidth/2. + double(j) * binWidth;
      if ( i == 0 && !resol ) band[j][0] = s1c;
      if ( fabs(S1s[i]) > (s1c-binWidth/2.) && fabs(S1s[i]) <= (s1c+binWidth/2.) ) {
	if ( S1s[i] >= 0. && S2s[i] >= 0. ) {
	  if ( resol ) {
	    signals[j].push_back(S2s[i]);
	  }
	  else {
	    if ( useS2 == 0 )
	      { if ( S1s[i] && S2s[i] && log10(S2s[i]/S1s[i])>logMin && log10(S2s[i]/S1s[i])<logMax ) signals[j].push_back(log10(S2s[i]/S1s[i])); else signals[j].push_back(0.); }
	    else if ( useS2 == 1 )
	      { if ( S1s[i] && S2s[i] && log10(S2s[i]) > logMin && log10(S2s[i]) < logMax ) signals[j].push_back(log10(S2s[i])); else signals[j].push_back(0.); }
	    else
	      { if ( S1s[i] && S2s[i] && log10(S1s[i]/S2s[i])>logMin && log10(S1s[i]/S2s[i])<logMax ) signals[j].push_back(log10(S1s[i]/S2s[i])); else signals[j].push_back(0.); }
	  }
	  band[j][2] += signals[j].back();
	  if ( resol )
	    band[j][0] += S1s[i];
	  else
	    band[j][1] += S1s[i];
	}
	else
	  reject[j]++;
	break; }
    }
  }
  
  for ( j = 0; j < numBins; j++ ) {
    if ( band[j][0] <= 0. && !resol ) band[j][0] = border + binWidth/2. + double(j) * binWidth;
    signals[j].erase(signals[j].begin());
    numPts = (double)signals[j].size();
    if ( numPts <= 0 && resol ) {
      for ( i = 0; i < S1s.size(); i++ ) band[j][0] += fabs(S1s[i]);
      numPts = S1s.size();
    }
    if (resol)
      band[j][0] /= numPts;
    band[j][1] /= numPts;
    band[j][2] /= numPts;
    for ( i = 0; i <(int)numPts; i++ ) {
      if ( signals[j][i] != -999. ) band[j][3] += pow(signals[j][i]-band[j][2],2.); //std dev calc
    }
    for ( i = 0; i < S1s.size(); i++ ) {
      if ( resol && S1s[i] > 0.0 && S2s[i] > 0.0 ) band[j][1] += pow(S1s[i]-band[j][0],2.); //std dev calc
    }
    band[j][3] /= numPts - 1.;
    band[j][3] = sqrt(band[j][3]);
    if ( resol ) {
      band[j][1] /= numPts - 1.;
      band[j][1] = sqrt(band[j][1]);
    }
    band[j][4] = band[j][3] / sqrt ( numPts );
    band[j][5] = numPts/
      (numPts+double(reject[j]));
  }
  
  return signals;
  
}

void GetEnergyRes ( vector<double> Es ) {
  
  int i, numPts = Es.size();
  double numerator = 0.;
  
  for ( i = 0; i < numPts; i++ ) {
    if ( Es[i] > 0. )
      { energies[0] += Es[i]; numerator++; }
  }
  
  energies[0] /= numerator;
  
  for ( i = 0; i < numPts; i++ ) {
    if ( Es[i] > 0. )
      energies[1] += pow(energies[0]-Es[i],2.);
  }
  
  energies[1] /= numerator - 1.;
  energies[1] = sqrt(energies[1]);
  
  energies[2] = numerator / double ( numPts );
  return;
  
}


// this part to be added for python wrapper --->
using namespace boost::python;
BOOST_PYTHON_MODULE(nest_py_interfaceQym10)
{
    class_<std::vector<boost::python::list> >("double_vector1")
    .def(vector_indexing_suite<std::vector<boost::python::list> >());
    
    class_<std::vector<int> >("double_vector2")
    .def(vector_indexing_suite<std::vector<int> >());
    
    class_<std::vector<double> >("double_vector3")
    .def(vector_indexing_suite<std::vector<double> >());
    
    nestOP (nest_py_wrapper_std::*testNEST1)(boost::python::list) = &nest_py_wrapper_std::testNEST;
    nestOP (nest_py_wrapper_mass::*testNEST2)(int) = &nest_py_wrapper_mass::testNEST;
    
    class_<nest_py_wrapper>("nest_py_wrapper")
    .def("print_result_or_not", &nest_py_wrapper::print_result_or_not)
    ;
    
    class_<nest_py_wrapper_std, bases<nest_py_wrapper>>("nest_py_wrapper_std")
    .def("testNEST", testNEST1)
    ;
    
    class_<nest_py_wrapper_mass, bases<nest_py_wrapper>>("nest_py_wrapper_mass")
    .def("testNEST", testNEST2)
    .def("reset_input",&nest_py_wrapper_mass::reset_input)
    .def("input_type_interaction",&nest_py_wrapper_mass::input_type_interaction)
    .def("input_energy_list",&nest_py_wrapper_mass::input_energy_list)
    .def("input_x_list",&nest_py_wrapper_mass::input_x_list)
    .def("input_y_list",&nest_py_wrapper_mass::input_y_list)
    .def("input_z_list",&nest_py_wrapper_mass::input_z_list)
    .def("input_field_list",&nest_py_wrapper_mass::input_field_list)
    .def("set_seed",&nest_py_wrapper_mass::set_seed)
    ;
    
    class_<nestOP>("nestOP")
    .def_readonly("nph",&nestOP::nph)
    .def_readonly("ne",&nestOP::ne)
    .def_readonly("s1_n_hits",&nestOP::s1_n_hits)
    .def_readonly("s1_n_phe",&nestOP::s1_n_phe)
    .def_readonly("s1_raw_area_phe",&nestOP::s1_raw_area_phe)
    .def_readonly("s1c_area_phe",&nestOP::s1c_area_phe)
    .def_readonly("s1_raw_area_phd",&nestOP::s1_raw_area_phd)
    .def_readonly("s1c_area_phd",&nestOP::s1c_area_phd)
    .def_readonly("s1_raw_spike",&nestOP::s1_raw_spike)
    .def_readonly("s1c_spike",&nestOP::s1c_spike)
    .def_readonly("s2_ne_extract",&nestOP::s2_ne_extract)
    .def_readonly("s2_n_ph",&nestOP::s2_n_ph)
    .def_readonly("s2_n_hits",&nestOP::s2_n_hits)
    .def_readonly("s2_n_phe",&nestOP::s2_n_phe)
    .def_readonly("s2_raw_area_phe",&nestOP::s2_raw_area_phe)
    .def_readonly("s2c_area_phe",&nestOP::s2c_area_phe)
    .def_readonly("s2_raw_area_phd",&nestOP::s2_raw_area_phd)
    .def_readonly("s2c_area_phd",&nestOP::s2c_area_phd)
    .def_readonly("deposite_energy",&nestOP::deposite_energy)
    .def_readonly("pos_x_mm",&nestOP::pos_x_mm)
    .def_readonly("pos_y_mm",&nestOP::pos_y_mm)
    .def_readonly("pos_z_mm",&nestOP::pos_z_mm)
    .def_readonly("field_V_per_cm",&nestOP::field_V_per_cm);
}
// this part to be added for python wrapper <---
