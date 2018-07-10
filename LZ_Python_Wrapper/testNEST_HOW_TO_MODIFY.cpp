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
 boost::python::list s1_nhits;
 boost::python::list s1_nphe;
 boost::python::list s1_area_raw_phe;
 boost::python::list s1c_area_phe;
 boost::python::list s1_area_raw_phd;
 boost::python::list s1c_area_phd;
 boost::python::list s1_raw_spike;
 boost::python::list s1c_spike;
 boost::python::list s2_ne_extract;
 boost::python::list s2_nph;
 boost::python::list s2_nhits;
 boost::python::list s2_nphe;
 boost::python::list s2_area_raw_phe;
 boost::python::list s2c_area_phe;
 boost::python::list s2_area_raw_phd;
 boost::python::list s2c_area_phd;
 */

...............

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

...............

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
    boost::python::list s1_nhits;
    boost::python::list s1_nphe;
    boost::python::list s1_area_raw_phe;
    boost::python::list s1c_area_phe;
    boost::python::list s1_area_raw_phd;
    boost::python::list s1c_area_phd;
    boost::python::list s1_raw_spike;
    boost::python::list s1c_spike;
    boost::python::list s2_ne_extract;
    boost::python::list s2_nph;
    boost::python::list s2_nhits;
    boost::python::list s2_nphe;
    boost::python::list s2_area_raw_phe;
    boost::python::list s2c_area_phe;
    boost::python::list s2_area_raw_phd;
    boost::python::list s2c_area_phd;
};

class nest_py_wrapper{
protected:
    int argc;
    char** argv;
    nestOP nest_op;
    int print_result = 1;
    int testNEST(int argc, char** argv, int print_result,
                 vector<int>& nph, vector<int>& ne, vector<int>& s1_nhits, vector<int>& s1_nphe,
                 vector<double>& s1_area_raw_phe, vector<double>& s1c_area_phe, vector<double>& s1_area_raw_phd, vector<double>& s1c_area_phd, vector<double>& s1_raw_spike, vector<double>& s1c_spike,
                 vector<int>& s2_ne_extract, vector<double>& s2_nph, vector<double>& s2_nhits, vector<double>& s2_nphe,
                 vector<double>& s2_area_raw_phe, vector<double>& s2c_area_phe, vector<double>& s2_area_raw_phd, vector<double>& s2c_area_phd);
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
    vector<int> s1_nhits;
    vector<int> s1_nphe;
    vector<double> s1_area_raw_phe;
    vector<double> s1c_area_phe;
    vector<double> s1_area_raw_phd;
    vector<double> s1c_area_phd;
    vector<double> s1_raw_spike;
    vector<double> s1c_spike;
    vector<int> s2_ne_extract;
    vector<double> s2_nph;
    vector<double> s2_nhits;
    vector<double> s2_nphe;
    vector<double> s2_area_raw_phe;
    vector<double> s2c_area_phe;
    vector<double> s2_area_raw_phd;
    vector<double> s2c_area_phd;
    
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
        nest_py_wrapper::testNEST(argc, argv, print_result, nph, ne, s1_nhits, s1_nphe,
                                  s1_area_raw_phe, s1c_area_phe, s1_area_raw_phd, s1c_area_phd, s1_raw_spike, s1c_spike,
                                  s2_ne_extract, s2_nph, s2_nhits, s2_nphe,
                                  s2_area_raw_phe, s2c_area_phe, s2_area_raw_phd, s2c_area_phd);
        
    }
    free(argv);
    nest_op.nph = (std_vector_to_py_list<int>(nph));
    nest_op.ne = (std_vector_to_py_list<int>(ne));
    nest_op.s1_nhits = (std_vector_to_py_list<int>(s1_nhits));
    nest_op.s1_nphe = (std_vector_to_py_list<int>(s1_nphe));
    nest_op.s1_area_raw_phe = (std_vector_to_py_list<double>(s1_area_raw_phe));
    nest_op.s1c_area_phe = (std_vector_to_py_list<double>(s1c_area_phe));
    nest_op.s1_area_raw_phd = (std_vector_to_py_list<double>(s1_area_raw_phd));
    nest_op.s1c_area_phd = (std_vector_to_py_list<double>(s1c_area_phd));
    nest_op.s1_raw_spike = (std_vector_to_py_list<double>(s1_raw_spike));
    nest_op.s1c_spike = (std_vector_to_py_list<double>(s1c_spike));
    nest_op.s2_ne_extract = (std_vector_to_py_list<int>(s2_ne_extract));
    nest_op.s2_nph = (std_vector_to_py_list<double>(s2_nph));
    nest_op.s2_nhits = (std_vector_to_py_list<double>(s2_nhits));
    nest_op.s2_nphe = (std_vector_to_py_list<double>(s2_nphe));
    nest_op.s2_area_raw_phe = (std_vector_to_py_list<double>(s2_area_raw_phe));
    nest_op.s2c_area_phe = (std_vector_to_py_list<double>(s2c_area_phe));
    nest_op.s2_area_raw_phd = (std_vector_to_py_list<double>(s2_area_raw_phd));
    nest_op.s2c_area_phd = (std_vector_to_py_list<double>(s2c_area_phd));
    
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
    vector<int> s1_nhits;
    vector<int> s1_nphe;
    vector<double> s1_area_raw_phe;
    vector<double> s1c_area_phe;
    vector<double> s1_area_raw_phd;
    vector<double> s1c_area_phd;
    vector<double> s1_raw_spike;
    vector<double> s1c_spike;
    vector<int> s2_ne_extract;
    vector<double> s2_nph;
    vector<double> s2_nhits;
    vector<double> s2_nphe;
    vector<double> s2_area_raw_phe;
    vector<double> s2c_area_phe;
    vector<double> s2_area_raw_phd;
    vector<double> s2c_area_phd;
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
    nest_py_wrapper::testNEST(argc, argv, print_result, nph, ne, s1_nhits, s1_nphe,
                              s1_area_raw_phe, s1c_area_phe, s1_area_raw_phd, s1c_area_phd, s1_raw_spike, s1c_spike,
                              s2_ne_extract, s2_nph, s2_nhits, s2_nphe,
                              s2_area_raw_phe, s2c_area_phe, s2_area_raw_phd, s2c_area_phd);
    
    
    free(argv);
    
    // define the print function
    boost::python::object print =
    boost::python::import("__main__").attr("__builtins__").attr("print");
    
    nest_op.nph = (std_vector_to_py_list<int>(nph));
    nest_op.ne = (std_vector_to_py_list<int>(ne));
    nest_op.s1_nhits = (std_vector_to_py_list<int>(s1_nhits));
    nest_op.s1_nphe = (std_vector_to_py_list<int>(s1_nphe));
    nest_op.s1_area_raw_phe = (std_vector_to_py_list<double>(s1_area_raw_phe));
    nest_op.s1c_area_phe = (std_vector_to_py_list<double>(s1c_area_phe));
    nest_op.s1_area_raw_phd = (std_vector_to_py_list<double>(s1_area_raw_phd));
    nest_op.s1c_area_phd = (std_vector_to_py_list<double>(s1c_area_phd));
    nest_op.s1_raw_spike = (std_vector_to_py_list<double>(s1_raw_spike));
    nest_op.s1c_spike = (std_vector_to_py_list<double>(s1c_spike));
    nest_op.s2_ne_extract = (std_vector_to_py_list<int>(s2_ne_extract));
    nest_op.s2_nph = (std_vector_to_py_list<double>(s2_nph));
    nest_op.s2_nhits = (std_vector_to_py_list<double>(s2_nhits));
    nest_op.s2_nphe = (std_vector_to_py_list<double>(s2_nphe));
    nest_op.s2_area_raw_phe = (std_vector_to_py_list<double>(s2_area_raw_phe));
    nest_op.s2c_area_phe = (std_vector_to_py_list<double>(s2c_area_phe));
    nest_op.s2_area_raw_phd = (std_vector_to_py_list<double>(s2_area_raw_phd));
    nest_op.s2c_area_phd = (std_vector_to_py_list<double>(s2c_area_phd));
    
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

//This is what used to be the method declaration of the main method. Simply change the name and arguments
// this part to be modified for python wrapper --->
int nest_py_wrapper::testNEST(int argc, char** argv, int print_result,
                              vector<int>& nph, vector<int>& ne, vector<int>& s1_nhits, vector<int>& s1_nphe,
                              vector<double>& s1_area_raw_phe, vector<double>& s1c_area_phe, vector<double>& s1_area_raw_phd, vector<double>& s1c_area_phd, vector<double>& s1_raw_spike, vector<double>& s1c_spike,
                              vector<int>& s2_ne_extract, vector<double>& s2_nph, vector<double>& s2_nhits, vector<double>& s2_nphe,
                              vector<double>& s2_area_raw_phe, vector<double>& s2c_area_phe, vector<double>& s2_area_raw_phd, vector<double>& s2c_area_phd){
    // this part to be modified for python wrapper <---
    
    ...............
  
    double rho = n.SetDensity ( detector->get_T_Kelvin(), detector->get_p_bar() ); //cout.precision(12);
  if ( rho < 1. ) detector->set_inGas(true);
  
	// Print g1 and g2 values, requiring S2 calculation for 1 SE and x,y,d = 0,0,0
    // this part to be modified for python wrapper --->
    if(print_result){
        // this part to be modified for python wrapper <---
        vector<double> secondary = n.GetS2 ( 1, 0., 0., 0., 1., 0, 1e2, false );
        //to be added --->
    }
    // to be added <---
	if ( atof(argv[5]) == -1. ) {
    vTable = n.SetDriftVelocity_NonUniform(rho, z_step);
    vD_middle = vTable[int(floor(.5*detector->get_TopDrift()/z_step))];
  }
    
  ...............
    
      // this part to be added for python wrapper --->  (print_result was set to be 0 in original NEST code)
      if ( print_result ) { //fabs(scint[7]) > DBL_MIN && fabs(scint2[7]) > DBL_MIN ) { //if you want to skip specific below-threshold events, then please comment in this if statement
          // this part to be added for python wrapper <----
      printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%d\t%d\t",keV,field,driftTime,pos_x,pos_y,pos_z,quanta.photons,quanta.electrons); //comment this out when below line in
      //printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%lf\t%lf\t",keV,field,driftTime,pos_x,pos_y,pos_z,yields.PhotonYield,yields.ElectronYield); //for when you want means
      if ( pos_z < detector->get_cathode() ) printf("g-X ");
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
      s1_nhits.push_back(scint[0]);
      s1_nphe.push_back(scint[1]);
      s1_area_raw_phe.push_back(scint[2]);
      s1c_area_phe.push_back(scint[3]);
      s1_area_raw_phd.push_back(scint[4]);
      s1c_area_phd.push_back(scint[5]);
      s1_raw_spike.push_back(scint[6]);
      s1c_spike.push_back(scint[7]);
      s2_ne_extract.push_back(scint2[0]);
      s2_nph.push_back(scint2[1]);
      s2_nhits.push_back(scint2[2]);
      s2_nphe.push_back(scint2[3]);
      s2_area_raw_phe.push_back(scint2[4]);
      s2c_area_phe.push_back(scint2[5]);
      s2_area_raw_phd.push_back(scint2[6]);
      s2c_area_phd.push_back(scint2[7]);
      // this part to be added for python wrapper <----
  }
    // this part to be added for python wrapper --->
    if(print_result){
        // this part to be added for python wrapper <---
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
	    if( ( (detector->get_g1()*yieldsMax.PhotonYield) < maxS1 || (g2*yieldsMax.ElectronYield) < maxS2 ) && j != 0)
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
      if ( band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 || band[j][3] <= 0.0 ||
	   std::isnan(band[j][0]) || std::isnan(band[j][1]) || std::isnan(band[j][2]) || std::isnan(band[j][3]) )
	cerr << "CAUTION: YOUR S1 and/or S2 MIN and/or MAX may be set to be too restrictive, please check.\n";
      else if ( energies[0] == eMin || energies[0] == eMax || energies[1] <= 0.0 )
	cerr << "If your energy resolution is 0% then you probably still have MC truth energy on." << endl;
      else ;
    }
  }
        //to be added ----->
    }
    // to be added <-----
  return 1;
  
}

...............
            
// this part to be added for python wrapper --->
using namespace boost::python;
BOOST_PYTHON_MODULE(nest_py_interface)
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
    .def_readonly("s1_nhits",&nestOP::s1_nhits)
    .def_readonly("s1_nphe",&nestOP::s1_nphe)
    .def_readonly("s1_area_raw_phe",&nestOP::s1_area_raw_phe)
    .def_readonly("s1c_area_phe",&nestOP::s1c_area_phe)
    .def_readonly("s1_area_raw_phd",&nestOP::s1_area_raw_phd)
    .def_readonly("s1c_area_phd",&nestOP::s1c_area_phd)
    .def_readonly("s1_raw_spike",&nestOP::s1_raw_spike)
    .def_readonly("s1c_spike",&nestOP::s1c_spike)
    .def_readonly("s2_ne_extract",&nestOP::s2_ne_extract)
    .def_readonly("s2_nph",&nestOP::s2_nph)
    .def_readonly("s2_nhits",&nestOP::s2_nhits)
    .def_readonly("s2_nphe",&nestOP::s2_nphe)
    .def_readonly("s2_area_raw_phe",&nestOP::s2_area_raw_phe)
    .def_readonly("s2c_area_phe",&nestOP::s2c_area_phe)
    .def_readonly("s2_area_raw_phd",&nestOP::s2_area_raw_phd)
    .def_readonly("s2c_area_phd",&nestOP::s2c_area_phd);
}
// this part to be added for python wrapper <---
