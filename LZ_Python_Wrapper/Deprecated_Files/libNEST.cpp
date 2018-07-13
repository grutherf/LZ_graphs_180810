#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "NEST.h"
using namespace boost::python;

template<class T> boost::python::list std_vector_to_py_list(const std::vector<T>& v)
{
    boost::python::object get_iter = boost::python::iterator<std::vector<T> >();
    boost::python::object iter = get_iter(v);
    boost::python::list l(iter);
    return l;
}

template< class T >
inline
std::vector< T > to_std_vector( const boost::python::object& iterable )
{
    return std::vector< T >( boost::python::stl_input_iterator< T >( iterable ),
                             boost::python::stl_input_iterator< T >( ) );
}


BOOST_PYTHON_MODULE(libNEST){
    
    //! python access to stl integer vectors
    class_< std::vector<float> >("vectorFloat")
        .def(vector_indexing_suite<std::vector<float> >())
    ;

    class_< std::vector<int> >("vectorInt")
        .def(vector_indexing_suite<std::vector<int> >())
    ;


    def("std_vector_to_py_list_float",std_vector_to_py_list<float>);
    def("std_vector_to_py_list_int",std_vector_to_py_list<int>);
    def("to_std_vector",to_std_vector<float>);

    // This will enable user-defined docstrings and python signatures,
    // while disabling the C++ signatures
    docstring_options local_docstring_options(true, true, false);



    class_<NEST>("NEST", init<int,float,float,float,float>())
        .def("SetParticleType",&NEST::SetParticleType)
        .def("GetParticleType",&NEST::GetParticleType)
        .def("SetEnergy",&NEST::SetEnergy)
        .def("GetEnergy",&NEST::GetEnergy)
        .def("SetElectricField",&NEST::SetElectricField)
        .def("GetElectricField",&NEST::GetElectricField)
        .def("SetXeDensity",&NEST::SetXeDensity)
        .def("GetXeDensity",&NEST::GetXeDensity)
        .def("GetNumExcitons",&NEST::GetNumExcitons)
        .def("GetNumIons",&NEST::GetNumIons)
        .def("GetNumQuanta",&NEST::GetNumQuanta)
        .def("GetNumPhotons",&NEST::GetNumPhotons)
        .def("GetNumElectrons",&NEST::GetNumElectrons)
        .def("GetS1",&NEST::GetS1)
        .def("GetS1c",&NEST::GetS1c)
        .def("GetS2",&NEST::GetS2)
        .def("GetS2c",&NEST::GetS2c)
        .def("SetDetectorParameters",&NEST::SetDetectorParameters)
        .def("GetDetectorParameters",&NEST::GetDetectorParameters)
        .def("SetDriftLocation",&NEST::SetDriftLocation)
        .def("GetDriftLocation",&NEST::GetDriftLocation)
        .def("SetRandomSeed",&NEST::SetRandomSeed)
        .def("GetRandomSeed",&NEST::GetRandomSeed)
        .def("ExcitonsIons",&NEST::ExcitonsIons)
        .def("PhotonsElectrons",&NEST::PhotonsElectrons)
        .def("DetectorResponse",&NEST::DetectorResponse)
        .def("SkinResponse",&NEST::SkinResponse)
        .def("GetSkinPDE",&NEST::GetSkinPDE)
        .def("GetSkinEField",&NEST::GetSkinEField)

        
    ;

    class_<Detector>("Detector")
        .def("LUXSettings",&Detector::LUXSettings)
        .def("LZSettings",&Detector::LZSettings)
        .def("LZProjectedSettings",&Detector::LZProjectedSettings)
        .def_readwrite("name",&Detector::name)
        .def_readwrite("cathodeVoltage",&Detector::cathodeVoltage)
        .def_readwrite("g1",&Detector::g1)
        .def_readwrite("sp",&Detector::sp)
        .def_readwrite("ee",&Detector::ee)
        .def_readwrite("g1_gas",&Detector::g1_gas)
        .def_readwrite("ff",&Detector::ff)
        .def_readwrite("s1polA",&Detector::s1polA)
        .def_readwrite("s1polB",&Detector::s1polB)
        .def_readwrite("s1polC",&Detector::s1polC)
        .def_readwrite("s1polD",&Detector::s1polD)
        .def_readwrite("s1polE",&Detector::s1polE)
        .def_readwrite("eDriftSpeed",&Detector::eDriftSpeed)
        .def_readwrite("liquidBorder",&Detector::liquidBorder)
        .def_readwrite("efpolA",&Detector::efpolA)
        .def_readwrite("efpolB",&Detector::efpolB)
        .def_readwrite("efpolC",&Detector::efpolC)
        .def_readwrite("efpolD",&Detector::efpolD)
        .def_readwrite("efpolE",&Detector::efpolE)
        .def_readwrite("efpolF",&Detector::efpolF)
        .def_readwrite("af",&Detector::af)
        .def_readwrite("nFold",&Detector::nFold)
        .def_readwrite("spEff",&Detector::spEff)
        .def_readwrite("spEff_direct",&Detector::spEff_direct)
        .def_readwrite("dtCntr",&Detector::dtCntr)
        .def_readwrite("dt_min",&Detector::dt_min)
        .def_readwrite("dt_max",&Detector::dt_max)
        .def_readwrite("e_life",&Detector::e_life)
        .def_readwrite("s2_thr",&Detector::s2_thr)
        .def_readwrite("s2_rat",&Detector::s2_rat)
        .def_readwrite("gasRho",&Detector::gasRho)
        .def_readwrite("gasGap",&Detector::gasGap)
        .def_readwrite("elGain",&Detector::elGain)
        .def_readwrite("elThld",&Detector::elThld)
        .def_readwrite("numberPMTs1",&Detector::numberPMTs)
        .def_readwrite("erf_mu",&Detector::erf_mu)
        .def_readwrite("erf_sg",&Detector::erf_sg)
        .def_readwrite("Xe_rho",&Detector::Xe_rho)
        .def_readwrite("P_dphe",&Detector::P_dphe)
        .def_readwrite("logMax",&Detector::logMax)
        .def_readwrite("logMin",&Detector::logMin)
        .def_readwrite("vuv1_1to1_2",&Detector::vuv1_1to1_2)
        .def_readwrite("coinWind",&Detector::coinWind)
        .def_readwrite("s1_min",&Detector::s1_min)
        .def_readwrite("s1_max",&Detector::s1_max)
    ;
    
}
