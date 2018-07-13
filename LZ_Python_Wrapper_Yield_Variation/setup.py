from distutils.core import setup, Extension
#from distutils.extension import Extension

nest_py_interface = Extension(
    'nest_py_interface',
    sources=['testNEST.cpp','TestSpectra.cpp','NEST.cpp','VDetector.cpp','RandomGen.cpp'],
                    extra_compile_args = ['-O2', '-std=c++11', '-stdlib=libc++', '-mmacosx-version-min=10.7'],
                    libraries=['boost_python3']#,
)

setup(
    name='nest_py_interface',
    version='0.1',
    ext_modules=[nest_py_interface])
