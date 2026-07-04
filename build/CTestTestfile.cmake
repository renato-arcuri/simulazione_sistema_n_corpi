# CMake generated Testfile for 
# Source directory: /home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi
# Build directory: /home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test([=[nbody.t]=] "/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/build/Debug/nbody.t")
  set_tests_properties([=[nbody.t]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/CMakeLists.txt;37;add_test;/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test([=[nbody.t]=] "/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/build/Release/nbody.t")
  set_tests_properties([=[nbody.t]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/CMakeLists.txt;37;add_test;/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test([=[nbody.t]=] "/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/build/RelWithDebInfo/nbody.t")
  set_tests_properties([=[nbody.t]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/CMakeLists.txt;37;add_test;/home/nicolo/Scrivania/Fisica/Programmazione/simulazione_sistema_n_corpi/CMakeLists.txt;0;")
else()
  add_test([=[nbody.t]=] NOT_AVAILABLE)
endif()
