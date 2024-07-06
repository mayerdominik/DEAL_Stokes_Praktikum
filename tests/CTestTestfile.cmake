# CMake generated Testfile for 
# Source directory: /home/dome/projects/DEAL_Stokes_Praktikum/tests
# Build directory: /home/dome/projects/DEAL_Stokes_Praktikum/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(tests/template.debug "/usr/bin/cmake" "-DTRGT=tests.template.debug.test" "-DTEST=tests/template.debug" "-DEXPECT=PASSED" "-DBINARY_DIR=/home/dome/projects/DEAL_Stokes_Praktikum" "-P" "/usr/local/share/deal.II/scripts/run_test.cmake")
set_tests_properties(tests/template.debug PROPERTIES  LABEL "tests" PROCESSORS "1" TIMEOUT "600" WORKING_DIRECTORY "/home/dome/projects/DEAL_Stokes_Praktikum/tests/template.debug" _BACKTRACE_TRIPLES "/usr/local/share/deal.II/macros/macro_deal_ii_add_test.cmake;577;add_test;/usr/local/share/deal.II/macros/macro_deal_ii_pickup_tests.cmake;343;deal_ii_add_test;/home/dome/projects/DEAL_Stokes_Praktikum/tests/CMakeLists.txt;2;DEAL_II_PICKUP_TESTS;/home/dome/projects/DEAL_Stokes_Praktikum/tests/CMakeLists.txt;0;")
