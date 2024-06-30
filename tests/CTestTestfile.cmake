# CMake generated Testfile for 
# Source directory: /mnt/c/users/korin/documents/uni heidelberg/deal.II/DEAL_Stokes_Praktikum/tests
# Build directory: /mnt/c/users/korin/documents/uni heidelberg/deal.II/DEAL_Stokes_Praktikum/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(tests/template.debug "/usr/bin/cmake" "-DTRGT=tests.template.debug.test" "-DTEST=tests/template.debug" "-DEXPECT=PASSED" "-DBINARY_DIR=/mnt/c/users/korin/documents/uni heidelberg/deal.II/DEAL_Stokes_Praktikum" "-P" "/usr/share/deal.ii//scripts/run_test.cmake")
set_tests_properties(tests/template.debug PROPERTIES  LABEL "tests" PROCESSORS "1" TIMEOUT "600" WORKING_DIRECTORY "/mnt/c/users/korin/documents/uni heidelberg/deal.II/DEAL_Stokes_Praktikum/tests/template.debug" _BACKTRACE_TRIPLES "/usr/share/deal.ii/macros/macro_deal_ii_add_test.cmake;577;add_test;/usr/share/deal.ii/macros/macro_deal_ii_pickup_tests.cmake;343;deal_ii_add_test;/mnt/c/users/korin/documents/uni heidelberg/deal.II/DEAL_Stokes_Praktikum/tests/CMakeLists.txt;2;DEAL_II_PICKUP_TESTS;/mnt/c/users/korin/documents/uni heidelberg/deal.II/DEAL_Stokes_Praktikum/tests/CMakeLists.txt;0;")
