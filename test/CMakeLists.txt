
add_pfunit_ctest (numerov_test
TEST_SOURCES test_integrate.pf
LINK_LIBRARIES sut
MAX_PES 1
)

add_pfunit_ctest(bessel_test
TEST_SOURCES test_special_functions.pf
LINK_LIBRARIES sut
MAX_PES 1 # what does this do again?
)
