exclude_from_headercheck(
  typetreetestswitch.hh
  typetreetestutility.hh)

dune_add_test(SOURCES testtypetree.cc
              COMPILE_DEFINITIONS TEST_TYPETREE)

dune_add_test(SOURCES testtreecontainer.cc)
