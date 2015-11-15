#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @SEMBLE_CXXFLAGS@ @JACKFITTER_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

AM_LDFLAGS  = -L@top_builddir@/lib/ \
               @LDFLAGS@ @SEMBLE_LDFLAGS@ @JACKFITTER_LDFLAGS@ @ADAT_LDFLAGS@ 

LDADD = -lradmat @SEMBLE_LIBS@ @JACKFITTER_LIBS@ @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ 

bin_PROGRAMS =   estimate_error\
								 grab_ff		\
								 grab_pol		\
								 grab_vector_mean_cmp		\
								 invert_subduced		\
								 mass_overlap_db_combine		\
								 mass_overlap_db_create		\
								 mass_overlap_db_util		\
								 pole_enhancement			\
								 radmat		\
								 radmat_util		\
								 simple_test		\
								 single_plot		\
								 subduce 				\
								 test_levi_civita \
								 test_redstar_ops		\
								 test_rotations \
								 test_tens		


estimate_error_SOURCES = estimate_error.cc
grab_ff_SOURCES = grab_ff.cc
grab_pol_SOURCES = grab_pol.cc
grab_vector_mean_cmp_SOURCES = grab_vector_mean_cmp.cc
invert_subduced_SOURCES = invert_subduced.cc
mass_overlap_db_combine_SOURCES = mass_overlap_db_combine.cc
mass_overlap_db_create_SOURCES = mass_overlap_db_create.cc
mass_overlap_db_util_SOURCES = mass_overlap_db_util.cc
pole_enhancement_SOURCES = pole_enhancement.cc
radmat_SOURCES = radmat.cc
radmat_util_SOURCES = radmat_util.cc
simple_test_SOURCES = simple_test.cc
single_plot_SOURCES = single_plot.cc
subduce_SOURCES = subduce.cc
test_levi_civita_SOURCES = test_levi_civita.cc
test_redstar_ops_SOURCES = test_redstar_ops.cc
test_rotations_SOURCES = test_rotations.cc
test_tens_SOURCES = test_tens.cc
#
# The folow rule will always try to rebuild libadat.a when you 
# try to compile a main prog
#
DEPENDENCIES = build_lib
${bin_PROGRAMS}: ${DEPENDENCIES}  

#${check_PROGRAMS}: ${DEPENDENCIES}
#${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
