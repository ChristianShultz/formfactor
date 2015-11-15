#!/bin/bash

a='
#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
								@CXXFLAGS@ @SEMBLE_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ @JACKFITTER_CXXFLAGS@

#utils
nobase_include_HEADERS =   \'

b=$(./regen_makefile_helper.csh 1)

c='
# the lib
lib_LIBRARIES = libradmat.a


libradmat_a_SOURCES =	\'

d=$(./regen_makefile_helper.csh 2)

f="Makefile.am"

printf '%s\n' "$a" > $f
printf '%s\n' "$b" >> $f
printf '%s\n' "$c" >> $f
printf '%s\n' "$d" >> $f

