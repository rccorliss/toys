AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(includedir) -I$(OFFLINE_MAIN)/include -I/usr/local/Cellar/root/6.18.04_1/include/root/

#-I$(ROOTSYS)/include 

lib_LTLIBRARIES = libfieldsim.la   

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

#libfieldsim_la_LIBADD = \
#-lfun4all \
#  -lphool \
#  -lSubsysReco

# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOTDICTS = \
  FieldSim_Dict.cc 
# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# MAKEROOT6 is set in the configure.ac, true for root6
if MAKEROOT6
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  FieldSim_rdict.pcm 
else
  ROOT5_DICTS = \
    FieldSim_Dict.cc 
endif

libfieldsim_la_SOURCES = \
  $(ROOTDICTS) \
  $(ROOT5_DICTS) \
  FieldSim.cc 

pkginclude_HEADERS = \
  FieldSim.h 

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libfieldsim.la

testexternals.cc:
	echo "this is a test.  $(ROOTSYS)/include"
	echo $(ROOTSYS)
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm 

