AM_CFLAGS = $(OPENMP_CFLAGS) @AM_CFLAGS@
AM_CXXFLAGS = $(OPENMP_CFLAGS) @AM_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@ @FLOAT_PF_FLAG@ @DEPRECATION_WARNING@

bin_PROGRAMS = RNAfold mmfold

LDADD = $(OPENMP_CFLAGS) $(top_builddir)/src/ViennaRNA/libRNA_conv.la


GENGETOPT_CMDL =  RNAfold_cmdl.c RNAfold_cmdl.h mmfold_cmdl.c mmfold_cmdl.h

GENGETOPT_FILES =  RNAfold.ggo  mmfold.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}


# how we make the commandline interface parser
%_cmdl.h %_cmdl.c: %.ggo
	$(AM_V_GEN)gengetopt -i $<

# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

RNAfold_SOURCES=RNAfold_cmdl.c RNAfold.c

mmfold_SOURCES=mmfold_cmdl.c mmfold.c

