AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = AnalyseDists AnalyseSeqs
man_MANS = AnalyseDists.1 AnalyseSeqs.1

AM_CFLAGS = $(OPENMP_CFLAGS) @AM_CFLAGS@
AM_CXXFLAGS = $(OPENMP_CFLAGS) @AM_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir) -I$(top_builddir) @AM_CPPFLAGS@ @FLOAT_PF_FLAG@ @DEPRECATION_WARNING@
AM_LDFLAGS = @AM_LDFLAGS@

LDADD = $(OPENMP_CFLAGS) $(top_builddir)/src/ViennaRNA/libRNA_conv.la

AnalyseDists_SOURCES = split.c \
          distance_matrix.c \
          cluster.c \
          treeplot.c \
          AD_main.c

AnalyseSeqs_SOURCES = statgeom.c \
          distance_matrix.c \
          cluster.c \
          AS_main.c \
          split.c \
          treeplot.c \
          PS3D.c

noinst_HEADERS = cluster.h split.h treeplot.h distance_matrix.h statgeom.h \
                PS3D.h StrEdit_CostMatrix.h  

EXTRA_DIST = $(man_MANS)
