# ViennaRNA.spec
#
# Copyright (c) 1994-2016 Ivo Hofacker, Peter Stadler, ivo@tbi.univie.ac.at
#

%{!?__python2: %global __python2 /usr/bin/python2}
%{!?python2_sitelib: %global python2_sitelib %(%{__python2} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())")}
%{!?python2_sitearch: %global python2_sitearch %(%{__python2} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1))")}
%global __provides_exclude_from ^%{python2_sitearch}/.*\\.so$

##%if 0%{?suse_version} || 0%{?fedora} || 0%{?rhel_version} || 0%{?centos_version}
%{!?_pkgdocdir: %global _pkgdocdir %%{_docdir}/%{name}}
##%else
##%{!?_pkgdocdir: %global _pkgdocdir %{buildroot}%{_docdir}/%{name}-%{version}}
##%endif

Name:           viennarna
Version:        2.2.4
Release:        1%{?dist}
Summary:        RNA Secondary Structure Prediction and Comparison
Provides:       ViennaRNA = %{version}-%{release}

Vendor:         Ivo Hofacker, TBI - University of Vienna
Packager:       Ronny Lorenz <ronny@tbi.univie.ac.at>
Distribution:   viennarna-package

Group:          Applications/Engineering
License:        Free for non commercial use.
URL:            http://www.tbi.univie.ac.at/RNA
Source0:        http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-%{version}.tar.gz
BuildRoot:      %(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)

BuildRequires:  autoconf
BuildRequires:  automake
BuildRequires:  libtool
BuildRequires:  rpm-devel

BuildRequires:  libstdc++-devel
BuildRequires:  gcc gcc-c++ glibc-devel info
BuildRequires:  gsl-devel 
BuildRequires:  perl(ExtUtils::Embed)

Requires:       perl libstdc++ glibc info gsl
Provides:       Kinfold

%if 0%{?suse_version}
BuildRequires:  python-devel 
%else
BuildRequires:  python2-devel 
%endif


%description
The ViennaRNA package consists of a library and several standalone
programs for RNA secondary structure analysis. It includes algorithms
for predicting optimal and suboptimal secondary structures, base pair
probabilities and partition functions, for comparing secondary
structures, and the design of RNA sequences with a desired structure.

%package devel
Summary:  Library and header files for ViennaRNA RNAlib
Group:    Development/Libraries
Provides: libRNA.a = %{version}-%{release}
Requires: %{name} = %{version}-%{release}
Requires: libstdc++-devel
Requires: pkgconfig

%description devel
Header files for ViennaRNA.

%package -n perl-rna
Summary:  Perl binding for ViennaRNA RNAlib
Group:    Development/Languages/Perl
Requires: %{name} = %{version}-%{release}
Requires: perl


%description -n perl-rna
Perl binding for ViennaRNA RNAlib.

%package -n python2-rna
Summary:  Python binding for ViennaRNA RNAlib
Group:    Libraries/Python
Requires: %{name} = %{version}-%{release}
Requires: python

%description -n python2-rna
Python binding for ViennaRNA RNAlib.

%package compat
Summary:    Transitional package for ViennaRNA to new package split
Group:      Applications/Engineering
Obsoletes:  ViennaRNA < 2.2.1
Requires:   viennarna = %{version}
Requires:   viennarna-devel = %{version}
Requires:   perl-rna = %{version}
Requires:   python2-rna = %{version}

%description compat
This package only exists to help transition ViennaRNA users to the new
package split. Older versions of ViennaRNA contained the entire package.
With version 2.2.1, we split the package into core, development, and
scripting language interface packages to conform to distribution requirements.
This compat package will be removed after one distribution release cycle,
please do not reference it or depend on it in any way.

%prep
%setup -n ViennaRNA-%{version} -q

%build
%if 0%{?fedora} && 0%{?fedora_version} < 22
%configure --with-cluster --without-forester --disable-lto --docdir=%{_pkgdocdir} INSTALLDIRS=vendor PERL=/usr/bin/perl
%else
%configure --with-cluster --without-forester --docdir=%{_pkgdocdir} INSTALLDIRS=vendor PERL=/usr/bin/perl
%endif
make %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT

%{__make} install DESTDIR=$RPM_BUILD_ROOT

# crude bugfix for info clash
rm -f $RPM_BUILD_ROOT/usr/share/info/dir

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
%{_bindir}/*
%{_datadir}/ViennaRNA/
%{_mandir}/man1/RNAplfold.1.gz
%{_mandir}/man1/AnalyseSeqs.1.gz
%{_mandir}/man1/Kinfold.1.gz
%{_mandir}/man1/RNApdist.1.gz
%{_mandir}/man1/RNAduplex.1.gz
%{_mandir}/man1/RNAup.1.gz
%{_mandir}/man1/RNAheat.1.gz
%{_mandir}/man1/RNAeval.1.gz
%{_mandir}/man1/RNALalifold.1.gz
%{_mandir}/man1/RNAaliduplex.1.gz
%{_mandir}/man1/RNAsubopt.1.gz
%{_mandir}/man1/RNALfold.1.gz
%{_mandir}/man1/RNAfold.1.gz
%{_mandir}/man1/RNApvmin.1.gz
%{_mandir}/man1/AnalyseDists.1.gz
%{_mandir}/man1/RNAsnoop.1.gz
%{_mandir}/man1/RNAplot.1.gz
%{_mandir}/man1/RNAalifold.1.gz
%{_mandir}/man1/RNAinverse.1.gz
%{_mandir}/man1/RNAPKplex.1.gz
%{_mandir}/man1/RNAparconv.1.gz
%{_mandir}/man1/RNAcofold.1.gz
%{_mandir}/man1/RNA2Dfold.1.gz
%{_mandir}/man1/RNApaln.1.gz
%{_mandir}/man1/RNAdistance.1.gz
%{_mandir}/man1/RNAplex.1.gz
%{_mandir}/man1/ct2db.1.gz

%files devel
%defattr(-,root,root,-)
%{_libdir}/libRNA.a
%{_libdir}/libRNA.la
%{_includedir}/ViennaRNA
%{_libdir}/pkgconfig/RNAlib2.pc
%doc %{_pkgdocdir}
%doc /usr/share/info/RNAlib.info.gz

%files -n perl-rna
%defattr(644,root,root,755)
%{perl_vendorlib}/RNA.pm
%{perl_vendorlib}/RNA
%{perl_vendorarch}/auto/RNA
%attr(755,root,root) %{perl_vendorarch}/auto/RNA/RNA.so

%files -n python2-rna
%defattr(644,root,root,755)
%{python2_sitearch}/RNA
%attr(755,root,root) %{python2_sitearch}/RNA/_RNA.so
%{python2_sitelib}/RNA

%files compat
