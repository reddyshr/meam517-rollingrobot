#-------------------------------------------------------------------------------
AC_PREREQ(2.69)
AC_INIT([SNOPT],[7.6],[optimizers@ccom.ucsd.edu],[snopt])

AC_CONFIG_SRCDIR([src/sn02lib.f])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config/m4])

VERSION="7.6.0";
AC_SUBST([VERSION])

AC_PREFIX_DEFAULT([`pwd`])

user_FCFLAGS=${FCFLAGS}
user_CFLAGS=${CFLAGS}
user_CPPFLAGS=${CPPFLAGS}
user_CXXFLAGS=${CXXFLAGS}
user_LDFLAGS=${LDFLAGS}


# Libtool
LT_INIT
AC_SUBST([LIBTOOL_DEPS])


# Check compilers.
AC_PROG_INSTALL
AC_PROG_FC([gfortran nagfor ifort])
AC_PROG_CC
AC_PROG_CXX

AC_FC_LIBRARY_LDFLAGS
AC_SUBST([FCLIBS])

AC_FC_MODULE_OUTPUT_FLAG([],
  [ case $FC in
      *nagfor* ) FC_MODOUT="-mdir ";;
    esac
  ])
AC_FC_MODULE_FLAG([],
  [ case $FC in
      *nagfor* ) FC_MODINC="-mdir ";;
    esac
  ])
AC_SUBST(FC_MODOUT)
AC_SUBST(FC_MODINC)

AC_PROG_MAKE_SET
AC_SUBST([SET_MAKE])


dnl ----------------------------------------------------------------------------

# BLAS library.
AC_ARG_WITH([blas],
	    [AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])],
	    [with_blas=${withval} ], [ blas=no ] )
case $with_blas in
	yes | "" | no) blas=no ;;
	-* | */* | *.a | *.so | *.so.* | *.o)
	  blas=yes; BLAS_LIBS="$with_blas" ;;
        matlab )
	  blas=yes; BLAS_LIBS="matlab" ;;
	*)
	  blas=yes; BLAS_LIBS="-l$with_blas" ;;
esac
AC_SUBST([blas])
AC_SUBST([BLAS_LIBS])

dnl ----------------------------------------------------------------------------

AC_ARG_WITH([matlab],
	    [AS_HELP_STRING([--with-matlab=DIR],[Provide location of MATLAB]) ],
            [with_matlab=yes; MATLABDIR=${withval};],
            [with_matlab=no; matlab="no"; MATLABDIR="";] )
MEX="";
if test x"${with_matlab}" == xyes; then
  AC_MSG_CHECKING([for Matlab])

  if test -e ${MATLABDIR}/bin/mex; then
    with_matlab=yes;
    matlab=yes;
    MEX="${MATLABDIR}/bin/mex";
    AC_MSG_RESULT([yes])
    AC_MSG_WARN([Please make sure to use a compiler supported by Matlab]);
  else
    with_matlab=no;
    matlab=no;
    MEX="";
    AC_MSG_RESULT([mex not found in ${MATLABDIR}/bin])
  fi
fi
AC_SUBST(matlab)
AC_SUBST(MATLABDIR)
AC_SUBST(MEX)

dnl ----------------------------------------------------------------------------

AC_ARG_WITH([c-cpp],
	    [AS_HELP_STRING([--with-c-cpp=yes/no],
                    [Compile both C and C++ interfaces for SNOPT (default: no)]) ],
            [with_bothc=yes;], [with_bothc=no;] )

if test x"${with_bothc}" == xyes; then
  AC_MSG_CHECKING([for C interface])
  AC_MSG_RESULT([${with_bothc}])
  AC_MSG_CHECKING([for C++ interface])
  AC_MSG_RESULT([${with_bothc}])
else
   AC_ARG_WITH([c],[AS_HELP_STRING([--with-c=yes/no],
       [Compile C interfaces for SNOPT (default: no)]) ],
       [with_c=yes;], [with_c=no;] )
   AC_MSG_CHECKING([for C interface])
   AC_MSG_RESULT([${with_c}])

  AC_ARG_WITH([cpp],
  	    [AS_HELP_STRING([--with-cpp=yes/no],
                    [Compile C++ interfaces for SNOPT (default: no)]) ],
            [with_cpp=yes;], [with_cpp=no;] )
  AC_MSG_CHECKING([for C++ interface])
  AC_MSG_RESULT([${with_cpp}])
fi

if test x"${with_bothc}" != xno; then
  with_c=yes;
  with_cpp=yes;
fi

api_modules="";
if test x"${with_c}" != xno -o x"${with_cpp}" != xno; then
   case $FC in
     *gfortran* ) AC_MSG_WARN([C/C++ support requires gfortran 4.4 or higher]);
  esac
  api_modules="interfaces";
fi
AC_SUBST(with_c)
AC_SUBST(with_cpp)

dnl ----------------------------------------------------------------------------

AC_ARG_WITH([cmex],
            [AS_HELP_STRING([--with-cmex=yes/no],
                [Compile C Mex-executable; requires Matlab]) ],
            [with_cmex=${withval};], [with_cmex=no;] )
if test x"${with_cmex}" != xno; then
  if test x"${with_matlab}" == xno; then
    with_cmex=no;
    AC_MSG_WARN([CMEX requires Matlab, reconfigure --with-matlab])
  fi
fi


AC_ARG_WITH([csrc],
	    [AS_HELP_STRING([--with-csrc=yes/no],
            [Compile f2c'd C libraries for SNOPT (default: no)]) ],
            [with_csrc=${withval};], [with_csrc=no;] )
AC_MSG_CHECKING([for f2c/C support])
AC_MSG_RESULT([${with_csrc}])
if test "x${with_csrc}" == xyes; then
  AC_MSG_WARN([C support needs an F2C library])
fi

c_modules="";
if test x"${with_csrc}" != xno; then
  c_modules="csrc cexamples";
fi

AC_ARG_ENABLE([f2c],
              [AS_HELP_STRING([--enable-f2c],
              [Force us to build our own f2c library]) ],
	      [build_f2c=${enableval};], [build_f2c=no;] )
if test x"${build_f2c}" == xyes; then
  AC_MSG_NOTICE([we are force-building f2c])

  f2cBLD="yes";
  F2C="\${top_builddir}/bin/f2c";
  f2cLIB="\${top_builddir}/lib/libf2c.la";
  f2cINC="-I\${top_builddir}/include";
  f2c_enabled=yes;
else
  f2cBLD="no";
  f2c_enabled=no;
fi

AC_ARG_WITH([f2c],
	    [AS_HELP_STRING([--with-f2c=DIR],
            [Provide location of f2c; --without-f2c or --with-f2c=no
             will disable f2c support, and C and C++ modules will not be compiled]) ],
            [ with_f2c=${withval} ], [] )

dnl If we're not forcing the f2c build, check if user wants f2c
if test x"${build_f2c}" != xyes; then

  if test x"${with_csrc}" == xno; then
      dnl  We don't need f2c so why bother looking for it or building it.
      f2CBLD="no"; f2cINC=""; f2cLIB="";
      f2c_enabled=no;
      AC_MSG_CHECKING([whether to build f2c])
      AC_MSG_RESULT([no])

  else
    dnl We need f2c so try to find it
    if test x"${with_f2c}" != xno; then
      dnl  Look for f2c and build if we don't find it.
      if test x"${with_f2c}" != xyes; then
        CHECK_F2C([${with_f2c}])
      else
        CHECK_F2C([])
      fi
      f2c_enabled=yes;

    else
      dnl User has disabled f2c
      f2CBLD="no"; f2cINC=""; f2cLIB="";
      AC_MSG_CHECKING([whether to build f2c])
      AC_MSG_RESULT([no])
      AC_MSG_WARN([f2c disabled])
      with_csrc=no;
      f2c_enabled=no;
    fi
  fi
fi

AC_SUBST([f2cBLD])
AC_SUBST([F2C])
AC_SUBST([f2c_enabled])
AC_SUBST([f2cINC])
AC_SUBST([f2cLIB])

dnl ----------------------------------------------------------------------------

# Checking default integer size
prec=32
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for default size of Fortran integer], [my_cv_intsize],[
AC_RUN_IFELSE([
  AC_LANG_PROGRAM([],
     [
        use iso_fortran_env, only : file_storage_size
        integer :: nbytes, nunits
        inquire(iolength=nunits) 10
        nbytes = nunits*file_storage_size/8
        open(unit=42,file='conftestval',status='new')
        write(42,'(i1)') nbytes
        close(42)
     ])
  ],
  [ my_cv_intsize=`cat conftestval`
    if test x"${my_cv_intsize}" == x8; then
      prec=64
    fi
  ],
  [ AC_MSG_WARN([could not determine integer size -- continuing with 32-bit integers])
    prec=32 ])
])
AC_LANG_POP([Fortran])


dnl ----------------------------------------------------------------------------

# User-defined option will override above test
AC_ARG_WITH([32],
            [AS_HELP_STRING([--with-32],
                            [32-bit integer support (default is 32-bit)]) ],
 	    [prec=32;
             AC_MSG_NOTICE([32-bit integers enabled]) ], [])

AC_ARG_WITH([64],
            [AS_HELP_STRING([--with-64],
                            [64-bit integer support (default is 32-bit)]) ],
 	    [prec=64;
             AC_MSG_NOTICE([64-bit integers enabled]) ], [])

#AC_ARG_WITH([128],
#            [AS_HELP_STRING([--with-128],
#                            [quad-precision support (64-bit ints, 128-bit reals)]) ],
# 	    [prec=128;
#             AC_MSG_NOTICE([quad-precision enabled (64-bit ints, 128-bit reals)]) ], [])

AC_ARG_WITH([debug],
            [AS_HELP_STRING([--with-debug],
                            [Compile with debug flags])],
 	    [DBGFLAG="-g"; OPTFLAG="";],
	    [OPTFLAG="-O"; DBGFLAG="";])

dnl ----------------------------------------------------------------------------

FCFLAGS="${user_FCFLAGS}"
CFLAGS="${user_CFLAGS}"
CXXFLAGS="${user_CXXFLAGS}"
LDFLAGS="${user_LDFLAGS}"

if test x"${prec}" = x64; then
   PC_FLAGS="-DD64"

   case "${FC}" in
   *gfortran*)
     PF_FLAGS="-fdefault-integer-8"
   ;;
   *nagfor*)
     PF_FLAGS="-i8"
   ;;
   *ifort*)
     PF_FLAGS="-i8"
   ;;
   esac
elif test x"${prec}" = x128; then
   PC_FLAGS=" -DD128"

   case "${FC}" in
   *gfortran*)
     PF_FLAGS=" -fdefault-integer-8"
   ;;
   *nagfor* )
     PF_FLAGS=" -i8"
   ;;
   *ifort* )
     PF_FLAGS=" -i8"
   ;;
   esac
fi

case "${FC}" in
  *nagfor* )
    PF_FLAGS+=" -w"
  ;;
esac

if test -z "${FCFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    FCFLAGS="${OPTFLAG} ${PF_FLAGS}"
  else
    FCFLAGS="${DBGFLAG} ${PF_FLAGS}"
  fi
else
  FCFLAGS="${FCFLAGS} ${PF_FLAGS}"
fi

if test -z "${CFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    CFLAGS="${OPTFLAG} ${PC_FLAGS}"
  else
    CFLAGS="${DBGFLAG} ${PC_FLAGS}"
  fi
else
  CFLAGS="${CFLAGS} ${PC_FLAGS}"
fi

if test -z "${CXXFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    CXXFLAGS="${OPTFLAG}"
  else
    CXXFLAGS="${DBGFLAG}"
  fi
fi

if test -z "${MEXFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    MEXFLAGS="${OPTFLAG}"
  else
    MEXFLAGS="${DBGFLAG}"
  fi
fi

if test x"${with_matlab}" == xyes -a x"${blas}" == xyes; then
  MEXFLAGS="${MEXFLAGS} -largeArrayDims"
fi

AC_SUBST([prec])
AC_SUBST([FCFLAGS])
AC_SUBST([CFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([MEXFLAGS])

dnl ----------------------------------------------------------------------------

# Set modules.
if test x"${f2c_enabled}" == xno; then
  possible_modules="src examples ${api_modules}"
else
  possible_modules="src examples ${c_modules} ${api_modules}"
fi
modules=""
for module in ${possible_modules}; do
	if test -d ${srcdir}/${module}; then
		modules="${modules} ${module}"
	fi
done

if test -d ${srcdir}/mex; then
  if test x"${with_matlab}" == xyes; then
    modules="${modules} mex"
  fi
fi

if test -d ${srcdir}/cmex; then
  if test x"${with_cmex}" == xyes; then
    modules="${modules} cmex"
  fi
fi

AC_SUBST([modules])

dnl ----------------------------------------------------------------------------
case $FC in
  *nagfor* ) enable_shared=no;
esac

AC_SUBST([enable_shared])
AC_SUBST([enable_static])

AC_OUTPUT([GNUmakefile])

echo ""
echo "SNOPT will compile the following modules:"
echo "    ${modules}"
echo ""
echo "  Using external BLAS library?      ${blas}    ${BLAS_LIBS} "
echo ""
echo "  C   support?                      ${with_c}               "
echo "  C++ support?                      ${with_cpp}             "
echo "    (Note C/C++ support requires gfortran 4.4 or higher)    "
echo ""
if test x"${with_csrc}" == xyes; then
  echo "  Creating a C SNOPT library?       ${with_csrc}            "
  echo "  Building F2C library?             ${f2cBLD}               "
  echo "  Using F2C library?                ${f2c_enabled} ${f2cLIB}"
  echo "  (F2C version is no longer supported --"
  echo "      if possible, please use the new C/C++ interface)      "
  echo ""
fi
echo "  Matlab mex-executable available?  ${matlab}  ${MATLABDIR} "
if test x"${matlab}" == xyes; then
  echo "    (Make sure to use a compiler compatible with Matlab to create mex-files)"
  echo ""
fi
echo ""
echo "Libraries will be installed into ${prefix}/lib"
echo "  Shared libraries will be created: ${enable_shared}"
echo "  Static libraries will be created: ${enable_static}"
echo ""
echo ""
echo "To build all libraries              'make'          "
echo "To build the SNOPT Fortran library  'make snopt'    "
echo "To install the libraries            'make install'  "
echo "To build the examples               'make examples' "
echo ""
echo "To build the Matlab mex-file        'make matlab'   "
echo ""
if test x"${with_c}" == xyes -o x"${with_cpp}" == xyes; then
   echo "To build the C and C++ interfaces   'make interface'"
   echo "To build the C interface            'make c'"
   echo "To build the C examples             'make cexamples' "
   echo "To build the C++ interface          'make cpp'"
   echo "To build the C++ examples           'make cppexamples' "
   echo ""
fi
echo "To check the build                  'make check'    "
echo ""
