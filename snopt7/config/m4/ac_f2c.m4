# Find f2c library and headers.
AC_DEFUN([CHECK_F2C],
[
  sav_LDFLAGS=${LDFLAGS}
  sav_CXXFLAGS=${CXXFLAGS}

  myLIB=""
  myINC=""
  if test -n "$1"; then
    myLIB="-L$1/lib"
    myINC="-I$1/include"

    LDFLAGS="${myLIB} ${LDFLAGS}"
    CXXFLAGS="${myINC} ${CXXFLAGS}"
  fi

  AC_CHECK_PROG( [F2C], [f2c], [f2c], [${top_builddir}/bin/f2c])

  AC_CHECK_LIB( [f2c], [main],
                [ f2cBLD="no"; f2cLIB="${myLIB} -lf2c";
  		  AC_CHECK_HEADER( [f2c.h], [ f2cBLD="no"; f2cINC="${myINC}" ],
                                   [ f2cBLD="yes";
				     F2C="\${top_builddir}/bin/f2c";
                                     f2cLIB="\${top_builddir}/lib/libf2c.la";
                                     f2cINC="-I\${top_builddir}/include";
                  AC_MSG_WARN([f2c.h not found -- will build]) ]) ],
                [ f2cBLD="yes";
                  F2C="\${top_builddir}/bin/f2c";
                  f2cLIB="\${top_builddir}/lib/libf2c.la";
                  f2cINC="-I\${top_builddir}/include";
		  AC_MSG_WARN([f2c library not found -- will build]) ])

  AC_MSG_CHECKING([whether to build f2c])
  AC_MSG_RESULT([${f2cBLD}])
  LDFLAGS=${sav_LDFLAGS}
  CXXFLAGS=${sav_CXXFLAGS}
])
