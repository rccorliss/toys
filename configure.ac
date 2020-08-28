AC_INIT(MYTREE,[2.0])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_F77(gfortran)
AC_PROG_CXX(CC g++)
AC_ENABLE_STATIC(no)

LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -std=c++11 -Wall -Werror"
fi

dnl test for root 6
if test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi
AM_CONDITIONAL([MAKEROOT6],[test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
