#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.6])
AC_INIT(meteoio, 1.1.2, bavay@slf.ch)

m4_include(m4/ax_cflags_warn_all.m4)
m4_include(m4/m4_ax_cflags_gcc_option.m4)
m4_include(m4/m4_ax_compiler_vendor.m4)
m4_include(m4/pkg.m4)
#m4_include(m4/m4_ax_check_compiler_flags.m4)
#m4_include(m4/m4_ax_gcc_x86_cpuid.m4)
#m4_include(m4/m4_ax_gcc_archflag.m4)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_CHECK_PROG([COLORGCC], [colorgcc], [yes], [no])
if test "$COLORGCC" == "yes"
then
	AC_MSG_NOTICE([Compiling using colorgcc])
	CC=colorgcc
	#CXX=colorgcc
fi

#AC_PROG_LIBTOOL
AC_LANG(C++)

#testing GCC flags, assuming that c and c++ compilers are the same
AX_CFLAGS_WARN_ALL(WARNING)
AC_SUBST(WARNING)
AC_ARG_WITH(gcc-arch, [AC_HELP_STRING([--with-gcc-arch=<arch>], [use architecture <arch> for gcc -march/-mtune, instead of guessing])], GCC_ARCH=$withval, GCC_ARCH=yes)
if test "x$GCC_ARCH" != xyes
then
		ARCH="-march=$GCC_ARCH"
fi

AX_COMPILER_VENDOR
if test "$ax_cv_cxx_compiler_vendor" == "gnu"
then
	AC_MSG_NOTICE([Checking specific GNU flags])
	AX_CFLAGS_GCC_OPTION(-Wextra, EXTRA_WARN)
	AC_SUBST(EXTRA_WARN)

	if test "x$GCC_ARCH" = xyes
	then
		AX_CFLAGS_GCC_OPTION(-march=native, ARCH)
	fi

	AX_CFLAGS_GCC_OPTION(-g, OPTIM)
	AX_CFLAGS_GCC_OPTION(-O3, OPTIM)
	AX_CFLAGS_GCC_OPTION(-ftree-vectorize, OPTIM)
	AC_SUBST(OPTIM)

	AX_CXXFLAGS_GCC_OPTION(-combine, MAKE_OPTIM)
	AX_CXXFLAGS_GCC_OPTION(-pipe, MAKE_OPTIM)
	AC_SUBST(MAKE_OPTIM)
fi
AC_SUBST(ARCH)

# Checks for libraries.
AC_CHECK_LIB([dl], [main],,[AC_MSG_ERROR([Error, the dl library is missing!])])
AC_CHECK_LIB([c], [main],,[AC_MSG_ERROR([Error, the c library is missing!])])
AC_CHECK_LIB([m], [main],,[AC_MSG_ERROR([Error, the math library is missing!])])

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([dlfcn.h memory.h stdlib.h string.h sys/stat.h math.h],,[AC_MSG_ERROR([The header files must be present!])])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
AC_C_CONST
#AC_C_INLINE
#AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_STAT

AC_CHECK_PROG([POPCC], [popcc], [yes], [no])
if test "$POPCC" == "no"
then
	AC_MSG_NOTICE([Pop-C++ not found, you won't be able to use the parallel version])
fi

echo
echo "================================="
echo "=         USER OPTIONS          ="
echo "================================="
######################## DEBUG MODE ######################
AC_ARG_ENABLE(debug, [  --enable-debug	  Enable debug mode], [debug=yes])
if test "$debug" == "yes"
then
   AC_SUBST(MODE,debug)
else
   AC_SUBST(MODE,release)
fi

AC_MSG_NOTICE([Compiling for mode '$MODE'])

######################## PLUGINS ##########################
CPPFLAGSOLD="$CPPFLAGS"
LDFLAGSOLD="$LDFLAGS"

echo ""
echo "================================="
echo "=       ACTIVATED PLUGINS       ="
echo "================================="
AC_ARG_WITH(proj4,	[  --with-proj4	  enable proj4 geographic coordinates conversions (libproj required)])
AC_SUBST(PROJ4,no)
AC_ARG_WITH(bormaio, [  --with-bormaio	  enable BormaIO (libxml++-2.6 required)])
AC_SUBST(BORMAIO,no)
AC_ARG_WITH(gsnio,      [  --with-gsnio	          enable GSNIO (gsoap++ required)])
AC_SUBST(GSNIO,no)
AC_ARG_WITH(imisio,	[  --with-imisio           enable IMISIO ($ORACLE_HOME needs to be set)])
AC_SUBST(IMISIO,no)
AC_ARG_WITH(geotopio,	[  --with-geotopio         enable GeotopIO ], [use_geotopio=$withval], [use_geotopio=no])
AC_SUBST(GEOTOPIO,$use_geotopio)
AC_ARG_WITH(arcio,	[  --without-arcio         disable ARCIO plugin (default: enabled)], [use_arcio=$withval], [use_arcio=yes])
AC_SUBST(ARCIO,$use_arcio)
AC_ARG_WITH(snio,	[  --without-snio         disable SNIO plugin (default: enabled)], [use_snio=$withval], [use_snio=yes])
AC_SUBST(SNIO,$use_snio)
AC_ARG_WITH(smetio,	[  --without-smetio         disable SMETIO plugin (default: enabled)], [use_smetio=$withval], [use_smetio=yes])
AC_SUBST(SMETIO,$use_snio)
AC_ARG_WITH(grassio,	[  --without-grassio         disable GrassIO plugin (default: enabled)], [use_grassio=$withval], [use_grassio=yes])
AC_SUBST(GRASSIO,$use_grassio)
AC_ARG_VAR(ORACLE_HOME, [The home of your Oracle client installation])

if test "$with_geotopio" == "yes"; then AC_MSG_NOTICE([GeotopIO enabled. No prerequisites...]); fi
if test "x$use_grassio" == "xyes"; then AC_MSG_NOTICE([GrassIO enabled. No prerequisites...]); fi
if test "x$use_arcio" == "xyes"; then AC_MSG_NOTICE([ARCIO enabled. No prerequisites...]); fi
if test "x$use_snio" == "xyes"; then AC_MSG_NOTICE([SNIO enabled. No prerequisites...]); fi
if test "x$use_smetio" == "xyes"; then AC_MSG_NOTICE([SMETIO enabled. No prerequisites...]); fi

if test "$with_proj4" == "yes"
then
	AC_MSG_NOTICE([PROJ4 enabled. Checking for prerequisites...])
	AC_CHECK_HEADER([proj_api.h],,[AC_MSG_ERROR([Error, the proj_api.h header file is missing!])])
	AC_CHECK_LIB(proj,main,,[AC_MSG_ERROR([Error, the proj4 library is missing!])])

	AC_SUBST(PROJ4_CFLAGS)
	AC_SUBST(PROJ4_LIBS)
	AC_SUBST(PROJ4,yes)
fi

if test "$with_gsnio" == "yes"
then
	AC_MSG_NOTICE([GSNIO enabled. Checking for prerequisites...])
	AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
	if test "$PKG_CONFIG" == "no"; then
		echo "*** The pkg-config script could not be found. Make sure it is"
		echo "*** in your path, or set the PKG_CONFIG environment variable"
		echo "*** to the full path to pkg-config."
		echo "*** Or see http://www.freedesktop.org/software/pkgconfig to get pkg-config."
		rm -f configure config.status
		exit 1
	fi
	PKG_CHECK_MODULES(GSOAPPP, gsoap++ >= 2.0)
	CPPFLAGS="$CPPFLAGS $GSOAPPP_CFLAGS"

	AC_SUBST(GSOAPPP_CFLAGS)
	AC_SUBST(GSOAPPP_LIBS)

	AC_SUBST(GSNIO,yes)
fi

if test "$with_bormaio" == "yes"
then
	AC_MSG_NOTICE([BormaIO enabled. Checking for prerequisites...])
	AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
	if test "$PKG_CONFIG" == "no"; then
		echo "*** The pkg-config script could not be found. Make sure it is"
		echo "*** in your path, or set the PKG_CONFIG environment variable"
		echo "*** to the full path to pkg-config."
		echo "*** Or see http://www.freedesktop.org/software/pkgconfig to get pkg-config."
		rm -f configure config.status
		exit 1
	fi
	PKG_CHECK_MODULES(XMLPP, libxml++-2.6 >= 2.8)
	CPPFLAGS="$CPPFLAGS $XMLPP_CFLAGS"

	AC_CHECK_HEADER([libxml++/libxml++.h],,[AC_MSG_ERROR([Error, the libxml++.h header file is missing!])])
	#AC_CHECK_LIB(xml2, xmlParseFile,,[AC_MSG_ERROR([Error, the xml++ library is missing!])])

	AC_SUBST(XMLPP_CFLAGS)
 	AC_SUBST(XMLPP_LIBS)

	AC_SUBST(BORMAIO,yes)
fi

if test "$with_imisio" == "yes"
then
	#use AX_LIB_ORACLE_OCI?
	echo ""
	AC_MSG_NOTICE([IMISIO enabled. Checking for prerequisites...])
	CPPFLAGS="$CPPFLAGS -I$ORACLE_HOME/client_1/include"
	LDFLAGS="$LDFLAGS -L$ORACLE_HOME/client_1/lib  -lclntsh -lnnz11"
	AC_CHECK_HEADER([occi.h],,[AC_MSG_ERROR([Error, the occi.h header file is missing!])])
	AC_CHECK_LIB([occi], [main],,[AC_MSG_ERROR([Error, the occi library is missing!])])
	AC_SUBST(IMISIO,yes)
fi

AC_SUBST(CPPFLAGS, $CPPFLAGSOLD)
AC_SUBST(LDFLAGS,$LDFLAGSOLD)
AC_SUBST(prefix)
AC_SUBST(execprefix)

echo
AC_OUTPUT(Makefile)
echo "Run 'make' to compile the meteoio library"