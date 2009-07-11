#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(meteoio, 0.1, bavay@slf.ch)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_LANG(C++)

# Checks for libraries.
AC_CHECK_LIB([dl], [main])

# Checks for header files.
AC_CHECK_HEADERS([occi.h dlfcn.h dirent.h sys/stat.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_REALLOC
#AC_CHECK_FUNCS([floor pow sqrt])


######################## DEBUG MODE ######################
AC_ARG_ENABLE(debug, [  --enable-debug	  Enable debug mode], [debug=yes])
if test "$debug" == "yes"
then
   AC_SUBST(MODE,debug)
else 
   AC_SUBST(MODE,release)
fi

######################## TARGET ##########################
AC_ARG_WITH(dest, [  --with-target=ARG	  where ARG is either 'safe' (default), 'grid' or 'zeus'])
AC_SUBST(DEST,safe)

if test "$with_target" == "safe" 
then
   AC_SUBST(DEST,safe)
elif test "$with_target" == "grid" 
then 
   AC_SUBST(DEST,grid)
elif test "$with_target" == "zeus" 
then 
   AC_SUBST(DEST,zeus)
elif test "$with_target" == "yes" 
then
   AC_MSG_ERROR([Please specify a valid argument for --with-target=ARG option])
fi


######################## PLUGINS ##########################
echo ""
echo "================================="
echo "=       ACTIVATED PLUGINS       ="
echo "================================="
AC_ARG_WITH(boschungio, [  --with-boschungio	  enable BoschungIO (libxml++-2.6 required)])
AC_SUBST(BOSCHUNGIO,no)
AC_ARG_WITH(imisio, [  --with-imisio	  enable IMISIO])
AC_SUBST(IMISIO,no)

if test "$with_boschungio" == "yes"
then
	AC_MSG_NOTICE([BoschungIO enabled. Checking for prerequisites...])
	PKG_CHECK_MODULES(XMLPP, libxml++-2.6 >= 2.8)
	CPPFLAGS="$CPPFLAGS $XMLPP_CFLAGS"
	AC_CHECK_HEADERS(libxml++/libxml++.h)
	AC_SUBST(XMLPP_CFLAGS)
 	AC_SUBST(XMLPP_LIBS)

	if test "$ac_cv_header_libxmlpp_libxmlpp_h" == "no"
	then
		AC_MSG_ERROR([Error, the libxml++.h header file is missing!])
	fi
	AC_SUBST(BOSCHUNGIO,yes)
fi

if test "$with_imisio" == "yes"
then
	AC_SUBST(IMISIO,yes)
fi


#AC_CONFIG_FILES([test/Makefile])
AC_OUTPUT(Makefile)
