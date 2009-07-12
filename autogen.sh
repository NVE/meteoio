#!/bin/sh

clean ()
{
  # remove autotools stuff
  rm -f aclocal.m4 configure config.log config.status
  rm -rf autom4te.cache
}


#
# option checking
#
if test "x$1" = "xclean"; then
  set -x
  clean
  set +x
  exit 0
fi


dir=`echo "$0" | sed 's,[^/]*$,,'`
test "x${dir}" = "x" && dir='.'

if test "x`cd "${dir}" 2>/dev/null && pwd`" != "x`pwd`"
then
    echo "This script must be executed directly from the source directory."
    exit 1
fi

rm -f config.cache acconfig.h 

echo "Running aclocal..."       && \
aclocal

if test "$?" != "0" 
then
    echo "Can't find program 'aclocal' - is the autoconf package installed?"
    exit 1
fi

echo "Running autoconf..."
autoconf

if test "$?" != "0" 
then
    echo "Can't find program 'autoconf' - is the autoconf package installed?"
    exit 1
fi

echo
echo "Please run './configure' to create the Makefile"
echo

exit 0