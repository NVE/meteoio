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
aclocal                         && \
echo "Running autoconf..."      && \
autoconf

echo "Please run './configure' to create the Makefile"

exit 0