#!/bin/bash
#This script creates a new release of Alpine3D and packages it.
#todo: use cmake?!
#prerequisites: the path to Alpine3D documentation, the path to A3D_view

#parsing command line arguments
if [ $# == 0 ]; then
	printf "$0 {op|r} {src?}\n"
	exit 0
fi
MODE=$1
case "${MODE}" in
	"r")	TYPE=0 ;;
	*)
		printf "Unknown option $1\n"
		exit 0
	;;
esac
if [ $# == 1 ]; then
	WITH_SRC=0
else
	WITH_SRC=1
fi

#source files for SNOWPACK
src_files="src/*.cc src/*.h src/*.ph Makefile"
test_files="test/*.cc test/Makefile"
tools_files="tools/*.cc"
others=""

MeteoIO_root="."

case "${MODE}" in 
	"r")	file_list=`eval ls ${src_files} ${test_files} ${tools_files} ${others}` ;;
esac

#destination directories
doc="doc"
libs="lib"

function make_temp {
#temporary storage (will be cleaned up at the end)
	tmp="/tmp/MeteoIO_Release.$$"
	if [ -d ${tmp} ]
	then
		rm -R "${tmp}"
	fi
	mkdir "${tmp}"
}

function make_directories {
#creates the directory structure
	release=$(./src/version.sh)
	dest="${tmp}/MeteoIO_${release}"
	mkdir "${dest}"
	mkdir "${dest}/${libs}"
	mkdir "${dest}/${doc}"
}

function make_src_directories {
 	#mkdir "${dest}/${src}"
	subdir_list=$(printf "${file_list}" | xargs -i dirname {} | sort -u)
	for src_dir in ${subdir_list}
	do
		if [ ${src_dir} = "." ]; then true; else
			mkdir ${dest}/${src_dir}
		fi
	done
}

function fill_up_src {
#filling up the src directory
	for fichier in ${file_list}
	do
		subdir=`dirname "${fichier}"`
		if [ ${subdir} = "." ]
		then
			cp ${fichier} ${dest}/${fichier}
		else
			#cp ${fichier} ${dest}/${subdir}/${fichier} #crappy HACK
			cp ${fichier} ${dest}/${fichier}
		fi
	done
}

function create_archives {
#creating the archives
	datum=$(date "+%Y%m%d")
	case "${MODE}" in 
		"r")	release_name="MeteoIO" ;;
	esac
	release_name="${release_name}_${datum}"
	if [ ${WITH_SRC} -eq 0 ]; then
		release_name="${release_name}-bin"
	fi
	MeteoIO_root=`pwd`
	cd ${tmp}
	tar -c `basename ${dest}` | gzip -9 > ${MeteoIO_root}/${release_name}.tgz
	cd ${MeteoIO_root}
}

function cleanup_tmp {
#removing temp files
	rm -Rf ${tmp}
}

################################# Main ################################
make_temp
make_directories
if [ ${WITH_SRC} -eq 1 ]; then
	make_src_directories
	fill_up_src
fi

create_archives
cleanup_tmp

#done!