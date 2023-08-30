// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <sys/stat.h>
#include <sstream>
#include <regex>
#include <fstream>
#include <iostream>

#include <meteoio/FileUtils.h>
#include <meteoio/FStream.h>

namespace mio {

std::string cutPathToCWD(const std::string &path){
	std::string outpath = "";
	if (FileUtils::isAbsolutePath(path)){
		std::string cwd = FileUtils::getCWD();
		std::stringstream cwds(cwd);
		std::vector<std::string> cwd_parts;
		std::string item;
		while (std::getline(cwds, item, '/')){
			cwd_parts.push_back(item);
		}

		std::stringstream ps(path);
		while (std::getline(ps,item,'/')){
			if (std::find(cwd_parts.begin(), cwd_parts.end(), item)!=cwd_parts.end()){
				continue;
			}
			else {
				outpath+=item;
				outpath+="/";
			}
		}
		if (outpath.empty()) return "./";
		else if (outpath.back()=='/') outpath.pop_back();
	}
	else {
		std::regex e("\\.\\.\\/");
		outpath = std::regex_replace(path, e, "");
		std::regex e2("\\.\\/");
		outpath = std::regex_replace(outpath, e2,"");
	}
	outpath = "./"+ outpath;
	return outpath;
}

std::string limitAccess(const char* filename, bool write_directories){
    std::string path = FileUtils::getPath(filename);
	std::string file = FileUtils::getFilename(filename);
	if (write_directories) {
#ifdef LIMIT_WRITE_ACCESS
		if (FileUtils::directoryExists(path)	) {
			path = cutPathToCWD(FileUtils::cleanPath(path,true));
		} else {
			path = cutPathToCWD(path);
		}
#endif
		createTree((path+"/"+file).c_str());
	} else {
#ifdef LIMIT_WRITE_ACCESS
		throw IOException("Write access was limited, but directories are not supposed to be created, which makes it impossible to limit the output directories");
#endif
	}

	std::string FILE = path+"/"+file;
    return FILE;
}

void createTree(const char* filename, bool verbose){
    if (FileUtils::directoryExists(FileUtils::getPath(filename,false)))
    {
    }
    else
    {
        FileUtils::createDirectories(FileUtils::getPath(filename,false), verbose);
    }
}

std::string ofilestream::initialize(const char* filename) {
	return limitAccess(filename);
}

std::string ofilestream::initialize(const char* filename, const Config& cfgreader) {
	bool write_directories = cfgreader.get("WRITE_DIRECTORIES", "Output", true);
	return initialize(filename, write_directories);
}

std::string ofilestream::initialize(const char* filename, bool write_directories) {
	std::cerr << "WRITE_DIRECTORIES: " << write_directories << std::endl;
	return limitAccess(filename, write_directories);
}

void ofilestream::open(const char* filename, std::ios_base::openmode mode) {
	std::string FILE = limitAccess(filename);
    std::ofstream::open(FILE.c_str(), mode);
}

ofilestream::ofilestream(const char* filename, std::ios_base::openmode mode) : std::ofstream(initialize(filename).c_str(), mode) {
}

ofilestream::ofilestream(const std::string filename, std::ios_base::openmode mode) : std::ofstream(initialize(filename.c_str()).c_str(), mode) {
}

ofilestream::ofilestream(const char* filename, const Config& cfgreader, std::ios_base::openmode mode) : std::ofstream(initialize(filename, cfgreader).c_str(),mode) {
}

ofilestream::ofilestream(const std::string filename, const Config& cfgreader, std::ios_base::openmode mode) : std::ofstream(initialize(filename.c_str(), cfgreader).c_str(), mode) {
}
ofilestream::ofilestream(const char* filename, bool write_directories, std::ios_base::openmode mode) : std::ofstream(initialize(filename, write_directories).c_str(),mode) {
}


}
