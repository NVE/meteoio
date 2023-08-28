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

std::string limitAccess(const char* filename){
    std::string path = FileUtils::getPath(filename);
	std::string file = FileUtils::getFilename(filename);
#ifdef LIMIT_WRITE_ACCESS
	std::cout << "from limited access" << std::endl;
	path = cutPathToCWD(path);
	std::cout << path << std::endl;
#endif
	createTree((path+"/"+file).c_str());
	std::string FILE = path+"/"+file;
    return FILE;
}

void createTree(const char* filename, bool verbose){
	struct stat sb;
    if (stat(mio::FileUtils::getPath(filename,false).c_str(), &sb) == 0)
    {
    }
    else
    {
        FileUtils::createDirectories(FileUtils::getPath(filename,false), verbose);
    }
}

std::string ofstream::initialize(const char* filename){
	return limitAccess(filename);
}

void ofstream::open(const char* filename, std::ios_base::openmode mode){
	std::string FILE = limitAccess(filename);
    std::ofstream::open(FILE.c_str(), mode);
}

ofstream::ofstream(const char* filename, std::ios_base::openmode mode) : std::ofstream(initialize(filename).c_str(), mode)
{
}

}
