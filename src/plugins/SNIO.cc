/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include "SNIO.h"

/**
 * @page snowpack SNIO
 * @section snowpack_format Format
 * This is for reading grid data in the SNOWPACK meteo format (HACK: read implementation still missing)
 *
 * @section snowpack_units Units
 * The distances are assumed to be in meters.
 *
 * @section snowpack_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDIN: input coordinate system (see Coordinate)
 * - COORDPARAM: extra input coordinates parameters (see Coordinate)
 * - METEOPATH: path to the output directory
 *
 */
//HACK: to do!!
using namespace std;
using namespace mio;

const int SNIO::sn_julian_offset = 2415021;
const double SNIO::plugin_nodata = 0.0; //plugin specific nodata value

SNIO::SNIO(void (*delObj)(void*), const std::string& filename) : IOInterface(delObj), cfg(filename)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

SNIO::SNIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

SNIO::SNIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

SNIO::~SNIO() throw()
{
	cleanup();
}

void SNIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void SNIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*filename*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readStationData(const Date& /*date*/, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
					 std::vector< std::vector<MeteoData> >& /*vecMeteo*/, 
					 std::vector< std::vector<StationData> >& /*vecStation*/,
					 const unsigned int& /*stationindex*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
                          const std::vector< std::vector<StationData> >& vecStation,
                          const std::string&)
{
	string path="";
	cfg.getValue("METEOPATH", "Output", path);

	for(unsigned int ii=0; ii<vecMeteo.size(); ii++) {
		if(vecStation[ii].size()>0) {
			const std::string station_name = vecStation[ii][0].stationName;
			const std::string output_name = path + "/" + station_name + ".inp";
			if( !IOUtils::fileExists(output_name) ) {
				fout.open(output_name.c_str());
				writeStationHeader(vecMeteo[ii], station_name);
			} else {
				fout.open(output_name.c_str());
			}
			writeStationMeteo(vecMeteo[ii], output_name);
			fout.close();
		}
	}
}

void SNIO::writeStationHeader(const std::vector<MeteoData>& Meteo, const std::string station_name)
{
	//writing the (very basic) metadata
	fout << "MTO <" << station_name << "> " << Meteo.size() << "\n";
}

void SNIO::writeStationMeteo(const std::vector<MeteoData>& Meteo, const std::string& file_name)
{ //write out the data for 1 station
	unsigned int failure_count = 0;
	unsigned int Dirichlet_failure_count = 0;

	for(unsigned int ii=0; ii<Meteo.size(); ii++) {
		int YYYY, MM, DD, HH, MI;
		Meteo[ii].date.getDate(YYYY, MM, DD, HH, MI);
		const double sn_julian = Meteo[ii].date.getJulianDate() - sn_julian_offset + 0.5;
		const double ta = Meteo[ii].ta;
		const double rh = Meteo[ii].rh;
		const double hnw = Meteo[ii].hnw;
		const double vw = Meteo[ii].vw;
		const double dw = Meteo[ii].dw;
		const double iswr = Meteo[ii].iswr;
		const double rswr = Meteo[ii].rswr;
		const double ilwr = Meteo[ii].ilwr;
		const double tss = Meteo[ii].tss;
		const double tsg = Meteo[ii].tsg;
		const double hs = Meteo[ii].hs;

		fout.fill('0');
		fout << "M " << setw(2) << DD << "." << setw(2) << MM << "." << setw(4) << YYYY << " " << setw(2) << HH << ":" << setw(2) << MI << " ";
		fout.flags ( ios::fixed );
		fout << setprecision(6) << setw(12) << sn_julian << " ";

		//default formatting parameters for the measurements
		fout.flags ( ios::fixed );
		fout.fill(' ');
		fout.width(6);

		//TA, RH, VW, DW
		if(ta==IOUtils::nodata) {
			failure_count++;
			fout << setw(6) << setprecision(2) << ta << " ";
		} else
			fout << setw(6) << setprecision(2) << K_TO_C(Meteo[ii].ta) << " ";
		if(rh==IOUtils::nodata) {
			failure_count++;
			fout << setw(5) << setprecision(1) << rh << " ";
		} else
			fout << setw(5) << setprecision(1) << rh * 100. << " ";
		if(vw==IOUtils::nodata)
			failure_count++;
		fout << setw(4) << setprecision(1) << vw << " ";
		if(dw==IOUtils::nodata)
			failure_count++;
		fout << setw(3) << setprecision(0) << dw << " ";
		
		//ISWR, RSWR
		if(iswr==IOUtils::nodata && rswr==IOUtils::nodata) {
			failure_count++;
			fout << setw(3) << setprecision(0) << iswr << " " << setprecision(0) << rswr << " ";
		} else {
			if(iswr==IOUtils::nodata)
				fout << setw(3) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(3) << setprecision(0) << iswr << " ";
			if(rswr==IOUtils::nodata)
				fout << setw(3) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(3) << setprecision(0) << rswr << " ";
		}

		//LWR
		if(ilwr==IOUtils::nodata)
			failure_count++;
		fout << setw(3) << setprecision(0) << Meteo[ii].ilwr << " ";

		//TSS, TSG (only required for Dirichlet)
		if(tsg==IOUtils::nodata) {
			Dirichlet_failure_count++;
			fout << setw(6) << setprecision(1) << "0.0" << " ";
		} else {
			fout << setw(6) << setprecision(2) << K_TO_C(tsg) << " ";
		}
		if(tss==IOUtils::nodata) {
			Dirichlet_failure_count++;
			fout << setw(6) << setprecision(1) << "0.0" << " ";
		} else {
			fout << setw(6) << setprecision(2) << K_TO_C(tss) << " ";
		}

		//HNW, HS
		if(hnw==IOUtils::nodata && hs==IOUtils::nodata) {
			failure_count++;
			fout << setw(5) << setprecision(2) << hnw << " " << setprecision(3) << hs << " ";
		} else {
			if(hnw==IOUtils::nodata)
				fout << setw(5) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(5) << setprecision(2) << hnw << " ";
			if(hs==IOUtils::nodata)
				fout << setw(5) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(5) << setprecision(3) << hs / 100. << " ";
		}

		//we don't write any snow depth temperatures.
		//we can not write wind velocity at the wind station, but since it is optional...

		fout << endl;
	}

	fout << "END" << endl;

	if(failure_count>0 || Dirichlet_failure_count>0) {
		std::cout << "[W] " << failure_count << " (and potentially " << Dirichlet_failure_count <<
		" more) errors found when writing " << file_name << std::endl;
	}
}


void SNIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

#ifndef _METEOIO_JNI
extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}
  
	void* loadObject(const string& classname, const string& filename) {
		if(classname == "SNIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new SNIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif
