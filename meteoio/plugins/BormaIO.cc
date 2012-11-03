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
#include "BormaIO.h"

using namespace std;

namespace mio {
/**
 * @page borma BORMA
 * @section borma_format Format
 * This plugin reads the XML files as generated by the Borma system. It requires libxml to compile and run. Since
 * years are expressed on two digits, all years greater than 80 are considered to belong to the XXth century while
 * years smaller belong to the XXIst century.
 *
 * @section borma_units Units
 * The units are assumed to be the following:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 *
 * @section borma_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - METEOPATH: string containing the path to the xml files
 * - STATION#: station id for the given number #
 */

const double BormaIO::plugin_nodata = -999.0; //plugin specific nodata value
const double BormaIO::default_tz = +1.; //default timezone
const double BormaIO::pivot_year = 80.; //pivot year for Y2K suppport
const std::string BormaIO::dflt_extension = ".xml";

BormaIO::BormaIO(const std::string& configfile) : cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = default_tz;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
}

BormaIO::BormaIO(const Config& cfgreader) : cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = default_tz;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
}

BormaIO::~BormaIO() throw()
{
	cleanup();
}

void BormaIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
}

void BormaIO::read2DGrid(Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::readDEM(DEMObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&,
                             const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::readStationData(const Date&, std::vector<StationData>&)
{
	//HACK: this method MUST be implemented for BufferedIOHandler to properly work
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                            std::vector< std::vector<MeteoData> >& vecMeteo,
                            const size_t& stationindex)
{
	if (vecStationName.size() == 0)
		readStationNames(); //reads station names into vector<string> vecStationName

	unsigned int indexStart=0, indexEnd=vecStationName.size();

	//The following part decides whether all the stations are rebuffered or just one station
	if (stationindex == IOUtils::npos){
		vecMeteo.clear();
		vecMeteo.insert(vecMeteo.begin(), vecStationName.size(), std::vector<MeteoData>());
	} else {
		if (stationindex < vecMeteo.size()){
			indexStart = stationindex;
			indexEnd   = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("", AT);
		}
	}

	for (unsigned int ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		bufferData(dateStart, dateEnd, vecMeteo, ii);
	}
}

void BormaIO::readStationNames()
{
	vecStationName.clear();

	size_t counter = 1;
	string stationname = "";

	do {
		stringstream ss;
		stationname = "";

		ss << "STATION" << counter;
		cfg.getValue(ss.str(), "Input", stationname, Config::nothrow);

		if (stationname != ""){
			vecStationName.push_back(stationname);
		}
		counter++;
	} while (stationname != "");

	nr_stations = counter - 1;
}

void BormaIO::getFiles(const std::string& stationname, const Date& start_date, const Date& end_date,
                       std::vector<std::string>& vecFiles, std::vector<Date>& vecDate)
{
	std::list<std::string> dirlist = std::list<std::string>();
	std::string xmlpath="";

	cfg.getValue("METEOPATH", "Input", xmlpath);
	vecFiles.clear();
	IOUtils::readDirectory(xmlpath, dirlist, "_" + stationname + dflt_extension);
	dirlist.sort();

	//Check date in every filename
	std::list<std::string>::iterator it = dirlist.begin();

	Date tmp_date(0., 0.);
	if (start_date > end_date){ //Special case return first data >= dateStart
		while ((it != dirlist.end())) {
			//check validity of filename
			if (validFilename(*it)) {
				const std::string filename_out = *it;
				tmp_date = stringToDate(filename_out);

				if (tmp_date > start_date) {
					vecFiles.push_back(xmlpath + "/" + filename_out);
					vecDate.push_back(tmp_date);
					return;
				}
			}

			it++;
		}
		return;
	}

	while ((it != dirlist.end()) && (tmp_date < end_date)) {
		//check validity of filename
		if (validFilename(*it)) {
			const std::string filename_out = *it;
			tmp_date = stringToDate(filename_out);

			if ((tmp_date >= start_date) && (tmp_date <= end_date)) {
				vecFiles.push_back(xmlpath + "/" + filename_out);
				vecDate.push_back(tmp_date);
			}
		}

		it++;
	}
}

bool BormaIO::bufferData(const Date& dateStart, const Date& dateEnd,
                         std::vector< std::vector<MeteoData> >& vecMeteo,
                         const unsigned int& stationnr)
{
	std::vector<std::string> vecFiles;
	std::vector<Date> vecDate;

	if (stationnr >= vecMeteo.size()) {
		throw IndexOutOfBoundsException("", AT);
	}

	vecMeteo[stationnr].clear();

	getFiles(vecStationName[stationnr], dateStart, dateEnd, vecFiles, vecDate);

	if (vecFiles.size()==0) { //No files in range between dateStart and dateEnd
		return false;
	}

	for (unsigned int ii=0; ii<vecFiles.size(); ii++) {
		MeteoData meteoData;
		StationData stationData;
		xmlExtractData(vecFiles[ii], vecDate[ii], meteoData, stationData);
		meteoData.meta = stationData;
		vecMeteo[stationnr].push_back(meteoData);
	}

	return true;
}


void BormaIO::checkForMeteoFiles(const std::string& xmlpath, const std::string& stationname, const Date& date_in,
                                 std::string& filename_out, Date& date_out)
{
	std::list<std::string> dirlist = std::list<std::string>();
	IOUtils::readDirectory(xmlpath, dirlist, "_" + stationname + ".xml");
	dirlist.sort();

	//Check date in every filename
	std::list<std::string>::iterator it = dirlist.begin();

	while ((it != dirlist.end()) && (date_out < date_in)) {
		//check validity of filename
		if (validFilename(*it)) {
			filename_out = *it;
			date_out = stringToDate(filename_out);
		}

		it++;
	}
}

void BormaIO::xmlExtractData(const std::string& filename, const Date& date_in, MeteoData& md, StationData& sd)
{
	double ta=IOUtils::nodata, iswr=IOUtils::nodata, vw=IOUtils::nodata, dw=IOUtils::nodata;
	double rh=IOUtils::nodata, ilwr=IOUtils::nodata, hnw=IOUtils::nodata, tsg=IOUtils::nodata;
	double tss=IOUtils::nodata, hs=IOUtils::nodata, rswr=IOUtils::nodata;
	double longitude=IOUtils::nodata, latitude=IOUtils::nodata, altitude=IOUtils::nodata;

	//Try to read xml file
	xmlpp::DomParser parser;
	//parser.set_validate(); //provide DTD to check syntax
	parser.set_substitute_entities(); //We just want the text to be resolved/unescaped automatically.
	parser.parse_file(filename);
	if(parser) {
		//Walk the tree: ROOT NODE
		xmlpp::Node* pNode = parser.get_document()->get_root_node(); //deleted by DomParser.

		//Read in StationData
		const std::string stationName = xmlGetNodeContent(pNode, "stationsId");
		const std::string str_long = xmlGetNodeContent(pNode, "stationsLon");
		const std::string str_lati = xmlGetNodeContent(pNode, "stationsLat");
		const std::string str_alti = xmlGetNodeContent(pNode, "stationsAlt");

		xmlParseStringToDouble(str_long, longitude, "stationsLon");
		xmlParseStringToDouble(str_lati, latitude, "stationsLat");
		xmlParseStringToDouble(str_alti, altitude, "stationsAlt");

		//HACK!! would it be possible for getValueForKey() to do this transparently? (with a user flag)
		latitude = IOUtils::standardizeNodata(latitude, plugin_nodata);
		longitude = IOUtils::standardizeNodata(longitude, plugin_nodata);
		altitude = IOUtils::standardizeNodata(altitude, plugin_nodata);

		//compute/check WGS coordinates (considered as the true reference) according to the projection as defined in cfg
		Coords location(coordin, coordinparam);
		location.setLatLon(latitude, longitude, altitude);
		sd.setStationData(location, stationName, stationName);

		//lt = ta
		const std::string str_lt = xmlGetNodeContent(pNode, "lt");
		xmlParseStringToDouble(str_lt, ta, "lt");

		//gs = iswr
		const std::string str_gs = xmlGetNodeContent(pNode, "gs");
		xmlParseStringToDouble(str_gs, iswr, "gs");

		//Wind velocity
		const std::string str_vw = xmlGetNodeContent(pNode, "wgm");
		xmlParseStringToDouble(str_vw, vw, "wgm");

		//rlf = rh
		const std::string str_rh = xmlGetNodeContent(pNode, "rlf");
		xmlParseStringToDouble(str_rh, rh, "rlf");

		//ni = hnw //TODO: not sure that this field really contains what we want...
		const std::string str_ns = xmlGetNodeContent(pNode, "ni");
		xmlParseStringToDouble(str_ns, hnw, "ni");

		//sb = ilwr
		const std::string str_sb = xmlGetNodeContent(pNode, "sb");
		xmlParseStringToDouble(str_sb, ilwr, "sb");

		//fbt = tss
		const std::string str_fbt = xmlGetNodeContent(pNode, "fbt");
		xmlParseStringToDouble(str_fbt, tss, "fbt");

		//sh = hs
		const std::string str_sh = xmlGetNodeContent(pNode, "sh");
		xmlParseStringToDouble(str_sh, hs, "sh");

		md.setDate(date_in);
		md(MeteoData::TA)   = ta;
		md(MeteoData::ISWR) = iswr;
		md(MeteoData::VW)   = vw;
		md(MeteoData::DW)   = dw;
		md(MeteoData::RH)   = rh;
		md(MeteoData::ILWR) = ilwr;
		md(MeteoData::HNW)  = hnw;
		md(MeteoData::TSG)  = tsg;
		md(MeteoData::TSS)  = tss;
		md(MeteoData::HS)   = hs;
		md(MeteoData::RSWR) = rswr;

		convertUnits(md);

	} else {
		throw IOException("Error parsing XML", AT);
	}
}

void BormaIO::xmlParseStringToDouble(const std::string& str_in, double& d_out, const std::string& parname)
{
	if (str_in!="") {//do nothing if empty content for a node was read
		if (!IOUtils::convertString(d_out, str_in, std::dec)) {//but if content of node is not empty, try conversion
			throw ConversionFailedException("Error while reading value for " + parname, AT);
		}
	}
}

std::string BormaIO::xmlGetNodeContent(xmlpp::Node* pNode, const std::string& nodename)
{
	xmlpp::Node* tmpNode= xmlGetNode(pNode, nodename);

	if (tmpNode!=NULL) {
		xmlpp::Node* tmpNode2= xmlGetNode(tmpNode, "text"); //Try to retrieve text content
		if (tmpNode2!=NULL) {
			const xmlpp::TextNode* nodeText = dynamic_cast<const xmlpp::TextNode*>(tmpNode2);
			return std::string(nodeText->get_content());
		}
	}

	return std::string("");
}


xmlpp::Node* BormaIO::xmlGetNode(xmlpp::Node* parentNode, const std::string& nodename)
{
	if (xmlGetNodeName(parentNode)==nodename) {
		return parentNode;
	}

	xmlpp::Node::NodeList list = parentNode->get_children();

	for(xmlpp::Node::NodeList::iterator iter = list.begin(); iter != list.end(); ++iter) {
		xmlpp::Node* tmpNode2= (xmlGetNode(*iter, nodename));
		if (tmpNode2!=NULL) {
			return tmpNode2;
		}
	}

	return NULL;
}

std::string BormaIO::xmlGetNodeName(xmlpp::Node* pNode)
{
	std::string nodename = pNode->get_name();
	return nodename;
}

Date BormaIO::stringToDate(const std::string& instr) const
{
	double year;
	IOUtils::convertString(year, instr.substr(0,2));

	Date date_out;
	//some Y2K support
	if(year<pivot_year)
		IOUtils::convertString(date_out, "20"+instr.substr(0,10) , in_tz);
	else
		IOUtils::convertString(date_out, "19"+instr.substr(0,10) , in_tz);

	return date_out;
}

bool BormaIO::validFilename(const std::string& tmp) const
{
	size_t pos = tmp.find_first_not_of("0123456789");//Filename must start with 10 numbers
	if (pos!=10) {
		return false;
	}

	return true;
}

void BormaIO::readAssimilationData(const Date&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::write2DGrid(const Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BormaIO::convertUnits(MeteoData& meteo)
{
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts RH to [0,1]
	double& ta = meteo(MeteoData::TA);
	if (ta != IOUtils::nodata)
		ta = C_TO_K(ta);

	double& tsg = meteo(MeteoData::TSG);
	if (tsg != IOUtils::nodata)
		tsg = C_TO_K(tsg);

	double& tss = meteo(MeteoData::TSS);
	if (tss != IOUtils::nodata)
		tss = C_TO_K(tss);

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata)
		rh /= 100.;
}

} //namespace
