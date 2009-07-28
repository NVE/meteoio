#include "BoschungIO.h"

using namespace std;

BoschungIO::BoschungIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){}

BoschungIO::BoschungIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	//Nothing else so far
}

//Copy constructor
//BoschungIO::BoschungIO(const BoschungIO& bio) : cfg(bio.cfg){
//  createBuffer();
//}

BoschungIO::BoschungIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	//Nothing else so far
}

BoschungIO::~BoschungIO() throw()
{
	cleanup();
}

void BoschungIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
}

//Clone function
//BoschungIO* BoschungIO::clone() const { return new BoschungIO(*this); }

void BoschungIO::read2DGrid(Grid2DObject&, const string& filename)
{
	//Nothing so far
	(void)filename;
	throw IOException("Nothing implemented here", AT);
}

void BoschungIO::readDEM(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BoschungIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BoschungIO::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  unsigned int stationindex)
{
	if (vecStationName.size() == 0)
		readStationNames(); //reads station names into vector<string> vecStationName

	unsigned int indexStart=0, indexEnd=vecStationName.size();

	//The following part decides whether all the stations are rebuffered or just one station
	if (stationindex == IOUtils::npos){
		vecMeteo.clear();
		vecStation.clear();
		
		vecMeteo.insert(vecMeteo.begin(), vecStationName.size(), vector<MeteoData>());
		vecStation.insert(vecStation.begin(), vecStationName.size(), vector<StationData>());
	} else {
		if ((stationindex < vecMeteo.size()) && (stationindex < vecStation.size())){
			indexStart = stationindex;
			indexEnd   = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("", AT);
		}
	}

	for (unsigned int ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		//cout << vecStationName[ii] << endl;
		bufferData(dateStart, dateEnd, vecMeteo, vecStation, ii);
	}
}

void BoschungIO::readStationNames()
{
	vecStationName.clear();

	//Read in the StationNames
	string xmlpath="", str_stations="";
	int stations=0;

	cfg.getValue("NROFSTATIONS", str_stations);

	if (!IOUtils::convertString(stations, str_stations, std::dec)) {
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);
	}
  
	for (int ii=0; ii<stations; ii++) {
		stringstream tmp_stream;
		string stationname="", tmp_file="";
		Date_IO tmp_date(0);
    
		tmp_stream << (ii+1); //needed to construct key name
		cfg.getValue(string("STATION"+tmp_stream.str()), stationname);

		vecStationName.push_back(stationname);
	}    
}

void BoschungIO::getFiles(const string& stationname, const Date_IO& start_date, const Date_IO& end_date, 
					 vector<string>& vecFiles, vector<Date_IO>& vecDate_IO)
{
	list<string> dirlist = list<string>();
	Date_IO tmp_date;
	string xmlpath="";

	cfg.getValue("XMLPATH", xmlpath);
	vecFiles.clear();
	IOUtils::readDirectory(xmlpath, dirlist, "_" + stationname + ".xml");

	//Sort dirlist in ascending order
	dirlist.sort();

	//Check date in every filename
	list<string>::iterator it = dirlist.begin(); 

	if (start_date > end_date){ //Special case return first data >= dateStart
		while ((it != dirlist.end())) {
			//check validity of filename
			if (validFilename(*it)) {
				string filename_out = *it;
				stringToDate_IO(filename_out, tmp_date);

				if (tmp_date > start_date) {
					vecFiles.push_back(xmlpath + "/" + filename_out);
					vecDate_IO.push_back(tmp_date);
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
			string filename_out = *it;
			stringToDate_IO(filename_out, tmp_date);

			if ((tmp_date >= start_date) && (tmp_date <= end_date)) {
				vecFiles.push_back(xmlpath + "/" + filename_out);
				vecDate_IO.push_back(tmp_date);
			}
		}

		it++;
	}
}

bool BoschungIO::bufferData(const Date_IO& dateStart, const Date_IO& dateEnd, 
					   std::vector< std::vector<MeteoData> >& vecMeteo, 
					   std::vector< std::vector<StationData> >& vecStation,					   
					   const unsigned int& stationnr)
{
	vector<string> vecFiles;
	vector<Date_IO> vecDate_IO;

	if (stationnr >= vecMeteo.size()) {
		throw IndexOutOfBoundsException("", AT);
	}

	vecMeteo[stationnr].clear();
	vecStation[stationnr].clear();

	getFiles(vecStationName[stationnr], dateStart, dateEnd, vecFiles, vecDate_IO);
	//cout << "[i] Buffering station number: " << vecStationName[stationnr] << "  " << vecFiles.size() << " files" << endl;

	if (vecFiles.size()==0) { //No files in range between dateStart and dateEnd
		return false;
	}

	for (unsigned int ii=0; ii<vecFiles.size(); ii++) {
		MeteoData meteoData;
		StationData stationData;
		xmlExtractData(vecFiles[ii], vecDate_IO[ii], meteoData, stationData);    
    
		vecMeteo[stationnr].push_back(meteoData);
		vecStation[stationnr].push_back(stationData);
	}

	return true;
}


void BoschungIO::checkForMeteoFiles(const string& xmlpath, const string& stationname, const Date_IO& date_in,
							 string& filename_out, Date_IO& date_out)
{
	list<string> dirlist = list<string>();
	IOUtils::readDirectory(xmlpath, dirlist, "_" + stationname + ".xml");

	//Sort dirlist in ascending order
	dirlist.sort();

	//Check date in every filename
	list<string>::iterator it = dirlist.begin(); 
  
	while ((it != dirlist.end()) && (date_out < date_in)) {
		//check validity of filename
		if (validFilename(*it)) {
			filename_out = *it;
			stringToDate_IO(filename_out, date_out);
			//cout << tmp_file << "  " << tmp_date.toString() << endl;
		}
    
		it++;
	}  
}

void BoschungIO::xmlExtractData(const string& filename, const Date_IO& date_in, MeteoData& md, StationData& sd)
{
	double ta=nodata, iswr=nodata, vw=nodata, dw=nodata, rh=nodata, lwr=nodata, nswc=nodata, tsg=nodata, tss=nodata, hs=nodata, rswr=nodata;
	double longitude=nodata, latitude=nodata, altitude=nodata;

	//Try to read xml file
	xmlpp::DomParser parser;
	//parser.set_validate(); provide DTD to check syntax
	parser.set_substitute_entities(); //We just want the text to be resolved/unescaped automatically.
	parser.parse_file(filename);
	if(parser) {
		//Walk the tree: ROOT NODE
		xmlpp::Node* pNode = parser.get_document()->get_root_node(); //deleted by DomParser.

		//Read in StationData
		sd.stationName = xmlGetNodeContent(pNode, "stationsId");
		string str_long = xmlGetNodeContent(pNode, "stationsLon");
		string str_lati = xmlGetNodeContent(pNode, "stationsLat");
		string str_alti = xmlGetNodeContent(pNode, "stationsAlt");

		xmlParseStringToDouble(str_long, longitude, "stationsLon");
		xmlParseStringToDouble(str_lati, latitude, "stationsLat");
		xmlParseStringToDouble(str_alti, altitude, "stationsAlt");
		sd.longitude = longitude;
		sd.latitude = latitude;
		sd.altitude = altitude;
		sd.eastCoordinate = nodata;
		sd.northCoordinate = nodata;

		//Air temperature
		string str_lt = xmlGetNodeContent(pNode, "lt");
		xmlParseStringToDouble(str_lt, ta, "lt");

		//gs = iswr
		string str_gs = xmlGetNodeContent(pNode, "gs");
		xmlParseStringToDouble(str_gs, iswr, "gs");

		//Wind velocity
		string str_vw = xmlGetNodeContent(pNode, "wgm");
		xmlParseStringToDouble(str_vw, vw, "wgm");

		//Relative humidity
		string str_rh = xmlGetNodeContent(pNode, "rlf");
		xmlParseStringToDouble(str_rh, rh, "rlf");

		//nswc
		string str_ns = xmlGetNodeContent(pNode, "ns");
		xmlParseStringToDouble(str_ns, nswc, "ns");

		//sb = lwr
		string str_sb = xmlGetNodeContent(pNode, "sb");
		xmlParseStringToDouble(str_sb, lwr, "sb");

		md.setMeteoData(date_in, ta, iswr, vw, dw, rh, lwr, nswc, tsg, tss, hs, rswr);
		convertUnits(md);
    
	} else {
		throw IOException("Error parsing XML", AT);
	}
}

void BoschungIO::xmlParseStringToDouble(const string& str_in, double& d_out, const string& parname)
{
	if (str_in!="") {//do nothing if empty content for a node was read
		if (!IOUtils::convertString(d_out, str_in, std::dec)) {//but if content of node is not empty, try conversion
			throw ConversionFailedException("Error while reading value for " + parname, AT);
		}
	}
}

std::string BoschungIO::xmlGetNodeContent(xmlpp::Node* pNode, const string& nodename)
{
	xmlpp::Node* tmpNode= xmlGetNode(pNode, nodename);
  
	if (tmpNode!=NULL) {
		xmlpp::Node* tmpNode2= xmlGetNode(tmpNode, "text"); //Try to retrieve text content
		if (tmpNode2!=NULL) {
			const xmlpp::TextNode* nodeText = dynamic_cast<const xmlpp::TextNode*>(tmpNode2);
			return string(nodeText->get_content());
		}
	}

	return string("");
}


xmlpp::Node* BoschungIO::xmlGetNode(xmlpp::Node* parentNode, const string& nodename)
{
	if (xmlGetNodeName(parentNode)==nodename) {
		return parentNode;
	}

	xmlpp::Node::NodeList list = parentNode->get_children();

	for(xmlpp::Node::NodeList::iterator iter = list.begin(); iter != list.end(); ++iter) {
		//xmlpp::Node* tmpNode= *iter;
		//cout << tmpNode->get_name() << endl;

		xmlpp::Node* tmpNode2= (xmlGetNode(*iter, nodename));
		if (tmpNode2!=NULL) {
			return tmpNode2;
		}
	}

	return NULL;
}

std::string BoschungIO::xmlGetNodeName(xmlpp::Node* pNode)
{
	string nodename = pNode->get_name();
	return nodename;
}

void BoschungIO::stringToDate_IO(const string& instr, Date_IO& date_out) const
{
	int tmp[5];

	string year = "20" + instr.substr(0,2);
	string month = instr.substr(2,2);
	string day = instr.substr(4,2);
	string hour = instr.substr(6,2);
	string minute = instr.substr(8,2);

	IOUtils::convertString(tmp[0], year, std::dec);
	IOUtils::convertString(tmp[1], month, std::dec);
	IOUtils::convertString(tmp[2], day, std::dec); 
	IOUtils::convertString(tmp[3], hour, std::dec);
	IOUtils::convertString(tmp[4], minute, std::dec);
  
	date_out.setDate(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4]);
}

bool BoschungIO::validFilename(const string& tmp) const
{
	size_t pos = tmp.find_first_not_of("0123456789");//Filename must start with 10 numbers
	if (pos!=10) {
		return false;
	}
  
	return true;
}

void BoschungIO::read2DMeteo(const Date_IO& date_in, vector<MeteoData>& meteo_out)
{
	vector<StationData> vecStation;
	read2DMeteo(date_in, meteo_out, vecStation);
}

void BoschungIO::read2DMeteo(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
{
	vecMeteo.clear();
	vecStation.clear();

	//Read in the StationNames
	string xmlpath="", str_stations="";
	int stations=0;

	cfg.getValue("NROFSTATIONS", str_stations);
	cfg.getValue("XMLPATH", xmlpath);

	if (!IOUtils::convertString(stations, str_stations, std::dec)) {
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);
	}
  
	for (int ii=0; ii<stations; ii++) {
		stringstream tmp_stream;
		string stationname="", tmp_file="";
		Date_IO tmp_date(0);
		MeteoData md;
		StationData sd;
    
		tmp_stream << (ii+1); //needed to construct key name
		cfg.getValue(string("STATION"+tmp_stream.str()), stationname);

		checkForMeteoFiles(xmlpath, stationname, date_in, tmp_file, tmp_date);
		//Check whether file was found
		if (tmp_date<date_in) {
			throw FileNotFoundException("No XML file in path '" + xmlpath 
								   + "' found for date " + date_in.toString() + " for station " + stationname, AT);
		}

		//Read in data from XML File
		xmlExtractData(xmlpath+"/"+tmp_file, tmp_date, md, sd);
    
		vecMeteo.push_back(md);
		vecStation.push_back(sd);
	}  
}

void BoschungIO::readAssimilationData(const Date_IO&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BoschungIO::readSpecialPoints(CSpecialPTSArray&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void BoschungIO::write2DGrid(const Grid2DObject&, const string& name)
{
	//Nothing so far
	(void)name;
	throw IOException("Nothing implemented here", AT);
}

void BoschungIO::convertUnits(MeteoData& meteo)
{
	//converts C to Kelvin, converts lwr to ea, converts RH to [0,1]
	if(meteo.ta==nodata) {
		meteo.ta=nodata;
	} else {
		meteo.ta=C_TO_K(meteo.ta);
	}
	
	if(meteo.tsg==nodata) {
		meteo.tsg=nodata;
	} else {
		meteo.tsg=C_TO_K(meteo.tss);
	}
	
	if(meteo.tss==nodata) {
		meteo.tss=nodata;
	} else {
		meteo.tss=C_TO_K(meteo.tss);
	}

	if(meteo.rh==nodata) {
		meteo.rh=nodata;
	} else {
		meteo.rh /= 100.;
	}
}

extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}
  
	void* loadObject(const string& classname, const string& filename) {
		if(classname == "BoschungIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new BoschungIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
