#include "GSNIO.h"

using namespace std;

GSNIO::GSNIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){
	initGSNConnection();
}

GSNIO::GSNIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	initGSNConnection();
}

GSNIO::GSNIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	initGSNConnection();
}

GSNIO::~GSNIO() throw(){}

void GSNIO::initGSNConnection(){
	hostname = port = userid = passwd="";
	proxyport = -1;

	//soap_init(&gsn);
	//soap_init2(&gsn, SOAP_IO_KEEPALIVE, SOAP_IO_KEEPALIVE);

	/*
	 * Trying to read proxy settings: 
	 * - Firstly the hostname and port (both have to be provided).
	 * - If this succeeds then the username and password will be read
	 * - parameters not set will be set to ""
	 */
	try {
		cfg.getValue("PROXY", hostname); 
		cfg.getValue("PROXYPORT", port); 

		if (!IOUtils::convertString(proxyport, port, std::dec))
			throw ConversionFailedException("", AT);
		if (proxyport < 1) 
			throw IOException("",AT);

		gsn.proxy_host = hostname.c_str();
		gsn.proxy_port = proxyport;

		cfg.getValue("PROXYUSER", userid); 
		gsn.proxy_userid = userid.c_str();
		cfg.getValue("PROXYPASS", passwd); 
		gsn.proxy_passwd = passwd.c_str(); 
	} catch(std::exception& e){}
}

void GSNIO::read2DGrid(Grid2DObject&, const string& filename)
{
	//Nothing so far
	(void)filename;
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readDEM(DEMObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex)
{
	if (vecStationName.size() == 0)
		readStationNames(); //reads station names into vector<string> vecStationName

	if (vecStationName.size() == 0) //if there are no stations -> return
		return;

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
			throw IndexOutOfBoundsException("You tried to access a stationindex in readMeteoData that is out of bounds", AT);
		}
	}

	for (unsigned int ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		StationData sd;
		readStationMetaData(sd, ii);

		readData(dateStart, dateEnd, vecMeteo, vecStation, sd, ii);
	}
}

void GSNIO::readStationMetaData(StationData& sd, const unsigned int& stationindex)
{
	_ns1__getSensorLocation sensorloc_req;
	_ns1__getSensorLocationResponse sensorloc;

	sensorloc_req.sensor = &vecStationName[stationindex];
	
	if (gsn.getSensorLocation(&sensorloc_req, &sensorloc) == SOAP_OK){
		if (sensorloc.return_.size() != 4) throw IOException("Not enough SensorLocation data received ...", AT);

		string *s0 = &sensorloc.return_[0]; //easier to type
		string *s1 = &sensorloc.return_[1]; //easier to type
		string *s2 = &sensorloc.return_[2]; //easier to type
		string *s3 = &sensorloc.return_[3]; //easier to type
		size_t sep0 = s0->find_first_of("=");
		size_t sep1 = s1->find_first_of("=");
		size_t sep2 = s2->find_first_of("=");
		size_t sep3 = s3->find_first_of("=");
		double latitude=IOUtils::nodata, longitude=IOUtils::nodata, altitude=IOUtils::nodata;
		string name = s0->substr((sep0+1), (s0->size()-sep0-2));

		convertStringToDouble(latitude, s1->substr((sep1+1), (s1->size()-sep1-2)), "Latitude");
		convertStringToDouble(longitude, s2->substr((sep2+1), (s2->size()-sep2-2)), "Longitude");
		convertStringToDouble(altitude, s3->substr((sep3+1), (s3->size()-sep3-2)), "Altitude");

		sd.setStationData(IOUtils::nodata, IOUtils::nodata, altitude, name, latitude, longitude);
	} else {
		soap_print_fault(&gsn, stdout);
		throw IOException("Error in communication with GSN",AT);
	}
}

void GSNIO::parseString(const std::string& _string, std::vector<std::string>& vecString, MeteoData& md){
	vecString.clear();

	//cout << _string << endl;
	
	stringstream ss(_string);
	string tmpstring;

	while (std::getline(ss, tmpstring, ';')){
		stringstream data(tmpstring);
		while (std::getline(data, tmpstring, '=')){
			string key = tmpstring;
			if (!(std::getline(data, tmpstring, '=')))
				throw InvalidFormatException("",AT);

			if (key == "LIGHT") convertStringToDouble(md.iswr, tmpstring, "ISWR");
			else if (key == "TEMPERATURE") convertStringToDouble(md.ta, tmpstring, "Air Temperature");				
			else if (key == "TIMED") {
				tmpstring = tmpstring.substr(0,tmpstring.length()-3); //cut away the seconds
				time_t measurementTime;
				if (!IOUtils::convertString(measurementTime, tmpstring, std::dec))
					throw ConversionFailedException("Conversion failed for value TIMED", AT);
				md.date.setDate(measurementTime);
			}
		}
	}
}

void GSNIO::convertStringToDouble(double& d, const string& _string, const string& _parname){
	if (!IOUtils::convertString(d, _string, std::dec))
		throw ConversionFailedException("Conversion failed for value " + _parname, AT);
}

void GSNIO::readData(const Date_IO& dateStart, const Date_IO& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, 
				 std::vector< std::vector<StationData> >& vecStation, const StationData& sd, const unsigned int& stationindex)
{
	_ns1__getMeteoData meteodata_req;
	_ns1__getMeteoDataResponse meteodata;
	vector<string> vecString;

	meteodata_req.sensor = &vecStationName[stationindex];
	/*
	  Date_IO dateStart1(time(NULL));
	  dateStart1 -= Date_IO(0.013);
	  Date_IO dateEnd1(time(NULL));

	  LONG64 l1(dateStart1.getEpochTime());
	  LONG64 l2(dateEnd1.getEpochTime());
	*/

	LONG64 l1(dateStart.getEpochTime());
	LONG64 l2(dateEnd.getEpochTime());

	l1*=1000; //GSN is using ms, not seconds
	l2*=1000; //GSN is using ms, not seconds

	//cout << dateStart1 << "  " << dateEnd1 << endl;
	//cout << dateStart1.getEpochTime() << "==" << l1 << endl; cout << dateEnd1.getEpochTime() << "==" << l2 << endl;
	
	meteodata_req.from = l1;
	meteodata_req.to = l2;
	
	//cout << std::dec << meteodata_req.param1 << "\t" << meteodata_req.param2 << endl;

	if (gsn.getMeteoData(&meteodata_req, &meteodata) == SOAP_OK){
		cout << "\t[D] GSN: nr of datasets received for station "<< stationindex << ": " << meteodata.return_.size() << endl;
		for (unsigned int jj=0; jj<meteodata.return_.size(); jj++){
			MeteoData md;
			parseString(meteodata.return_[jj], vecString, md);
			convertUnits(md);

			//cout << md.toString() << endl;
			vecMeteo[stationindex].push_back(md);
			vecStation[stationindex].push_back(sd);
		}
	} else {
		soap_print_fault(&gsn, stderr);
		throw IOException("Error in communication with GSN",AT);
	}
}

void GSNIO::readStationNames()
{
	vecStationName.clear();
	_ns1__getSensorsResponse sensors;
	if (gsn.getSensors(&sensors) == SOAP_OK){
		for (unsigned int ii=0; ii<sensors.return_.size(); ii++){
			//cout << "[d] Sensor " << ii << " Name: " << sensors.return_[ii] << endl;
			vecStationName.push_back(sensors.return_[ii]);
		}
	} else {
		soap_print_fault(&gsn, stderr);
		throw IOException("Error in communication with GSN",AT);
	}
}

void GSNIO::readAssimilationData(const Date_IO&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readSpecialPoints(CSpecialPTSArray&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::write2DGrid(const Grid2DObject&, const string& name)
{
	//Nothing so far
	(void)name;
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::convertUnits(MeteoData& meteo)
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
		if(classname == "GSNIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new GSNIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
