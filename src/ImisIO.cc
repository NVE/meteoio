#include "ImisIO.h"

using namespace oracle;
using namespace oracle::occi;


ImisIO::ImisIO(const string& configfile) : IOInterface(NULL), cfg(configfile)
{
	mbImis.clear();
	vecStationName.clear();
}

ImisIO::ImisIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){}

ImisIO::~ImisIO() throw()
{
	cleanup();
}

void ImisIO::cleanup() throw()
{
	mbImis.clear();
	vecStationName.clear();
}

void ImisIO::createBuffer()
{
	//Clear the buffers
	mbImis.clear();
	const unsigned int stations = (unsigned int)vecStationName.size();

	//Allocate one MeteoBuffer per station
	for (unsigned int ii=0; ii<stations; ii++) {
		mbImis.push_back(MeteoBuffer(1000));
	}
	cout << "[I] "<<AT<<": Created Buffer for " << stations << " stations" << endl;
}

ConfigReader ImisIO::getCfg()
{
	return cfg;
}

vector<string> ImisIO::getVecStationName()
{
	return vecStationName;
}

vector<MeteoBuffer> ImisIO::getMbImis()
{
	return mbImis;
}

void ImisIO::get2DGridSize(int& nx, int& ny)
{
	//Nothing so far
	(void)nx;
	(void)ny;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::read2DGrid(Grid2DObject& grid_out, const string& parameter)
{
	//Nothing so far
	(void)parameter;
	(void)grid_out;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readDEM(Grid2DObject& dem_out)
{
	//Nothing so far
	(void)dem_out;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readLanduse(Grid2DObject& landuse_out)
{
	//Nothing so far
	(void)landuse_out;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	//Nothing so far
	(void)date_in;
	(void)da_out;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readSpecialPoints(CSpecialPTSArray& pts)
{
	//Nothing so far
	(void)pts;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::write2DGrid(const Grid2DObject& grid_in, const string& options)
{
	//Nothing so far
	(void)grid_in;
	(void)options;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo)
{
	vector<StationData> vecStation;
	readMeteoData(date_in, vecMeteo, vecStation);
}

void ImisIO::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
{	
	if (mbImis.size() == 0) {
		getStationName();
		//createBuffer();	
		setMbImis(date_in);
		resampleMbImis(vecMeteo, vecStation, date_in);
	} else {
		if (mbImis[0].getMeteoData(0).date <= date_in && mbImis[0].getMeteoData(mbImis[0].size()-1).date >= date_in) {
			resampleMbImis(vecMeteo, vecStation, date_in);
		} else {
			setMbImis(date_in);
			resampleMbImis(vecMeteo, vecStation, date_in);
		}
	}

	if (vecMeteo.size() == 0) {//No data found
		THROW IOException("[E] No data for any station for date " + date_in.toString() + " found", AT);
	}

}

void ImisIO::createData(vector< vector<string> >& meteo_in, vector<string>& station_in, MeteoBuffer& mb)
{
	MeteoData md;
	StationData sd;
	Date_IO tmpDate;
	
	double east, north, lat, lon, alt;
	convertString(east, station_in[1], dec);
	convertString(north, station_in[2], dec);
	convertString(alt, station_in[3], dec);
	string sName = "";
	if (station_in[0].size() == 0) {
		sName = vecStationName[mbImis.size()];
	} else {
		sName = station_in[0];
	}
	CH1903_to_WGS84(east, north, lat, lon);
	sd.setStationData(east, north, alt, sName, lat, lon);
	
	double ta, iswr, vw, dw, rh, lwr, nswc, tsg, tss, hs, rswr;
	const unsigned int size = meteo_in.size();
	for (unsigned int i=0; i<size; i++) {
		ImisIO::stringToDate(meteo_in[i][0], tmpDate);
		convertString(ta, meteo_in[i][1], dec);
		convertString(iswr, meteo_in[i][2], dec);
		convertString(vw, meteo_in[i][3], dec);
		convertString(dw, meteo_in[i][4], dec);
		convertString(rh, meteo_in[i][5], dec);
		convertString(lwr, meteo_in[i][5], dec);
		convertString(nswc, meteo_in[i][7], dec);
		convertString(tsg, meteo_in[i][8], dec);
		convertString(tss, meteo_in[i][9], dec);
		convertString(hs, meteo_in[i][10], dec);
		convertString(rswr, meteo_in[i][11], dec);
		md.setMeteoData(tmpDate, ta, iswr, vw, dw, rh, lwr, nswc, tsg, tss, hs, rswr);
		
		mb.put(md, sd);
	}
}

void ImisIO::getStation2Data(const string stat_abk, unsigned int stao_nr, vector<string>& data2S)
{
	const string userName = "slf";
	const string password = "sdb+4u";
	const string dbName = "sdbo";
	unsigned int timeOut = 0, seconds = 60;

	Environment *env = Environment::createEnvironment();// static OCCI function
	{
		Connection *conn;
		Statement *stmt;
		ResultSet *rs;
		while (timeOut != 3) {
			timeOut = 0;
			try {
				conn = env->createConnection(userName, password, dbName);
				timeOut++;
			} catch (SQLException &connex) {
				cout <<"getStation2Data : Connection failed, please verify if userName, password and dbName are correct........"<< endl;
				cout << connex.getMessage();
				exit(1);
			}
			try {
				stmt = conn->createStatement("select stao_name,stao_x,stao_y,stao_h from station2.standort                              								where STAT_ABK =: 1 AND STAO_NR =: 2");
				stmt->setString(1, stat_abk); // set 1st variable's value
				stmt->setInt(2, stao_nr); // set 2nd variable's value 		
				rs = stmt->executeQuery(); // execute the statement stmt
				timeOut++;
			} catch (SQLException &stmtex) {
				cout <<"getStation2Data : Statement failed, please verify if it is correctly written............"<< endl;
				cout << stmtex.getMessage();
				exit(1);
			}
			try {			
				while (rs->next() == true) {
					for (int i=0; i<4; i++) {
						data2S.push_back(rs->getString(i+1));
					}
				}
				timeOut++;
			} catch (SQLException &rsex) {
				cout <<"getStation2Data : ResultSet manipulation failed, please verify if there is no mistake............."<< endl;
				cout << rsex.getMessage();
				exit(1);
			}catch (exception &cppex) { // C++ exception
				cout<< "Error "<< cppex.what()<<endl;
			}
			if (timeOut != 3 && seconds <= 27*60) {
				sleep(seconds);
				seconds *= 3;
			} else if (seconds > 27*60) {
				break;
			}	
		}   	   
		stmt->closeResultSet(rs);
		conn->terminateStatement(stmt);
		env->terminateConnection(conn);
	}
	Environment::terminateEnvironment(env); // static OCCI function
		
}

void ImisIO::getImisData (const string &stat_abk, const unsigned int &stao_nr, vector<int> date_in, vector< vector<string> >& dataImis)
{
	const string userName = "slf";
	const string password = "sdb+4u";
	const string dbName = "sdbo"; //or sdbo
	vector<string> vec;
	unsigned int timeOut = 0, seconds = 60;

	Environment *env = Environment::createEnvironment();// static OCCI function
	{
		Connection *conn;
		Statement *stmt;
		ResultSet *rs;
		while (timeOut != 3) {
			timeOut = 0;
			try {
				conn = env->createConnection(userName, password, dbName);
				timeOut++;
			} catch (SQLException &connex) {
				cout <<"getImisData : Connection failed, please verify if userName, password and dbName are correct........."<< endl;
				cout << connex.getMessage();
				exit(1);
			}
			try {
				if (stao_nr != 0) {
					stmt = conn->createStatement("select to_char(datum, 'YYYY/MM/DD HH24:MI') as datum,ta,iswr,vw,dw,rh,lwr,nswc, 									    tsg,tss,hs,rswr from ams.v_amsio where STAT_ABK =: 1 AND STAO_NR =: 2             									    and DATUM >=: 3 and rownum<=1000");
					Date edate(env, date_in[0], date_in[1], date_in[2], date_in[3], date_in[4]); // year, month, day, hour, minutes
					stmt->setString(1, stat_abk); // set 1st variable's value
					stmt->setInt(2, stao_nr); // set 2nd variable's value
					stmt->setDate(3, edate); // set 3rd variable's value
				} else {
					string sql = "select to_char(datum, 'YYYY/MM/DD HH24:MI') as datum,ta,iswr,vw,dw,rh,lwr,nswc,tsg,tss,hs,rswr 								from ams.v_amsio where STAT_ABK=:1 AND STAO_NR is null and DATUM>=:2 and rownum<=1000";
					stmt = conn->createStatement(sql);
					Date edate(env, date_in[0], date_in[1], date_in[2], date_in[3], date_in[4]); // year, month, day, hour, minutes
					stmt->setString(1, stat_abk); // set 1st variable's value
					stmt->setDate(2, edate); // set 2nd variable's value				
				}
				rs = stmt->executeQuery(); // execute the statement stmt
				timeOut++;
			} catch (SQLException &stmtex) {
				cout <<"getImisData : Statement failed, please verify if it is correctly written............"<< endl;
				cout << stmtex.getMessage();
				exit(1);
			}
			try {			
				rs->setMaxColumnSize(7,22);
				while (rs->next() == true) {
					vec.clear();
					for (int i=1; i<=12; i++) { // 12 columns 
						vec.push_back(rs->getString(i));
					}
					dataImis.push_back(vec);
				}
				timeOut++;
			} catch (SQLException &rsex) {
				cout <<"getImisData : ResultSet manipulation failed, please verify if there is no mistake............."<< endl;
				cout << rsex.getMessage();
				exit(1);
			}catch (exception &cppex) { // C++ exception
				cout<< "Error "<< cppex.what()<<endl;
			}
			if (timeOut != 3 && seconds <= 27*60) {
				sleep(seconds);
				seconds *= 3;
			} else if (seconds > 27*60) {
				break;
			}
		}   	   
		stmt->closeResultSet(rs);
		conn->terminateStatement(stmt);
		env->terminateConnection(conn);
	}
	Environment::terminateEnvironment(env); // static OCCI function
}

void ImisIO::getStationName()
{
	vecStationName.clear();
	
	string str_stations="";
	int stations=0;
	
	cfg.getValue("NROFSTATIONS", str_stations);
	
	if (!convertString(stations, str_stations, dec)) {
		THROW ConversionFailedException("Error while reading value for NROFSTATIONS", AT);
	}
	for (int ii=0; ii<stations; ii++) {
		stringstream tmp_stream;
		string stationname="";
		
		tmp_stream << (ii+1); //needed to construct key name
		cfg.getValue(string("STATION"+tmp_stream.str()), stationname);
		
		vecStationName.push_back(stationname);
	}
}
		
void ImisIO::setMbImis(Date_IO date_in)
{
	mbImis.clear();
	int date[5], stao;
	date_in -= 1./24.; // one hour before
	date_in.getDate(date[0],date[1],date[2],date[3],date[4]);
	vector<int> date_io(date, date+sizeof(date)/sizeof(int));
	string station, name, number;
	const unsigned int size = vecStationName.size();
	
	for (unsigned int i=0; i<size; i++) {
		vector<string>* data2s = new vector<string>;
		vector< vector<string> >* data_imis = new vector< vector<string> >;
		MeteoBuffer* mb = new MeteoBuffer(1000);
		station = vecStationName[i];
		number = station.substr(station.length()-1);
		name = station.substr(0, station.length()-1);
	
		if (!convertString(stao, number, dec)) {
			THROW ConversionFailedException("Error while reading station number in readMeteoData(...) ", AT);
		}
	
		getStation2Data(name,stao,*data2s);
		getImisData(name,stao,date_io,*data_imis);
		createData(*data_imis,*data2s,*mb);
		mbImis.push_back(*mb);
		free(data2s);
		free(data_imis);
		free(mb);
	}
}

void ImisIO::resampleMbImis(vector<MeteoData>& vecMeteo, vector<StationData>& vecStation, const Date_IO& date_in)
{
	vecMeteo.clear();
	vecStation.clear();
	const unsigned int size = mbImis.size();
	
	for (unsigned int ii=0; ii<size; ii++) {
		unsigned int index = mbImis[ii].seek(date_in);
		if (index != MeteoBuffer::npos) {
			if (mbImis[ii].getMeteoData(index).date == date_in) {
				vecMeteo.push_back(mbImis[ii].getMeteoData(index));
				vecStation.push_back(mbImis[ii].getStationData(index));
			} else {
				Meteo1DResampler mresampler;
				mresampler.resample(index, date_in, mbImis[ii]);
				if (index != MeteoBuffer::npos) {
					vecMeteo.push_back(mbImis[ii].getMeteoData(index));
					vecStation.push_back(mbImis[ii].getStationData(index));
				} else {
					cout << "[i] Buffering data for Station " << vecStationName[ii] << " at date " 
						<< date_in.toString() << " failed" << endl;
				}
			}
		}
	}
}

void ImisIO::stringToDate(const string& instr, Date_IO& date_out)
{
	int tmp[5];
	
	string year = instr.substr(0,4);
	string month = instr.substr(5,2);
	string day = instr.substr(8,2);	
	string hour = instr.substr(11,2);
	string min = instr.substr(14,2);
	
	convertString(tmp[0], year, dec);
	convertString(tmp[1], month, dec);
	convertString(tmp[2], day, dec);
	convertString(tmp[3], hour, dec);
	convertString(tmp[4], min, dec);
	
	date_out.setDate(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4]);
}

double ImisIO::strToDouble(const string &str)
{ //HACK: this method should be replaced by a call to a standard one
	int length = str.size();
	double result;
	if (length == 0) {
		return nodata;
	} else {
		istringstream ss(str);
		ss >> result;		
		return result;
	}
}

extern "C"
{
	//using namespace MeteoIO;
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}
	
	void* loadObject(const string& classname, const string& filename) {
		if(classname == "ImisIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new ImisIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
