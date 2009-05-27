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

ConfigReader ImisIO::getCfg()
{
	return cfg;
}

vector<string> ImisIO::getVecStationName()
{
	return vecStationName;
}

void ImisIO::get2DGridSize(int& nx, int& ny)
{
	//Nothing so far
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::read2DGrid(Grid2DObject& grid_out, const string& parameter)
{
	//Nothing so far
	(void)parameter;
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readDEM(Grid2DObject& dem_out)
{
	//Nothing so far
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readLanduse(Grid2DObject& landuse_out)
{
	//Nothing so far
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	//Nothing so far
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::readSpecialPoints(CSpecialPTSArray& pts)
{
	//Nothing so far
	THROW IOException("Nothing implemented here", AT);
}

void ImisIO::write2DGrid(const Grid2DObject& grid_in, const string& options)
{
	//Nothing so far
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
	vecMeteo.clear();
	vecStation.clear();
	int date[5];
	date_in.getDate_IO(date[0],date[1],date[2],date[3],date[4]);
	vector<int> date_io(date, date+sizeof(date)/sizeof(int));
	
	getStationName();
	for (unsigned int i=0; i< vecStationName.size(); i++) {
		vector<string> data2s;
		vector< vector<string> > data_imis;
		string station = vecStationName[i].substr(0,3);
		int stao;
		if (!convertString(stao, vecStationName[i].substr(3), dec)) {
			THROW ConversionFailedException("Error while reading station number in readMeteoData(...) ", AT);
		}
		getStation2Data(station,stao,data2s);
		getImisData(station,stao,date_io,data_imis);		
		MeteoBuffer mb(1000);
		createData(vecMeteo,vecStation,data_imis,data2s,mb);
		mbImis.push_back(mb);
	}	
}

void ImisIO::createData(vector<MeteoData>& vecMeteo, vector<StationData>& vecStation,
			vector< vector<string> >& meteo_in, vector<string>& station_in, MeteoBuffer& mb)
{
	MeteoData md;
	StationData sd;
	Date_IO tmpDate;
	
	double east, north, lat, lon, alt;
	convertString(lat, station_in[1], dec);
	convertString(lon, station_in[2], dec);
	convertString(alt, station_in[3], dec);
	string sName = "";
	if (station_in[0].size() == 0) {
		sName = "Unnamed Station";
	} else {
		sName = station_in[0];
	}
	WGS84_to_CH1903(lat, lon, east, north);
	sd.setStationData(east, north, alt, sName, lat, lon);
	vecStation.push_back(sd);
	
	double ta, iswr, vw, rh, lwr, nswc, ts0, hs, rswr;
	for (unsigned int i=0; i<meteo_in.size(); i++) {
		ImisIO::stringToDate(meteo_in[i][0], tmpDate);
		/*ta = strToDouble(meteo_in[i][1]);*/convertString(ta, meteo_in[i][1], dec);
		/*iswr = strToDouble(meteo_in[i][2]);*/convertString(iswr, meteo_in[i][2], dec);
		/*vw = strToDouble(meteo_in[i][3]);*/convertString(vw, meteo_in[i][3], dec);
		/*rh = strToDouble(meteo_in[i][4]);*/convertString(rh, meteo_in[i][4], dec);
		/*lwr = strToDouble(meteo_in[i][5]);*/convertString(lwr, meteo_in[i][5], dec);
		/*nswc = strToDouble(meteo_in[i][6]);*/convertString(nswc, meteo_in[i][6], dec);
		/*ts0 = strToDouble(meteo_in[i][7]);*/convertString(ts0, meteo_in[i][7], dec);
		/*hs = strToDouble(meteo_in[i][8]);*/convertString(hs, meteo_in[i][8], dec);
		/*rswr = strToDouble(meteo_in[i][9]);*/convertString(rswr, meteo_in[i][9], dec);
		md.setMeteoData(tmpDate, ta, iswr, vw, rh, lwr, nswc, ts0, hs, rswr);
		vecMeteo.push_back(md);
		mb.put(md, sd);
	}	
}

void ImisIO::getStation2Data(const string stat_abk, unsigned int stao_nr, vector<string>& data2S)
{
	const string userName = "slf";
	const string password = "sdb+4u";
	const string dbName = "sdbt";
	int timeOut = 0, seconds = 60;

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
				cout <<"Connection failed, please verify if userName, password and dbName are correct........."<< endl;
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
				cout <<"Statement failed, please verify if it is correctly written............"<< endl;
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
				cout <<"ResultSet manipulation failed, please verify if there is no mistake............."<< endl;
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
	const string dbName = "sdbt"; //or sdbo
	vector<string> vec;
	int timeOut = 0, seconds = 60;

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
				cout <<"Connection failed, please verify if userName, password and dbName are correct........."<< endl;
				cout << connex.getMessage();
				exit(1);
			}
			try {
				stmt = conn->createStatement("select to_char(datum, 'YYYY/MM/DD HH24:MI') as datum,ta,iswr,vw,rh,lwr,nswc,ts0,hs,rswr   						     	     from ams.v_amsio where STAT_ABK =: 1 AND STAO_NR =: 2 and DATUM >=: 3 and rownum<=100");
				Date edate(env, date_in[0], date_in[1], date_in[2], date_in[3], date_in[4]); // year, month, day, hour, minutes
				stmt->setString(1, stat_abk); // set 1st variable's value
				stmt->setInt(2, stao_nr); // set 2nd variable's value 		   
				stmt->setDate(3, edate); // set 3rd variable's value
				rs = stmt->executeQuery(); // execute the statement stmt
				timeOut++;
			} catch (SQLException &stmtex) {
				cout <<"Statement failed, please verify if it is correctly written............"<< endl;
				cout << stmtex.getMessage();
				exit(1);
			}
			try {			
				rs->setMaxColumnSize(6,22);
				while (rs->next() == true) {
					vec.clear();
					for (int i=1; i<=10; i++) { // 10 columns 
						vec.push_back(rs->getString(i));
					}
					dataImis.push_back(vec);
				}
				timeOut++;
			} catch (SQLException &rsex) {
				cout <<"ResultSet manipulation failed, please verify if there is no mistake............."<< endl;
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
	
	date_out.setDate_IO(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4]);
}

double ImisIO::strToDouble(const string &str)
{
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

void ImisIO::displayData()
{
	cout<<endl <<"Contenu de mbImis : " <<endl;
	cout<<"----------------------------------------------------------------------------------------------------------" <<endl;
	cout<<" N° |    Station    |          Date          | ta  | iswr |  vw |  rh  | lwr  | nswc | ts0  |  hs  | rswr |" <<endl;
	cout<<"----------------------------------------------------------------------------------------------------------" <<endl;
	int rows = 0;
	for (unsigned int i=0; i<mbImis.size(); i++) {
		deque<MeteoData> mbMet = mbImis[i].getMeteobuffer();
		deque<StationData> mbSta = mbImis[i].getStationbuffer();
		for (unsigned int ii=0; ii<mbMet.size(); ii++) {
			if (rows<9) {
				cout<<"00" <<rows+1 <<" | ";
			} else if (rows<99){
				cout<<"0" <<rows+1 <<" | ";
			} else {
				cout<<rows+1 <<" | ";
			}
			cout<<mbSta[ii].stationName <<" | ";
			cout<<mbMet[ii].date <<" | ";
			cout<<mbMet[ii].ta <<" | ";
			cout<<mbMet[ii].iswr <<" | ";
			cout<<mbMet[ii].vw <<" | ";
			cout<<mbMet[ii].rh <<" | ";
			cout<<mbMet[ii].lwr <<" | ";
			cout<<mbMet[ii].nswc <<" | ";
			cout<<mbMet[ii].ts0 <<" | ";
			cout<<mbMet[ii].hs <<" | ";
			cout<<mbMet[ii].rswr <<" | ";
			cout<<endl;
			rows++;
		}
	}
	cout<<"----------------------------------------------------------------------------------------------------------" <<endl;
	cout<<"rows : " <<rows <<endl;
}

void ImisIO::test(vector<int> date) {
	vector<MeteoData> vecMeteo;
	vector<StationData> vecStation;
		
	Date_IO date_in(date[0],date[1],date[2],date[3],date[4]);
	readMeteoData(date_in, vecMeteo);
	displayData();

}

int main(int argc, char** argv) {
	if(argc<6) {
		cout<< " Pas assez d'arguments, tu t'es planté mon pote ahahahahahahaha........:-)E) "<< endl;
		exit(1);
	} else {
		vector<int> date;
		for (int i=0; i<argc-1; i++) {
			date.push_back(atoi(argv[i+1]));
		}
		ImisIO imis("io.ini");
		
		imis.test(date);
	}
	return EXIT_SUCCESS;	
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
