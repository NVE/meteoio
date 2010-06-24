#include <fstream>
#include <iostream>
#include <string>
#include "MeteoIO.h"

using namespace std;
using namespace IOUtils;


void create1DFile(const vector< vector<MeteoData> >& data, const vector<StationData>& stations)
{
	unsigned int sta_nr = stations.size(), size = data.size();
	for(unsigned int ii=0; ii<sta_nr; ii++) {
		string file = "meteo1D_"+stations[ii].getStationID()+".txt";
		ofstream flux(file.c_str(), ios::out | ios::trunc);
		if (flux) {
			flux<<"Name = " <<stations[ii].getStationID() <<endl;
			flux<<"Latitude = " <<stations[ii].latitude <<endl;
			flux<<"Longitude = " <<stations[ii].longitude <<endl;
			flux<<"X_Coord = " <<stations[ii].eastCoordinate <<endl;
			flux<<"Y_Coord = " <<stations[ii].northCoordinate <<endl;
			flux<<"Altitude = " <<stations[ii].altitude <<endl;
			flux<<"YYYY MM DD HH ta iswr vw rh ea nswc" <<endl;
			for(unsigned int j=0; j<size; j++) {
				int yyyy, mm, dd, hh;
				if (data[j][ii].iswr!=nodata || data[j][ii].lwr!=nodata || data[j][ii].nswc!=nodata) {
					data[j][ii].date.getDate(yyyy, mm, dd, hh);
					flux<<yyyy <<" " <<mm <<" " <<dd <<" " <<hh <<" ";
					flux<<data[j][ii].ta <<" ";
					flux<<data[j][ii].iswr <<" ";
					flux<<data[j][ii].vw <<" ";
					flux<<data[j][ii].rh <<" ";
					flux<<data[j][ii].lwr <<" ";
					flux<<data[j][ii].nswc <<endl;
				}
			}
			flux.close();
		} else { 
			cerr << "Erreur à l'ouverture !" << endl;
	        }
	}
}

void writeHeader(ofstream &file, const vector<StationData>& stations, const string parameter_name)
{
	ostringstream str_altitudes;
	ostringstream str_eastings;
	ostringstream str_northings;
	unsigned int sta_nr = stations.size();
	
	file<<"X:\\filepath " << parameter_name <<endl;
	for(unsigned int ii=0;ii<sta_nr;ii++) {
		str_altitudes << stations[ii].altitude << " ";
		str_eastings << stations[ii].eastCoordinate << " ";
		str_northings << stations[ii].northCoordinate << " ";
	}
	file<<"YY MM DD HH "<< str_altitudes.str() <<endl; //altitudes
	file<<"YY MM DD HH "<< str_eastings.str() <<endl; //easting
	file<<"YY MM DD HH "<< str_northings.str() <<endl; //northing
	file<<"YYYY MM DD HH";
	for(unsigned int i=0; i<sta_nr; i++) {
		file<<" " <<stations[i].stationID.substr(0,3);
	}
	file<<endl;
}

void create2DFile(const string& type, const vector< vector<MeteoData> >& data, const vector<StationData>& stations)
{
	int year, month, day, hour;
	unsigned int sta_nr = stations.size(), size = data.size();
	data[0][0].date.getDate(year, month, day);
	ostringstream out;
	out<<year; 
	if (type == "nswc") {
		string name = "prec"+out.str()+".txt";
		ofstream file(name.c_str(), ios::out | ios::trunc);
		if (file) {
			writeHeader(file, stations, "precipitations");
			for(unsigned int ii=0; ii<size; ii++) {
				data[ii][0].date.getDate(year, month, day, hour);
				if (month<10) {
					if (day<10) {
						if (hour<10) {
							file<<year <<" 0" <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" 0" <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" " <<day <<" " <<hour;
						}
					}
				} else {
					if (day<10) {
						if (hour<10) {
							file<<year <<" " <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" " <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" " <<day <<" " <<hour;
						}
					}
				}
				for(unsigned int j=0; j<sta_nr; j++) {
					file<<" " <<data[ii][j].nswc;
				}
				file<<endl;
			}
			file.close();
		} else { 
			cerr << "Erreur à l'ouverture !" << endl;
	        }
	        	        
	} else if (type == "rh") {
		string name = "rhum"+out.str()+".txt";
		ofstream file(name.c_str(), ios::out | ios::trunc);
		if (file) {
			writeHeader(file, stations, "relative humidity");
			
			for(unsigned int ii=0; ii<size; ii++) {
				data[ii][0].date.getDate(year, month, day, hour);
				if (month<10) {
					if (day<10) {
						if (hour<10) {
							file<<year <<" 0" <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" 0" <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" " <<day <<" " <<hour;
						}
					}
				} else {
					if (day<10) {
						if (hour<10) {
							file<<year <<" " <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" " <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" " <<day <<" " <<hour;
						}
					}
				}
				for(unsigned int j=0; j<sta_nr; j++) {
					file<<" " <<data[ii][j].rh;
				}
				file<<endl;
			}
			file.close();
		} else { 
			cerr << "Erreur à l'ouverture !" << endl;
	        }
	
	} else if (type == "ta") {
		string name = "tair"+out.str()+".txt";
		ofstream file(name.c_str(), ios::out | ios::trunc);
		if (file) {
			writeHeader(file, stations, "air temperature");
			
			for(unsigned int ii=0; ii<size; ii++) {
				data[ii][0].date.getDate(year, month, day, hour);
				if (month<10) {
					if (day<10) {
						if (hour<10) {
							file<<year <<" 0" <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" 0" <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" " <<day <<" " <<hour;
						}
					}
				} else {
					if (day<10) {
						if (hour<10) {
							file<<year <<" " <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" " <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" " <<day <<" " <<hour;
						}
					}
				}
				for(unsigned int j=0; j<sta_nr; j++) {
					file<<" " <<data[ii][j].ta;
				}
				file<<endl;
			}
			file.close();
		} else { 
			cerr << "Erreur à l'ouverture !" << endl;
	        }
	
	} else if (type == "vw") {
		string name = "wspd"+out.str()+".txt";
		ofstream file(name.c_str(), ios::out | ios::trunc);
		if (file) {
			writeHeader(file, stations, "wind velocity");
			
			for(unsigned int ii=0; ii<size; ii++) {
				data[ii][0].date.getDate(year, month, day, hour);
				if (month<10) {
					if (day<10) {
						if (hour<10) {
							file<<year <<" 0" <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" 0" <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" " <<day <<" " <<hour;
						}
					}
				} else {
					if (day<10) {
						if (hour<10) {
							file<<year <<" " <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" " <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" " <<day <<" " <<hour;
						}
					}
				}
				for(unsigned int j=0; j<sta_nr; j++) {
					file<<" " <<data[ii][j].vw;
				}
				file<<endl;
			}
			file.close();
		} else { 
			cerr << "Erreur à l'ouverture !" << endl;
	        }
	
	} else if (type == "dw") {
		string name = "wdir"+out.str()+".txt";
		ofstream file(name.c_str(), ios::out | ios::trunc);
		if (file) {
			writeHeader(file, stations, "wind direction");
			
			for(unsigned int ii=0; ii<size; ii++) {
				data[ii][0].date.getDate(year, month, day, hour);
				if (month<10) {
					if (day<10) {
						if (hour<10) {
							file<<year <<" 0" <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" 0" <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" 0" <<month <<" " <<day <<" " <<hour;
						}
					}
				} else {
					if (day<10) {
						if (hour<10) {
							file<<year <<" " <<month <<" 0" <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" 0" <<day <<" " <<hour;
						}
					} else {
						if (hour<10) {
							file<<year <<" " <<month <<" " <<day <<" 0" <<hour;
						} else {
							file<<year <<" " <<month <<" " <<day <<" " <<hour;
						}
					}
				}
				for(unsigned int j=0; j<sta_nr; j++) {
					file<<" " <<data[ii][j].dw;
				}
				file<<endl;
			}
			file.close();
		} else {
			cerr << "Erreur à l'ouverture !" << endl;
	        }
	}
}

int main(int argc, char *argv[])
{
	if (argc!=3) {
		cout<<"Two (2) date are needed in this format : YYYY-MM-DDTHH:mm:ss" <<endl;
		exit(1);
	} else {
		Date d1;
		Date d2;
		convertString(d1,string(argv[1]));
		convertString(d2,string(argv[2]));
	
		vector<MeteoData> vecMeteo;
		vector<StationData> vecStation;
		vector<StationData> tmpStation;
		vector< vector<MeteoData> > data;
		unsigned int k = 0;
	
		IOInterface *ioTest=NULL; //Initialization vital!
	
		try {
			ioTest = new IOHandler("io.ini");
		} catch (exception& e){
			cout << "Problem with IOHandler creation, cause: " << e.what() << endl;
		}
	
		try {
			while (d1 < d2) {
				ioTest->readMeteoData(d1, vecMeteo, tmpStation);
				data.push_back(vecMeteo);
				if (tmpStation.size() > k) {
					k = tmpStation.size();
					vecStation = tmpStation;
				}
				d1 += 1./24.; // incremented by 1 hour
			}
		} catch (exception& e){
			cout << "Problem when reading data, cause: " << e.what() << endl;
		}
		
		create1DFile(data, vecStation);
		string type[] = {"nswc","rh","ta","vw","dw"};
		for(int i=0; i<5; i++) {
			create2DFile(type[i], data, vecStation);
		}
	}
	
	return 0;
}



//compile with: g++ createFiles.cc ../src/Laws.c -I ../src/ -I ../src/filter/ -L ../lib/ -lmeteoio -ldl -lm -o create_demo -rdynamic

