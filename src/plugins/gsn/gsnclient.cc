#include "soapA3DWebServiceSoap12BindingProxy.h"
#include "A3DWebServiceSoap12Binding.nsmap"

#include <iostream>
#include <sstream>
#include <ctime>

using namespace std;

int main(){

	A3DWebServiceSoap12BindingProxy gsn;
	_ns1__getSensorsResponse sensors;
	_ns1__getSensorLocation sensorloc_req;
	_ns1__getSensorLocationResponse sensorloc;
	_ns1__getSensorInfo sensorinfo_req;
	_ns1__getSensorInfoResponse sensorinfo;
	_ns1__getMeteoData meteodata_req;
	_ns1__getMeteoDataResponse meteodata;

	//_ns1__getSensorInfoResponse sensorinfo;

	/*
	gsn.proxy_host = "77.244.247.232"; // IP or domain
	gsn.proxy_port = 3128;
	gsn.proxy_userid = "username";
	gsn.proxy_passwd = "secret"; 
	*/
	double d1 = 1249550075500.0;
	double d2 = 1258887378500.0;

	stringstream ss1, ss2;
	ss1 << d1;
	ss2 << d2;

	meteodata_req.from = (LONG64)(d1);
	meteodata_req.to = (LONG64)(d2);

	cout << time(NULL) << "==" << d1 << "==" << meteodata_req.from <<endl;

	cout << "TEST Webservice" << endl;

	if (gsn.getSensors(&sensors) == SOAP_OK){
		cout << "Number of sensors accessible thorugh GSN: " << sensors.return_.size() << endl;
		for (unsigned int ii=0; ii<sensors.return_.size(); ii++){
			cout << "\tSensor " << ii << " Name: " << sensors.return_[ii] << endl;
		}

		for (unsigned int ii=0; ii<sensors.return_.size(); ii++){
			cout << "Sensor " << ii << " Name: " << sensors.return_[ii] << endl;
			//get Meta information
			sensorloc_req.sensor = new std::string(sensors.return_[ii]);
			sensorinfo_req.sensor = new std::string(sensors.return_[ii]);
			meteodata_req.sensor = new std::string(sensors.return_[ii]);

			if (gsn.getSensorLocation(&sensorloc_req, &sensorloc) == SOAP_OK){
				if (sensorloc.return_.size() == 3){
					cout << "\t" << sensorloc.return_[1] << endl;
					cout << "\t" << sensorloc.return_[2] << endl;
				}
			}
			
			if (gsn.getSensorInfo(&sensorinfo_req, &sensorinfo) == SOAP_OK){
				for (unsigned int jj=0; jj<sensorinfo.return_.size(); jj++){
					cout << "\t" << sensorinfo.return_[jj] << endl;					
				}
			}

			if (gsn.getMeteoData(&meteodata_req, &meteodata) == SOAP_OK){
				cout << "\t" << meteodata.return_.size() << endl;
				for (unsigned int jj=0; jj<meteodata.return_.size(); jj++){
					cout << "\t" << meteodata.return_[jj] << endl;					
				}
			} else {
				cout << "NOT OK" << endl;
			}

			delete sensorloc_req.sensor;
			delete sensorinfo_req.sensor;
			delete meteodata_req.sensor;
		}
	}


	return 0;
}
