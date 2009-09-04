#include "BufferedIOHandler.h"

using namespace std;

#ifdef _POPC_
BufferedIOHandler::BufferedIOHandler(IOInterface& _iohandler, const ConfigReader& _cfg) 
	: iohandler(_iohandler), cfg(_cfg), meteoBuffer(), stationBuffer(), mapBufferedGrids(), resample(true)
#else
BufferedIOHandler::BufferedIOHandler(IOInterface& _iohandler, const ConfigReader& _cfg) 
	: IOInterface(NULL), iohandler(_iohandler), cfg(_cfg), meteoBuffer(), stationBuffer(), mapBufferedGrids(), resample(true)
#endif
{
	//Nothing else so far
}

#ifdef _POPC_
BufferedIOHandler::~BufferedIOHandler()
#else
BufferedIOHandler::~BufferedIOHandler() throw()
#endif
{
	//Nothing else so far
}

void BufferedIOHandler::read2DGrid(Grid2DObject& _grid2Dobj, const string& _filename)
{
	map<string, Grid2DObject>::iterator it = mapBufferedGrids.find(_filename);
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.read2DGrid(tmpgrid2D, _filename);
	mapBufferedGrids[_filename] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readDEM(Grid2DObject& _grid2Dobj)
{
	map<string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:DEM");
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readDEM(tmpgrid2D);
	mapBufferedGrids["/:DEM"] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readLanduse(Grid2DObject& _grid2Dobj)
{
	map<string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:LANDUSE");
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readLanduse(tmpgrid2D);
	mapBufferedGrids["/:LANDUSE"] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readAssimilationData(const Date_IO& _date, Grid2DObject& _grid2Dobj)
{
	map<string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:ASSIMILATIONDATA" + _date.toString());
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readAssimilationData(_date, tmpgrid2D);
	mapBufferedGrids["/:ASSIMILATIONDATA" + _date.toString()] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, vector< vector<MeteoData> >& vecMeteo)
{
	vector< vector<StationData> > vecStation;
	readMeteoData(dateStart, dateEnd, vecMeteo, vecStation);
}

void BufferedIOHandler::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation){
	/* For every station: 
	 * 1) See whether data is already buffered
	 * 2) Filter - not implemented yet
	 * 3) Resample if necessary
	 * 4) Filter resampled value
	 * 5) Return the values
	 */

	vecMeteo.clear();
	vecStation.clear();

	if (meteoBuffer.size() == 0){ //init
		bufferAllData(date_in);
	}

	//loop through all meteo buffers
	for (unsigned int ii=0; ii<meteoBuffer.size(); ii++) {
		unsigned int index = BufferedIOHandler::npos;

		StationData sd;
		string stationName = "";
		if (stationBuffer[ii].size() > 0){
			stationName = stationBuffer[ii][0].stationName;
			sd = stationBuffer[ii][0];
		}

		if (meteoBuffer[ii].size() > 0) {//check whether meteo data for the date exists in buffer
			index = seek(date_in, meteoBuffer[ii]);
		}

		if (index == BufferedIOHandler::npos) { //not in buffer
			cout << "[I] Station " << ii << "(" << stationName 
				<< ") data for date " << date_in.toString() << " not in buffer ..." << endl;

			bool dataexists = bufferData(date_in, ii);
			if (dataexists) {//date_in is contained in buffer
				index = seek(date_in, meteoBuffer[ii]);
			}
		} /*else {
			cout << "[I] Found data for station " << stationName << " and date " << date_in.toString() << " in buffer" << endl;
		}*/
    
		// RESAMPLING
		if ((index != BufferedIOHandler::npos) && (meteoBuffer[ii][index].date != date_in)) {
			if (resample){
				cerr << "[I] Resampling required for date: " << date_in.toString() << endl;
				
				Meteo1DResampler mresampler;
				mresampler.resample(index, date_in, meteoBuffer[ii], stationBuffer[ii]);
			} else {
				//insert a nodata touple
				meteoBuffer[ii].insert(meteoBuffer[ii].begin()+index, MeteoData());
				meteoBuffer[ii][index].date = date_in; //set correct date
				stationBuffer[ii].insert(stationBuffer[ii].begin()+index, sd);
			}
		}
    
		if (index != BufferedIOHandler::npos) {
			vecMeteo.push_back(meteoBuffer[ii][index]);
			vecStation.push_back(stationBuffer[ii][index]);
		} else {
			cout << "[I] Buffering data for Station " << stationName << " at date " 
				<< date_in.toString() << " failed" << endl;
			vecMeteo.push_back(MeteoData());
			vecMeteo[ii].date = date_in; //set correct date

			vecStation.push_back(StationData());
		}
	}

	if (vecMeteo.size() == 0) {//No data found - return one object set to date_in and nodata in all other fields
		vecMeteo.push_back(MeteoData());
		vecMeteo[0].date = date_in; //set correct date

		vecStation.push_back(StationData());
		//throw IOException("[E] No data for any station for date " + date_in.toString() + " found", AT);
	}	
}

void BufferedIOHandler::getNextMeteoData(const Date_IO& _date, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
{
	//TODO: check whether there is something in the buffer!
	//Try to rebuffer!

	vecMeteo.clear();
	vecStation.clear();
	
	vector< vector<MeteoData> > meteoTmpBuffer;
	vector< vector<StationData> > stationTmpBuffer;
	readMeteoData(_date, (_date-Date_IO(1900,1,2)), meteoTmpBuffer, stationTmpBuffer);	

	unsigned int emptycounter = 0;
	for (unsigned int ii=0; ii<meteoTmpBuffer.size(); ii++){//stations
		if ((meteoTmpBuffer[ii].size() > 0) && (stationTmpBuffer[ii].size() > 0)){
			vecMeteo.push_back(meteoTmpBuffer[ii][0]);
			vecStation.push_back(stationTmpBuffer[ii][0]);
		} else {
			emptycounter++;
		}
	}

	if (emptycounter == meteoTmpBuffer.size()){
		vecMeteo.clear();
		vecStation.clear();
	}
}

void BufferedIOHandler::bufferAllData(const Date_IO& _date){
	Date_IO fromDate = _date - (Date_IO(1900,1,5));  //minus 5 days
	Date_IO toDate   = _date + (Date_IO(1900,1,20)); //plus 20 days

	readMeteoData(fromDate, toDate, meteoBuffer, stationBuffer);
}

bool BufferedIOHandler::bufferData(const Date_IO& _date, const unsigned int& stationindex)
{
	Date_IO fromDate = _date - (Date_IO(1900,1,5));  //minus 5 days
	Date_IO toDate   = _date + (Date_IO(1900,1,20)); //plus 20 days

	readMeteoData(fromDate, toDate, meteoBuffer, stationBuffer, stationindex);

	if (meteoBuffer.size() == 0) {
		return false;
	}
	
	if (meteoBuffer[stationindex].size() == 0) {
		return false;
	}

	if ((!((_date >= meteoBuffer[stationindex][0].date) 
		  && (meteoBuffer[stationindex][meteoBuffer[stationindex].size()-1].date >= _date)))) {
		meteoBuffer[stationindex].clear();
		return false;
	}

	//If we reach this point: Date_IO is definitely covered
	return true;
}


void BufferedIOHandler::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							   std::vector< std::vector<MeteoData> >& vecMeteo, 
							   std::vector< std::vector<StationData> >& vecStation,
							   const unsigned int& stationindex)
	
{
	iohandler.readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);
	if ((&meteoBuffer != &vecMeteo) && (&stationBuffer != &vecStation)){
		meteoBuffer = vecMeteo;      //copy by value
		stationBuffer = vecStation;  //copy by value
	}
}

void BufferedIOHandler::readSpecialPoints(CSpecialPTSArray& _cpa)
{
	iohandler.readSpecialPoints(_cpa);
}

void BufferedIOHandler::write2DGrid(const Grid2DObject& _grid2Dobj, const string& _name)
{
	iohandler.write2DGrid(_grid2Dobj, _name);
}

void BufferedIOHandler::enableResampling(const bool& _resample)
{
	resample = _resample;
}

void BufferedIOHandler::clearBuffer(){
	meteoBuffer.clear();
	stationBuffer.clear();
	mapBufferedGrids.clear();
}

unsigned int BufferedIOHandler::seek(const Date_IO& date_in, vector<MeteoData>& vecM){ //TODO: binary search
	//returns index of element, if element does not exist it returns closest index
	//the element needs to be an exact hit or embedded between two measurments

	unsigned int ii = 1;

	if (vecM.size() <= 0) {//no elements in buffer
		return BufferedIOHandler::npos;
	}

	//if we reach this point: at least one element in buffer
	if (vecM[0].date > date_in) {
		return BufferedIOHandler::npos;
	}

	if (vecM[vecM.size()-1].date < date_in) {//last element is earlier, return npos
		return BufferedIOHandler::npos;
	}

	if (vecM[0].date == date_in) {//closest element
		return 0;
	}

	//if we reach this point: the date is spanned by the buffer and there are at least two elements
	while ((ii < vecM.size())) {
		//cerr << "in search loop" << vecM[ii].date.toString() << endl;
		if ((vecM[ii].date >= date_in) && (vecM[ii-1].date < date_in)) {
			return ii;
		}

		ii++;
	}

	return BufferedIOHandler::npos;
}
