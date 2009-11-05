#include "FilterAlgorithms.h"

using namespace std;

std::map<std::string, FilterProperties> FilterAlgorithms::filterMap;
const bool FilterAlgorithms::__init = FilterAlgorithms::initStaticData();

bool FilterAlgorithms::initStaticData()
{
	filterMap["rate"]     = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::RateFilter);
	filterMap["resample"] = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::ResamplingFilter);
	filterMap["min_max"]  = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MinMaxFilter);
	filterMap["min"]      = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MinValueFilter);
	filterMap["max"]      = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MaxValueFilter);

	return true;
}

const FilterProperties& FilterAlgorithms::filterProperties(const string& filtername)
{
	std::map<string, FilterProperties>::const_iterator it;
	it = filterMap.find(filtername);

	if (it==filterMap.end())
		throw UnknownValueException("Unknown filter called: " + filtername, AT);

	return it->second;
}


/******************************************************************************
 * The following functions are implementations of different filter algorithms *
 ******************************************************************************/

bool FilterAlgorithms::RateFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
						    const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
						    const unsigned int& paramindex,
						    vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS)
{
	(void)date; (void)vecS; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	vector<double> vecArgs; 
	parseFilterArguments("rate", _vecArgs, 1, 1, isSoft, vecArgs);
	
	Date_IO deltatime(vecArgs[0]/1440.0);

	unsigned int startposition = pos;
	for (unsigned int ii=vecFilteredM.size()-1; ii>0; ii--){
		unsigned int mypos = startposition;
		
		double sum = 0.0;
		while((vecM[mypos].date + deltatime) >= vecM[startposition].date){
			const double& val = vecM[mypos].param(paramindex);
			if (val != IOUtils::nodata)
				sum += vecM[mypos].param(paramindex);

			if (mypos>0) mypos--;
			else break;
		}
		
		//cout << "Accumulation for element (" << ii << "): " << sum << endl;
		vecFilteredM[ii].param(paramindex) = sum;

		if (startposition>0) startposition--; 
		else break;
	}

	return true;
}

bool FilterAlgorithms::ResamplingFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
							const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
							const unsigned int& paramindex,
							vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)_vecArgs;
	if ((vecFilteredM.size()==1) &&(date==vecFilteredM[0].date)){//Nothing to do
		return true;
	} else if ((vecFilteredM.size() < 2) || (vecFilteredM.size() > 3)){
		throw IOException("Not enough data to do resampling or index out of bounds", AT);
	} else if (vecFilteredM.size()==2){
		//add another element
		MeteoData newmd;
		newmd.date = date;
		newmd.setResampled(true);
		vecFilteredM.insert(vecFilteredM.begin() + 1, newmd);
		vecFilteredS.insert(vecFilteredS.begin() + 1, vecFilteredS[0]);
	}

	const MeteoData& tmpmd1 = vecFilteredM[0];
	MeteoData& resampledmd  = vecFilteredM[1];
	const MeteoData& tmpmd2 = vecFilteredM[2];

	double weight = (date.getJulian() - tmpmd1.date.getJulian()) / (tmpmd2.date.getJulian() - tmpmd1.date.getJulian());

	const double& val1 = tmpmd1.param(paramindex);
	const double& val2 = tmpmd2.param(paramindex);
	
	if ((val1 == IOUtils::nodata) || (val2 == IOUtils::nodata)){
		resampledmd.param(paramindex) = IOUtils::nodata;
	} else {
		resampledmd.param(paramindex) = Interpol1D::linearInterpolation(val1, val2, weight);
	}

	return true;
}

bool FilterAlgorithms::MinMaxFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
						 const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
						 const unsigned int& paramindex,
						 vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	vector<double> vecArgs; 
	parseFilterArguments("min_max", _vecArgs, 2, 2, isSoft, vecArgs);

	sort(vecArgs.begin(), vecArgs.end());

	//Run actual MinMax filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		double& value = vecFilteredM[ii].param(paramindex);

		if (value == IOUtils::nodata) continue;

		if (value<vecArgs[0]){
			if (isSoft) value=vecArgs[0];
			else value=IOUtils::nodata;				
			//cout << "Changed: " << value << endl;
		}
		if (value>vecArgs[1]){
			if (isSoft) value=vecArgs[1];
			else value=IOUtils::nodata;			
			//cout << "Changed: " << value << endl;	
		}
	}

	return true;
}

bool FilterAlgorithms::MinValueFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
						   const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
						   const unsigned int& paramindex,
						   vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	vector<double> vecArgs; 
	parseFilterArguments("min", _vecArgs, 1, 1, isSoft, vecArgs);

	//Run actual MinValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		double& value = vecFilteredM[ii].param(paramindex);
		if (value == IOUtils::nodata) continue;
		if (value<vecArgs[0]){
			if (isSoft) value=vecArgs[0];
			else value=IOUtils::nodata;				
		}
	}

	return true;
}

bool FilterAlgorithms::MaxValueFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
						   const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
						   const unsigned int& paramindex,
						   vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	vector<double> vecArgs; 
	parseFilterArguments("max", _vecArgs, 1, 1, isSoft, vecArgs);

	//Run actual MaxValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		double& value = vecFilteredM[ii].param(paramindex);
		if (value == IOUtils::nodata) continue;
		if (value>vecArgs[0]){
			if (isSoft) value=vecArgs[0];
			else value=IOUtils::nodata;				
		}
	}

	return true;
}

void FilterAlgorithms::parseFilterArguments(const string& filtername, const vector<string>& vecArgs_in,
							    const unsigned int& minArgs, const unsigned int& maxArgs, 
							    bool& isSoft, vector<double>& vecArgs_out)
{
	isSoft = false;
	unsigned int argindex = 0;
	vecArgs_out.clear(); 
	double tmp;

	//Test for softness
	if (vecArgs_in.size() > 0)
		if (vecArgs_in[0] == "soft"){ isSoft = true; argindex=1;}

	try {
		for (unsigned int ii=argindex; ii<vecArgs_in.size(); ii++){
			IOUtils::convertString(tmp, vecArgs_in[ii]);
			vecArgs_out.push_back(tmp);
		}

		if ((vecArgs_out.size() < minArgs) || (vecArgs_out.size() > maxArgs))
			throw InvalidArgumentException("Wrong number of arguments for filter " + filtername, AT);
	} catch(std::exception& e){
		cerr << "[E] While processing arguments for filter " << filtername << endl;
		throw;
	}
}
