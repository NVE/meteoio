#include "FilterAlgorithms.h"

 /**
 * @page filters Filters overview
 * The filtering infrastructure is described in FilterAlgorithms (for its API). The goal of this page is to give an overview of the available filters and their usage.
 * 
 * @section filters_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
 *
 * @section filters_available Available filters
 * The filters that are currently available are the following:
 * - rate: calculates the derivative at the current point and reject the point if the derivative is above a given threshold.
 * - min_max: filter all values outside of the interval [min, max], soft operation supported
 * - min: filter out all values less than the given parameter, soft operation supported
 * - max: filter out all values greater than the given parameter, soft operation supported
 *
 * A few data transformations are also supported besides filtering:
 * - resample: resamples (if necessary) the data so that if a required time stamp is between two data points, an interpolated value will be returned (instead of nodata).
 * - accumulation: accumulates the data on a given period. A practical use is to return hourly precipitations from a sensor measuring precipitation on a 10 minutes interval.
 */

std::map<std::string, FilterProperties> FilterAlgorithms::filterMap;
const bool FilterAlgorithms::__init = FilterAlgorithms::initStaticData();

bool FilterAlgorithms::initStaticData()
{
	filterMap["rate"]     = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::RateFilter);
	filterMap["resample"] = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::ResamplingFilter);
	filterMap["median_avg"]  = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MedianAvgFilter);
	filterMap["min_max"]  = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MinMaxFilter);
	filterMap["min"]      = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MinValueFilter);
	filterMap["max"]      = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MaxValueFilter);

	return true;
}

const FilterProperties& FilterAlgorithms::filterProperties(const std::string& filtername)
{
	std::map<std::string, FilterProperties>::const_iterator it;
	it = filterMap.find(filtername);

	if (it==filterMap.end())
		throw UnknownValueException("Unknown filter called: " + filtername, AT);

	return it->second;
}

void FilterAlgorithms::parseFilterArguments(const std::string& filtername, const std::vector<std::string>& vecArgs_in,
							    const unsigned int& minArgs, const unsigned int& maxArgs, 
							    bool& isSoft, std::vector<double>& vecArgs_out)
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
		std::cerr << "[E] While processing arguments for filter " << filtername << std::endl;
		throw;
	}
}


/******************************************************************************
 * The following functions are implementations of different filter algorithms *
 ******************************************************************************/

bool FilterAlgorithms::RateFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						    const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						    const unsigned int& paramindex,
						    std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)date; (void)vecS; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
	parseFilterArguments("rate", _vecArgs, 1, 1, isSoft, vecArgs);
	
	Date_IO deltatime(vecArgs[0]/1440.0); //making a julian date out of the argument given in minutes

	int startposition = pos;
	for (int ii=vecFilteredM.size()-1; ii>=0; ii--){
		unsigned int mypos = startposition;

		double sum = 0.0;
		while((vecM[mypos].date + deltatime) > vecM[startposition].date){
			const double& val = vecM[mypos].param(paramindex);
			if (val != IOUtils::nodata)
				sum += vecM[mypos].param(paramindex);

			//cout << vecM[mypos].date << " HNW:" << vecM[mypos].param(paramindex) << "  SUM:" << sum << endl;		

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

bool FilterAlgorithms::ResamplingFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
							const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
							const unsigned int& paramindex,
							std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)_vecArgs;
	if ((vecFilteredM.size()==1) &&(date==vecFilteredM[0].date)){//Nothing to do
		return false; //Interpretation: filter not applied
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

bool FilterAlgorithms::MinMaxFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						 const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						 const unsigned int& paramindex,
						 std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
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

bool FilterAlgorithms::MinValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						   const unsigned int& paramindex,
						   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
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

bool FilterAlgorithms::MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						   const unsigned int& paramindex,
						   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
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

bool FilterAlgorithms::MedianAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
				   const unsigned int& paramindex,
				   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecS; (void)date; (void)vecFilteredS;
	//Remarks:
	//1) nodata values are not considered when calculating the median
	//2) if there is a even number of window elements the arithmetic mean is used to calculate the median
	//3) isSoft is ignored
	//4) Two arguments expected (both have to be fullfilled for the filter to start operating): 
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	bool isSoft = false;
	std::vector<double> vecArgs; 
	parseFilterArguments("median_avg", _vecArgs, 2, 2, isSoft, vecArgs);

	if (vecArgs[0] < 1) //the window size has to be at least 1
		throw InvalidArgumentException("Parameter 1 for median_avg filter cannot be < 1", AT);
	if (vecArgs[1] < 0) //the time window has to be at least 0 minutes
		throw InvalidArgumentException("Parameter 2 for median_avg filter cannot be < 0", AT);

	//Deal with first parameter: minimal number of data points
	unsigned int windowSize = (unsigned int)vecArgs[0];
	unsigned int increment = (unsigned int)(windowSize/2);
	if (increment > 0) increment--;
	unsigned int startposition = pos + increment;

	if (startposition > (vecM.size()-1))
		startposition = (vecM.size()-1);

	unsigned int endposition = 0;
	if (startposition < (windowSize-1)){
		return false; //not enough data points available
	} else {
		endposition = startposition - (windowSize - 1);
	}

	//Now deal with the second argument: the time window
	Date_IO deltatime(vecArgs[1]/1440.0); //making a julian date out of the argument given in minutes	
	while(deltatime > (vecM[startposition].date - vecM[endposition].date)){
		if (endposition > 0) endposition--;
		else return false; //time window not big enough
	} 
	
	Date_IO gap(vecM[startposition].date - vecM[endposition].date);
	//cout << "The final gap is " << gap.toString() << endl;

	//Run actual median_avg filter over all relevant meteo data
	vector<double> vecWindow;
	//cout << "Start: " << startposition << "  End: " << endposition << endl;
	for (unsigned int ii=startposition; ii>=endposition; ii--){
		const double& tmp = vecM[ii].param(paramindex);
		if (tmp != IOUtils::nodata) vecWindow.push_back(tmp);

		//cout << ii << ": pushed at vecM[" <<  ii << "] " << vecM[ii].date << " : " << tmp << endl;
	}

	//Calculate median
	double median=IOUtils::nodata;
	unsigned int vecSize = vecWindow.size();
	unsigned int middle = (unsigned int)(vecSize/2);
	sort(vecWindow.begin(), vecWindow.end());

	if (vecSize == 0){
		return false; //only nodata values detected or other problem
	} else if ((vecSize % 2) == 1){ //uneven
		median = vecWindow.at(middle);
	} else { //use arithmetic mean of element n/2 and n/2-1
		median = Interpol1D::linearInterpolation(vecWindow.at(middle-1), vecWindow.at(middle), 0.5);
	}

	//cout << "Median value: " << median << endl;
	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		vecFilteredM[ii].param(paramindex) = median;
	}

	return true;
}

