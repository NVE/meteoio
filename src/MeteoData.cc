#include "MeteoData.h"

using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/
const unsigned int MeteoData::nrOfParameters =  MeteoData::lastparam - MeteoData::firstparam + 1;
map<unsigned int, string> MeteoData::meteoparamname;
const bool MeteoData::__init = MeteoData::initStaticData();

bool MeteoData::initStaticData()
{
	//This function should only be executed once for all MeteoData instances
	//Associate unsigned int value and a string representation of a meteo parameter
	meteoparamname[TA]   = "TA";
	meteoparamname[ISWR] = "ISWR";
	meteoparamname[VW]   = "VW";
	meteoparamname[DW]   = "DW";
	meteoparamname[RH]   = "RH";
	meteoparamname[LWR]  = "LWR";
	meteoparamname[HNW]  = "HNW";
	meteoparamname[TSG]  = "TSG";
	meteoparamname[TSS]  = "TSS";
	meteoparamname[HS]   = "HS";
	meteoparamname[RSWR] = "RSWR";
	meteoparamname[P]    = "P";

	return true;
}

const string& MeteoData::getParameterName(const unsigned int& parindex)
{
	if (parindex >= MeteoData::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);

	return MeteoData::meteoparamname[parindex];
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

void MeteoData::initParameterMap()
{
	//The following mapping needs to be done for every instance of MeteoData
	meteoparam[TA]       = &ta;
	meteoparam[ISWR]     = &iswr;
	meteoparam[VW]       = &vw;
	meteoparam[DW]       = &dw;
	meteoparam[RH]       = &rh;
	meteoparam[LWR]      = &lwr;
	meteoparam[HNW]      = &hnw;
	meteoparam[TSG]      = &tsg;
	meteoparam[TSS]      = &tss;
	meteoparam[HS]       = &hs;
	meteoparam[RSWR]     = &rswr;
	meteoparam[P]        = &p;

	//Check for inconsistency between enum Parameters and the two maps meteoparam and meteoparamname
	if ((meteoparam.size() != meteoparamname.size()) || (meteoparam.size() != MeteoData::nrOfParameters))
		throw IOException("Inconsistency within class MeteoData: Check function initMaps()", AT);
}

MeteoData::MeteoData() : resampled(false)
{
	setMeteoData(Date_IO(0.0), nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata);
	initParameterMap();
}

MeteoData::MeteoData(const Date_IO& date_in, const double& ta_in, const double& iswr_in, 
				 const double& vw_in, const double& dw_in, const double& rh_in,
				 const double& lwr_in, const double& hnw_in, const double& tsg_in, 
				 const double& tss_in, const double& hs_in, const double& rswr_in, const double& _p) : resampled(false)
{
	setMeteoData(date_in, ta_in, iswr_in, vw_in, dw_in, rh_in, lwr_in, hnw_in, tsg_in, tss_in, hs_in, rswr_in, _p);
	initParameterMap();
}

MeteoData::MeteoData(const MeteoData& md)
{
	*this = md;
}

MeteoData& MeteoData::operator=(const MeteoData& rhs)
{
	if (this == &rhs) //Test self assignment
		return *this;

	date = rhs.date;
	resampled = rhs.resampled;

	initParameterMap();

	std::map<unsigned int, double*>::const_iterator it;
	for (it=rhs.meteoparam.begin(); it!=rhs.meteoparam.end(); it++){
		*meteoparam[it->first] = *(it->second);
	}

	return *this;
}

void MeteoData::setMeteoData(const Date_IO& date_in, const double& ta_in, const double& iswr_in, const double& vw_in,
					    const double& dw_in, const double& rh_in, const double& lwr_in, const double& hnw_in,
					    const double& tsg_in, const double& tss_in, const double& hs_in, const double& rswr_in, 
					    const double& _p)
{
	date = date_in;
	ta = ta_in;
	iswr = iswr_in;
	vw = vw_in;
	dw = dw_in;
	rh = rh_in;
	lwr = lwr_in;
	hnw = hnw_in;
	tsg = tsg_in;
	tss = tss_in;
	hs = hs_in;
	rswr = rswr_in;
	p=_p;
}

bool MeteoData::isResampled()
{
	return resampled;
}

void MeteoData::setResampled(const bool& _resampled)
{
	resampled = _resampled;
}

/*
void MeteoData::Check_min_max(double& param, const double low_hard, const double low_soft, const double high_soft, const double high_hard)
{
	//This is the embryo of what would become the filtering library...
	//HACK: TODO: have the limits being passed as parameters from a config file
	//(so the user could adjust them to his context and units!)
	if(param<=nodata || param<=low_hard || param>=high_hard) {
		param=nodata;
	} else {
		if(param<=low_soft) {
			param=low_soft;
		} else if(param>=high_soft) {
			param=high_soft;
		}
	}
}

void MeteoData::cleanData()
{
	//HACK: TODO: have the limits being passed as parameters from a config file
	//HACK TODO: get rid of this method (as well as Check_min_max) and replace it by calls to the filters

	Check_min_max(ta, -50., -50., 50., 50.);

	Check_min_max(iswr, -50., 0., 2000., 2000.);

	Check_min_max(vw, -60., -60., 60., 60.);

	Check_min_max(dw, -180., -180., 180., 180.);

	Check_min_max(rh, -50., 0., 100., 150.);

	Check_min_max(hnw, -5., 0., 100., 100.);

	Check_min_max(tsg, -70., -70., 70., 70.);

	Check_min_max(tss, -70., -70., 70., 70.);

	Check_min_max(lwr, -50., 0., 2000., 2000.);
	
	Check_min_max(hs, -5., 0., 350., 350.);
	
	Check_min_max(rswr, -50., 0., 2000., 2000.);
}
*/

bool MeteoData::operator==(const MeteoData& in) const
{
	//An object is equal if the date is equal and all meteo parameters are equal
	bool eval = (date==in.date);

	for (unsigned int ii=0; ii<MeteoData::nrOfParameters; ii++){
		eval &= (param(ii) == in.param(ii));
	}

	return eval;
}

bool MeteoData::operator!=(const MeteoData& in) const
{
	return !(*this==in);
}

double& MeteoData::param(const unsigned int& parindex)
{
#ifndef NOSAFECHECKS
	if (parindex >= MeteoData::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif	
	return *(meteoparam[parindex]);
}

const double& MeteoData::param(const unsigned int& parindex) const
{
	std::map<unsigned int, double*>::const_iterator it;
	it = meteoparam.find(parindex);

#ifndef NOSAFECHECKS
	if (it == meteoparam.end())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return *(it->second);
}

const string MeteoData::toString() const
{
	stringstream tmpstr;

	tmpstr << setprecision(10) << "Date_IO: " << date.toString() << endl; 

	std::map<unsigned int, double*>::const_iterator it1;
	for (it1=meteoparam.begin(); it1 != meteoparam.end(); it1++){
		tmpstr << setw(7) << MeteoData::getParameterName(it1->first) << ":" << setw(15) << *it1->second << endl;
	}	

	return tmpstr.str();
}

#ifdef _POPC_
void MeteoData::Serialize(POPBuffer &buf, bool pack)
{//HACK TODO: adjust the serialize to the latest changes!! (call InitStaticData(), etc)
	if (pack){
		buf.Pack(&resampled,1);
		date.Serialize(buf,true);
		buf.Pack(&ta,1);
		buf.Pack(&iswr,1);
		buf.Pack(&vw,1);
		buf.Pack(&dw,1);
		buf.Pack(&rh,1);
		buf.Pack(&lwr,1);
    //buf.Pack(&ea,1);
		buf.Pack(&hnw,1);
		buf.Pack(&tsg,1);
		buf.Pack(&tss,1);
		buf.Pack(&hs,1);
		buf.Pack(&rswr,1);
	}else{
		buf.UnPack(&resampled,1);
		date.Serialize(buf,false);
		buf.UnPack(&ta,1);
		buf.UnPack(&iswr,1);
		buf.UnPack(&vw,1);
		buf.UnPack(&dw,1);
		buf.UnPack(&rh,1);
		buf.UnPack(&lwr,1);
    //buf.UnPack(&ea,1);
		buf.UnPack(&hnw,1);
		buf.UnPack(&tsg,1);
		buf.UnPack(&tss,1);
		buf.UnPack(&hs,1);
		buf.UnPack(&rswr,1);
		initParameterMap();
	}
}
#endif

