#include "MeteoData.h"

using namespace std;

MeteoData::MeteoData()
{
	setMeteoData(Date_IO(0.0), nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata, nodata);
}

MeteoData::MeteoData(const Date_IO& date_in, const double& ta_in, const double& iswr_in, const double& vw_in, const double& dw_in, const double& rh_in,
		const double& lwr_in, const double& nswc_in, const double& tsg_in, const double& tss_in, const double& hs_in, const double& rswr_in)
{
	setMeteoData(date_in, ta_in, iswr_in, vw_in, dw_in, rh_in, lwr_in, nswc_in, tsg_in, tss_in, hs_in, rswr_in);
}

void MeteoData::setMeteoData(const Date_IO& date_in, const double& ta_in, const double& iswr_in, const double& vw_in,
				const double& dw_in, const double& rh_in, const double& lwr_in, const double& nswc_in,
				const double& tsg_in, const double& tss_in, const double& hs_in, const double& rswr_in)
{

	date = date_in;
	ta = ta_in;
	iswr = iswr_in;
	vw = vw_in;
	dw = dw_in;
	rh = rh_in;
	lwr = lwr_in;
	nswc = nswc_in;
	tsg = tsg_in;
	tss = tss_in;
	hs = hs_in;
	rswr = rswr_in;

}

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

	Check_min_max(nswc, -5., 0., 100., 100.);

	Check_min_max(tsg, -70., -70., 70., 70.);

	Check_min_max(tss, -70., -70., 70., 70.);

	Check_min_max(lwr, -50., 0., 2000., 2000.);
	
	Check_min_max(hs, -5., 0., 350., 350.);
	
	Check_min_max(rswr, -50., 0., 2000., 2000.);

}

bool MeteoData::operator==(const MeteoData& in) const
{
	return ((date==in.date) && (ta==in.ta) 
				&& (iswr==in.iswr) && (vw==in.vw) 
				&& (rh==in.rh) && (lwr==in.lwr) 
				&& (dw==in.dw) && (tss==in.tss)
				&& (nswc==in.nswc) && (tsg==in.tsg)
				&& (hs==in.hs) && (rswr==in.rswr));
}

bool MeteoData::operator!=(const MeteoData& in) const
{
	return !(*this==in);
}

const string MeteoData::toString() const
{
	stringstream tmpstr;

	tmpstr << setprecision(10) << "Date_IO: " << date.toString() << endl 
				   << setw(6) << "ta: " << setw(15) << ta << endl
				   << setw(6) << "iswr: " << setw(15) << iswr << endl
				   << setw(6) << "vw: " << setw(15) << vw  << endl
				   << setw(6) << "dw: " << setw(15) << dw  << endl
				   << setw(6) << "rh: " << setw(15) << rh << endl
				   << setw(6) << "lwr: " << setw(15) << lwr << endl
				   << setw(6) << "nswc: " << setw(15) << nswc << endl
				   << setw(6) << "tsg: " << setw(15) << tsg << endl
				   << setw(6) << "tss: " << setw(15) << tss << endl
				   << setw(6) << "hs: " << setw(15) << hs << endl
				   << setw(6) << "rswr: " << setw(15) << rswr << endl;

	return tmpstr.str();
}

#ifdef _POPC_
void MeteoData::Serialize(paroc_buffer &buf, bool pack)
{
	if (pack){
		date.Serialize(buf,true);
		buf.Pack(&ta,1);
		buf.Pack(&iswr,1);
		buf.Pack(&vw,1);
		buf.Pack(&dw,1);
		buf.Pack(&rh,1);
		buf.Pack(&lwr,1);
    //buf.Pack(&ea,1);
		buf.Pack(&nswc,1);
		buf.Pack(&tsg,1);
		buf.Pack(&tss,1);
		buf.Pack(&hs,1);
		buf.Pack(&rswr,1);
	}else{
		date.Serialize(buf,false);
		buf.Pack(&ta,1);
		buf.Pack(&iswr,1);
		buf.Pack(&vw,1);
		buf.Pack(&dw,1);
		buf.Pack(&rh,1);
		buf.Pack(&lwr,1);
    //buf.Pack(&ea,1);
		buf.Pack(&nswc,1);
		buf.Pack(&tsg,1);
		buf.Pack(&tss,1);
		buf.Pack(&hs,1);
		buf.Pack(&rswr,1);
	}

  
  
}
#endif

