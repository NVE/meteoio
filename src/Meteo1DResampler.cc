#include "Meteo1DResampler.h"

Meteo1DResampler::Meteo1DResampler(){

}

void Meteo1DResampler::resample(const unsigned int& index_in, const Date_IO& date_in, MeteoBuffer& mbuffer_out){
  if ((index_in == 0) && (index_in >= mbuffer_out.size()))
    THROW IOException("[e] Not enough data to do resampling or index out of bounds", AT);

  MeteoData newmd;
  double ta, iswr, vw, rh, lwr, nswc, ts0; 

  MeteoData& tmpmd1 = mbuffer_out.getMeteoData(index_in);
  MeteoData& tmpmd2 = mbuffer_out.getMeteoData(index_in-1);

  double weight = (date_in.getJulian() - tmpmd2.date.getJulian()) / (tmpmd1.date.getJulian() - tmpmd2.date.getJulian());
  //cout << "WEIGHT" << weight << endl;


  if (( tmpmd1.date > date_in) && (tmpmd2.date < date_in)){
    //cout << "[i]-- Resampling date " << date_in.toString()<< endl;
    //cout << "[i]-- Resample date between " << tmpmd1.date.toString() << " and " << tmpmd2.date.toString() << endl;

    
    newmd.date = date_in;

    /*unsigned int rh_left=index_in-1, rh_right=index_in,
      ta_left=index_in-1, ta_right=index_in,
      vw_left=index_in-1, vw_right=index_in,
      nswc_left=index_in-1, nswc_right=index_in;
    */
    //checkForRH(mbuffer_out, rh_left, rh_right);

    /*if ((rh_left==MeteoBuffer::npos) || (rh_right==MeteoBuffer::npos)){
      cerr << "[i] No resampling of rh possible, due to lack of data" << endl;
      
    } else {
      //newmd.rh = 
      }*/

    if ((tmpmd1.ta == MeteoData::nodata) || (tmpmd2.ta == MeteoData::nodata)){
      ta = MeteoData::nodata;
    } else {
      ta = Interpol1D::linearInterpolation(tmpmd2.ta, tmpmd1.ta, weight);
    }

    if ((tmpmd1.rh == MeteoData::nodata) || (tmpmd2.rh == MeteoData::nodata)){
      rh = MeteoData::nodata;
    } else {
      rh = Interpol1D::linearInterpolation(tmpmd2.rh, tmpmd1.rh, weight);
    }
    
    if ((tmpmd1.iswr == MeteoData::nodata) || (tmpmd2.iswr == MeteoData::nodata)){
      iswr = MeteoData::nodata;
    } else {
      iswr = Interpol1D::linearInterpolation(tmpmd2.iswr, tmpmd1.iswr, weight);
    } 

    if ((tmpmd1.vw == MeteoData::nodata) || (tmpmd2.vw == MeteoData::nodata)){
      vw = MeteoData::nodata;
    } else {
      vw = Interpol1D::linearInterpolation(tmpmd2.vw, tmpmd1.vw, weight);
    }

    if ((tmpmd1.lwr == MeteoData::nodata) || (tmpmd2.lwr == MeteoData::nodata)){
      lwr = MeteoData::nodata;
    } else {
      lwr = Interpol1D::linearInterpolation(tmpmd2.lwr, tmpmd1.lwr, weight);
    } 
    
    if ((tmpmd1.nswc == MeteoData::nodata) || (tmpmd2.nswc == MeteoData::nodata)){
      nswc = MeteoData::nodata;
    } else {
      nswc = Interpol1D::linearInterpolation(tmpmd2.nswc, tmpmd1.nswc, weight);
    } 

    if ((tmpmd1.ts0 == MeteoData::nodata) || (tmpmd2.ts0 == MeteoData::nodata)){
      ts0 = MeteoData::nodata;
    } else {
      ts0 = Interpol1D::linearInterpolation(tmpmd2.ts0, tmpmd1.ts0, weight);
    }

    newmd.setMeteoData(date_in, ta, iswr, vw, rh, lwr, nswc, ts0); 

    mbuffer_out.insert(index_in, newmd, mbuffer_out.getStationData(index_in));
  } else {
    THROW IOException("[e] index of meteobuffer not sensible", AT);
  }


}
//void Meteo1DResampler::seekIndices(MeteoBuffer& mbuffer, const string& parameter, unsigned int& leftindex, unsigned int& rightindex){
void Meteo1DResampler::seekIndices(MeteoBuffer&, const string&, unsigned int&, unsigned int&){
  /*  if (parameter == "rh"){
    while ((rightindex < mbuffer.size()) && (mbuffer.getMeteoData(rightindex).rh == MeteoData::nodata)){
      rightindex++;
    }
    while ((leftindex >= 0) && (mbuffer.getMeteoData(rightindex).rh == MeteoData::nodata)){
      leftindex--;
    }

    if ((mbuffer.getMeteoData(rightindex).rh == MeteoData::nodata))
      rightindex == MeteoBuffer::npos;
    if ((mbuffer.getMeteoData(leftindex).rh == MeteoData::nodata))
      leftindex == MeteoBuffer::npos;
  } else if (parameter =="ta"){
    while ((rightindex < mbuffer.size()) && (mbuffer.getMeteoData(rightindex).ta == MeteoData::nodata)){
      rightindex++;
    }
    while ((leftindex >= 0) && (mbuffer.getMeteoData(rightindex).ta == MeteoData::nodata)){
      leftindex--;
    }

    if ((mbuffer.getMeteoData(rightindex).ta == MeteoData::nodata))
      rightindex == MeteoBuffer::npos;
    if ((mbuffer.getMeteoData(leftindex).ta == MeteoData::nodata))
      leftindex == MeteoBuffer::npos;
  } else if (parameter =="iswr"){
    while ((rightindex < mbuffer.size()) && (mbuffer.getMeteoData(rightindex).iswr == MeteoData::nodata)){
      rightindex++;
    }
    while ((leftindex >= 0) && (mbuffer.getMeteoData(rightindex).iswr == MeteoData::nodata)){
      leftindex--;
    }

    if ((mbuffer.getMeteoData(rightindex).iswr == MeteoData::nodata))
      rightindex == MeteoBuffer::npos;
    if ((mbuffer.getMeteoData(leftindex).iswr == MeteoData::nodata))
      leftindex == MeteoBuffer::npos;
  } else if (parameter =="vw"){
    while ((rightindex < mbuffer.size()) && (mbuffer.getMeteoData(rightindex).vw == MeteoData::nodata)){
      rightindex++;
    }
    while ((leftindex >= 0) && (mbuffer.getMeteoData(rightindex).vw == MeteoData::nodata)){
      leftindex--;
    }

    if ((mbuffer.getMeteoData(rightindex).vw == MeteoData::nodata))
      rightindex == MeteoBuffer::npos;
    if ((mbuffer.getMeteoData(leftindex).vw == MeteoData::nodata))
      leftindex == MeteoBuffer::npos;
  } else if (parameter =="lwr"){
    while ((rightindex < mbuffer.size()) && (mbuffer.getMeteoData(rightindex).lwr == MeteoData::nodata)){
      rightindex++;
    }
    while ((leftindex >= 0) && (mbuffer.getMeteoData(rightindex).lwr == MeteoData::nodata)){
      leftindex--;
    }

    if ((mbuffer.getMeteoData(rightindex).lwr == MeteoData::nodata))
      rightindex == MeteoBuffer::npos;
    if ((mbuffer.getMeteoData(leftindex).lwr == MeteoData::nodata))
      leftindex == MeteoBuffer::npos;
  } else if (parameter =="ts0"){
    while ((rightindex < mbuffer.size()) && (mbuffer.getMeteoData(rightindex).ts0 == MeteoData::nodata)){
      rightindex++;
    }
    while ((leftindex >= 0) && (mbuffer.getMeteoData(rightindex).ts0 == MeteoData::nodata)){
      leftindex--;
    }

    if ((mbuffer.getMeteoData(rightindex).ts0 == MeteoData::nodata))
      rightindex == MeteoBuffer::npos;
    if ((mbuffer.getMeteoData(leftindex).ts0 == MeteoData::nodata))
      leftindex == MeteoBuffer::npos;
  } else if (parameter =="nswc"){
        while ((rightindex < mbuffer.size()) && (mbuffer.getMeteoData(rightindex).nswc == MeteoData::nodata)){
      rightindex++;
    }
    while ((leftindex >= 0) && (mbuffer.getMeteoData(rightindex).nswc == MeteoData::nodata)){
      leftindex--;
    }

    if ((mbuffer.getMeteoData(rightindex).nswc == MeteoData::nodata))
      rightindex == MeteoBuffer::npos;
    if ((mbuffer.getMeteoData(leftindex).nswc == MeteoData::nodata))
      leftindex == MeteoBuffer::npos;
  }  
*/
}

