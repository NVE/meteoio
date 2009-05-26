#ifndef __A3DIO_H__
#define __A3DIO_H__

#include "IOHandler.h"
#include "ASCIIFileIO.h"
#include "LegacyIO.ph"
#include "slfexceptions.h"

/*#include "ConfigReader.h"
#include "StationData.h"
#include "MeteoData.h"
#include "Grid2DObject.h"
#include "Date.h"
#include "Alpine3D.h"*/

#include "marshal_Alpine3D.h"

parclass A3DIO;

parclass A3DIO{ // Note : No heritage here for POPC++ : a parclass cannot herit from a class
  classuid(1003);
 public:
  A3DIO(const string& configfile) @{ power=100 ?: 50; };
  ~A3DIO();

  virtual void get2DGridSize(int& nx, int& ny);
  virtual void read2DGrid([out]Grid2DObject& dem_out, const string& parameter="");

  virtual void readDEM([out]Grid2DObject& dem_out);
  virtual void readLanduse([out]Grid2DObject& landuse_out);

  virtual void readMeteoData([in]const Date& date_in,
                             [out, proc=marshal_vector_MeteoData]vector<MeteoData>& vecMeteo);
  
  virtual void readMeteoData([in]const Date& date_in,
                             [out, proc=marshal_vector_MeteoData]vector<MeteoData>& vecMeteo,
                             [out, proc=marshal_vector_StationData]vector<StationData>& vecStation);
  
  virtual void readAssimilationData([in] const Date&,[out] Grid2DObject& da_out);
  virtual void readSpecialPoints([out,proc=marshal_CSpecialPTSArray]CSpecialPTSArray& pts);

  virtual void write2DGrid([in]const Grid2DObject& grid_in, const string& filename);

 private:
  string ascii_src;
  string boschung_src;

 private:
  void cleanup();// throw();
  void loadDynamicPlugins();
  
  ConfigReader cfg;
  ASCIIFileIO fileio;
  DynamicLibrary* dynLibraryBoschung;
  IOHandler* boschungio;
};

#endif
