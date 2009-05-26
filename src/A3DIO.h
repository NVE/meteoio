#ifndef __A3DIO_H__
#define __A3DIO_H__

#ifdef _PAROC_
#error
#endif

#include "IOHandler.h"
#include "ASCIIFileIO.h"
#include "slfexceptions.h"

class A3DIO : public IOHandler {
 public:
  // virtual A3DIO* clone() const; // lwk : not used yet

  A3DIO(const std::string& configfile);
  A3DIO(const A3DIO&);
  A3DIO(const ConfigReader&);
  ~A3DIO() throw();

  virtual void get2DGridSize(int& nx, int& ny);
  virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");

  virtual void readDEM(Grid2DObject& dem_out);
  virtual void readLanduse(Grid2DObject& landuse_out);

  virtual void readMeteoData(const Date& date_in, vector<MeteoData>& vecMeteo);
  virtual void readMeteoData(const Date& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);

  virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
  virtual void readSpecialPoints(CSpecialPTSArray& pts);

  virtual void write2DGrid(const Grid2DObject& grid_in, const string& filename);

  static const string ascii_src;
  static const string boschung_src;

 private:
  void cleanup() throw();
  void loadDynamicPlugins();

  ConfigReader cfg;
  ASCIIFileIO fileio;
  DynamicLibrary* dynLibraryBoschung;
  IOHandler* boschungio;
};

#endif
