#ifdef _PAROC_
#include "A3DIO.ph"
#else
#include "A3DIO.h"
#endif

using namespace std;

#ifndef _PAROC_
const string A3DIO::ascii_src = "FILE";
const string A3DIO::boschung_src = "BOSCHUNG";
#endif

#ifdef _PAROC_
A3DIO::A3DIO(const std::string& configfile) :  cfg(configfile), fileio(cfg){
  A3DIO::ascii_src = "FILE";
  A3DIO::boschung_src = "BOSCHUNG";
  //load all dynamic plugins
  loadDynamicPlugins();
#else
A3DIO::A3DIO(const std::string& configfile) : IOHandler(NULL), cfg(configfile), fileio(cfg){
  //load all dynamic plugins
  loadDynamicPlugins();
#endif

}

//Copy constructor
#ifdef _PAROC_
//A3DIO::A3DIO(const A3DIO& aio) : cfg(aio.cfg), fileio(cfg), boschungio(cfg){
  //Nothing else so far
//}
#else
A3DIO::A3DIO(const A3DIO& aio) : IOHandler(NULL), cfg(aio.cfg), fileio(cfg) {
  //Nothing else so far
  //TODO: Deal with the IOHandler* pointers, e.g. boschungio
}
#endif

#ifdef _PAROC_
/*A3DIO::A3DIO(const ConfigReader& cfgreader) : cfg(cfgreader), fileio(cfg), boschungio(cfg){
  A3DIO::ascii_src = "FILE";
  A3DIO::boschung_src = "BOSCHUNG";
  loadDynamicPlugins();
}*/
#else
A3DIO::A3DIO(const ConfigReader& cfgreader) : IOHandler(NULL), cfg(cfgreader), fileio(cfg) {
    //Nothing else so far
  loadDynamicPlugins();
}
#endif

#ifdef _PAROC_
A3DIO::~A3DIO(){
#else
A3DIO::~A3DIO() throw(){
#endif
  // Get rid of the object
   
  if (dynLibraryBoschung != NULL) {
    if (boschungio != NULL) {
      boschungio->deleteSelf();
      boschungio = NULL;
    }
    // Close the dynamic library
    delete dynLibraryBoschung;
  }
   
  cleanup();
}

//Clone function
//A3DIO* A3DIO::clone() const { return new A3DIO(*this); }
#ifdef _PAROC_
void A3DIO::cleanup(){
#else
void A3DIO::cleanup() throw(){
#endif
}

void A3DIO::loadDynamicPlugins(){
  string pluginpath = "";

  cout << "[i] " << AT << ": Loading dynamic plugins:" << endl;
  try {
    cfg.getValue("PLUGINPATH", pluginpath); 
    
    //BoschungIO dynamic library needs to be loaded
    cout << "\t" << "Trying to load libBoschungIO.so ... ";
    string filename = pluginpath + "/libBoschungIO.so";
    dynLibraryBoschung = DynamicLoader::loadObjectFile(filename, RTLD_NOW);
    
    if(dynLibraryBoschung == NULL) {
      cout << "failed\n\tCouldn't load the dynamic library " << filename << "\n\t" << dlerror() << endl;
    } else {
      boschungio = dynamic_cast<IOHandler*>(dynLibraryBoschung->newObject("BoschungIO", cfg.getFileName()));
      
      if(boschungio == NULL) {
	cout << "failed" << endl;
	//delete dynLibraryBoschung; This causes a segfault !!
      } else {
	cout << "success" << endl;
      }
    }
  } catch (exception& e){
    if (dynLibraryBoschung != NULL)
      delete dynLibraryBoschung;
    cerr << "\t" << e.what() << endl;
  }
}
 
void A3DIO::get2DGridSize(int& nx, int& ny){
  string demsource="";

  try {
    cfg.getValue("DEMSRC", demsource);

    if (demsource==ascii_src){
      //ASCIIFileIO fileio(cfg.getFileName());
      fileio.get2DGridSize(nx, ny);

    } else if (demsource==boschung_src){
      //Nothing so far
      THROW SLFException("Nothing implemented here", AT);

    } else {
      THROW SLFException("DEMSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
    }

  } catch (...){
    throw;
  }
}

void A3DIO::read2DGrid(Grid2DObject&, const string& filename){
  //Nothing so far
  (void)filename;
  THROW SLFException("Nothing implemented here", AT);
}

void A3DIO::readDEM(Grid2DObject& dem_out){
  string demsource="";

  try {
    cfg.getValue("DEMSRC", demsource); // cout << tmp << endl;

    if (demsource==ascii_src){
      //ASCIIFileIO fileio(cfg.getFileName());
      fileio.readDEM(dem_out);

    } else if (demsource==boschung_src){
      //Nothing so far
      THROW SLFException("Nothing implemented here", AT);

    } else {
      THROW SLFException("DEMSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
    }

  } catch (...){
    throw;
  }
}

void A3DIO::readLanduse(Grid2DObject& landuse_out){
  string landusesource="";

  try {
    cfg.getValue("LANDUSESRC", landusesource); 

    if (landusesource==ascii_src){
      //ASCIIFileIO fileio(cfg.getFileName());
      fileio.readLanduse(landuse_out);

    } else if (landusesource==boschung_src){
      //Nothing so far
      THROW SLFException("Nothing implemented here", AT);

    } else {
      THROW SLFException("LANDUSESRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
    }

  } catch (...){
    throw;
  }
}


void A3DIO::readMeteoData(const Date& date_in, vector<MeteoData>& meteo_out){
  vector<StationData> vecStation;
  readMeteoData(date_in, meteo_out, vecStation);
}

void A3DIO::readMeteoData(const Date& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation){
  //See whether data is already buffered
  //Filter
  //Resample if necessary
  //Filter resampled value
  //return that value
  string meteo2dsource="";

  try {
    cfg.getValue("METEOSRC", meteo2dsource); 

    if (meteo2dsource==ascii_src){
      //ASCIIFileIO fileio(cfg.getFileName());
      fileio.readMeteoData(date_in, vecMeteo, vecStation);

    } else if (meteo2dsource==boschung_src){
      if (boschungio != NULL) {
	boschungio->readMeteoData(date_in, vecMeteo, vecStation);
      } else {
	THROW SLFException("Requesting to read data with plugin libBoschungIO.so, but plugin is not loaded", AT);
      }
    } else {
      THROW SLFException("METEOSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
    }

  } catch (...){
    throw;
  }

}

void A3DIO::readAssimilationData(const Date& date_in, Grid2DObject& da_out){
  string dasrc="";

  try {
    cfg.getValue("DASRC", dasrc); 

    if (dasrc==ascii_src){
      //ASCIIFileIO fileio(cfg.getFileName());
      fileio.readAssimilationData(date_in, da_out);

    } else if (dasrc==boschung_src){
      //Nothing so far
      THROW SLFException("Nothing implemented here", AT);

    } else {
      THROW SLFException("DASRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
    }

  } catch (...){
    throw;
  }
}

void A3DIO::readSpecialPoints(CSpecialPTSArray& pts){
  string specialptssrc="";

  try {
    cfg.getValue("SPECIALPTSSRC", specialptssrc); 

    if (specialptssrc==ascii_src){
      //ASCIIFileIO fileio(cfg.getFileName());
      fileio.readSpecialPoints(pts);

    } else if (specialptssrc==boschung_src){
      //Nothing so far
      THROW SLFException("Nothing implemented here", AT);

    } else {
      THROW SLFException("SPECIALPTSSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
    }

  } catch (...){
    throw;
  }
}
void A3DIO::write2DGrid(const Grid2DObject& grid_in, const string& filename){
  string output="";

  try {
    cfg.getValue("OUTPUT", output); 

    if (output==ascii_src){
      //ASCIIFileIO fileio(cfg.getFileName());
      fileio.write2DGrid(grid_in, filename);
    } else if (output==boschung_src){
      //Nothing so far
      THROW SLFException("Nothing implemented here", AT);

    } else {
      THROW SLFException("OUTPUT does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
    }

  } catch (...){
    throw;
  }
}


