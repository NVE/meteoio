#ifndef LEGACYIO_H
#define LEGACYIO_H

#include "timer.h"

#include "Alpine3D.h"
#include "DriftData.h"
#include "Snowpack.h"

#include "marshal_Alpine3D.h"


parclass LegacyIO
{
 public:
  LegacyIO( [in, proc=marshalstring, size=256] char *meteopath);
  ~LegacyIO();

  virtual void GetGridSize([out] int &nx, [out] int &ny, [out] int &nz);
  virtual void GetGridPoints([out, proc=marshal_CDoubleArray] CDoubleArray &x,[out, proc=marshal_CDoubleArray]  CDoubleArray &y,[out, proc=marshal_CDoubleArray]  CDoubleArray &z);
  virtual void GetGridData([out, proc=marshal_input_CNodeArray] CNodeArray &data, [in, proc=marshalstring, size=256] char *hour);

  async virtual void PrepareNextWindField([in, proc=marshalstring, size=256] char *hour);
  
  classuid(1002);

 private:
  char demfilename[MAX_STRING_LENGTH];
  char meteopathname[MAX_STRING_LENGTH];
  int dimx, dimy, dimz;

  //For caching data
  char cache_Hour[MAX_STRING_LENGTH];
  CNodeArray cache_WindField;

  Timer timer;
};

#endif
