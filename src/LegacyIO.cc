#ifdef _POPC_
#include "LegacyIO.ph"
#else
#include "LegacyIO.h"
#endif


LegacyIO::LegacyIO(const std::string &meteopath)
{
	meteopathname = meteopath;
	dimx=dimy=dimz=0;
}

LegacyIO::~LegacyIO()
{
}

void LegacyIO::GetGridSize(const std::string& grid_name, int &nx, int &ny, int &nz)
{
	if (dimx>0 && dimy>0 && dimz>0) {
		nx=dimx;
		ny=dimy;
		nz=dimz;
	} else {
		std::string file_name=meteopathname + grid_name + ".asc";
		FILE *fp;

		/* Open the Input file */
		if((fp=fopen(file_name.c_str(),"r")) == 0) {
			throw FileAccessException("Can not open file "+file_name, AT);
		}

		/* Read the dimensions */
		char dummy[MAX_STRING_LENGTH];
		for (int j=0; j<12; j++) {
			fgets(dummy,MAX_STRING_LENGTH,fp);
		}
		if (fscanf(fp,"%d %d %d\n",&nx,&ny,&nz)!=3) {
			fclose(fp);
			throw InvalidFormatException("Can not read nx, ny, nz", AT);
		}

		dimx=nx;
		dimy=ny;
		dimz=nz;
		fclose(fp);
	}
	
}

void LegacyIO::moveToMarker(FILE *fp, const std::string& file_name, const std::string& marker)
{
	char dummy[MAX_LINE_LENGTH];
	do {
		fscanf(fp," %s ",dummy);
	} while (!feof(fp) && strcmp(dummy,marker.c_str()) != 0);
	if(feof(fp)) {
		const std::string message = "End of file "+file_name+" should NOT have been reached when looking for "+marker;
		throw InvalidFormatException(message, AT);
	}
}

void LegacyIO::GetGridPoints(const std::string& grid_name, CDoubleArray &x, CDoubleArray &y, CDoubleArray &z )
{
	int nx, ny, nz;
	GetGridSize(grid_name, nx,ny,nz);

	x.resize(nx);
	y.resize(ny);
	z.resize(nz*ny*nx);

	std::string file_name = meteopathname + grid_name + ".asc";
	FILE *fp;

	/* Open the Input file */
	if((fp=fopen(file_name.c_str(),"r")) == 0) {
		throw FileAccessException("Can not open file "+file_name, AT);
	}

	/* Read the x- and y-coordinates */
	moveToMarker(fp, file_name, "x_coordinate");
	for (int j=0; j<nx; j++) {
		fscanf(fp," %16lf%*[\n]",&x[j]);
	}
	
	moveToMarker(fp, file_name, "y_coordinate");
	for (int j=0; j<ny; j++) {
		fscanf(fp," %16lf%*[\n]",&y[j]);
	}
  
	/* Read until zp-coordinate found */
	moveToMarker(fp, file_name, "zp_coordinat");
	for (int j = 0; j < nx*ny*nz; j++) {
		fscanf(fp," %16lf%*[\n]",&z[j]);
	}
	
	fclose(fp);
}

void LegacyIO::GetGridData(CNodeArray &nodes, const std::string& hour)
{
	std::string file_name = meteopathname + hour + ".dat";
	FILE *fp;

	int nx, ny, nz;
	GetGridSize(file_name, nx,ny,nz);

	nodes.resize(nx*ny*nz);

	if ((fp=fopen(file_name.c_str(),"r")) == 0) {
		throw FileAccessException("Can not open file "+file_name, AT);
	}

	DEBUG("Reading Wind Data from File %s",file_name.c_str());
	
	/* Read u */
	moveToMarker(fp, file_name, "u");
	for (int iz=0; iz<nz; iz++) for (int iy=ny-1; iy>=0; iy--) for (int ix=0; ix<nx; ix++) {
		const int j = iz*nx*ny + iy*nx + ix;
		fscanf(fp," %16lf%*[\n]",&nodes[j].u);
	}
	
	/* Read v */
	moveToMarker(fp, file_name, "v");
	for (int iz=0; iz<nz; iz++) for (int iy=ny-1; iy>=0; iy--) for (int ix=0; ix<nx; ix++) {
		const int j = iz*nx*ny + iy*nx + ix;
		fscanf(fp," %16lf%*[\n]",&nodes[j].v);
	}
	
	/* Read w */
	moveToMarker(fp, file_name, "w");
	for (int iz=0; iz<nz; iz++) for (int iy=ny-1; iy>=0; iy--) for (int ix=0; ix<nx; ix++) {
		const int j = iz*nx*ny + iy*nx + ix;
		fscanf(fp," %16lf%*[\n]",&nodes[j].w);
	}

	fclose(fp);

#ifdef _POPC_
	DEBUG("Read wind field for hour %s OK",hour.c_str());
#endif
}

