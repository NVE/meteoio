#ifdef _POPC_
#include "LegacyIO.ph"
#else
#include "LegacyIO.h"
#endif


LegacyIO::LegacyIO(char *meteopath)
{
	strcpy(meteopathname,meteopath);
	dimx=dimy=dimz=0;

#ifdef _POPC_
	*cache_Hour=0;
#endif
}

LegacyIO::~LegacyIO()
{
#ifdef _POPC_
	printf("Total time for reading grid data: %g seconds\n", timer.Elapsed());
#endif
}

void LegacyIO::GetGridSize(int &nx, int &ny, int &nz)
{
	if (dimx>0 && dimy>0 && dimz>0) {
		nx=dimx;
		ny=dimy;
		nz=dimz; 
	} else {
		char dummy[MAX_STRING_LENGTH], fname[MAX_STRING_LENGTH];
		FILE *fp;

		strcpy(fname,meteopathname);
		strcat(fname,".dat");
		/* Open the Input file */
		if((fp=fopen(fname,"r")) == 0) {
			DEBUG("cannot open file %s",fname);
			throw ERROR_INPUT_GRIDSIZE;
		}

		/* Read the dimensions */
		do {
			fscanf(fp," %s ",dummy);;
		} while (strcmp(dummy,"nx") != 0);
		
		if (fscanf(fp," = %d, ny = %d, nz = %d",&nx,&ny,&nz)!=3) {
			fclose(fp);
			throw ERROR_INPUT_GRIDSIZE; 
		}

		dimx=nx;
		dimy=ny;
		dimz=nz;
		fclose(fp);
	}
}

void LegacyIO::GetGridPoints(CDoubleArray &x, CDoubleArray &y, CDoubleArray &z )
{
	int nx, ny, nz;
	GetGridSize(nx,ny,nz);

	x.SetSize(nx);
	y.SetSize(ny);
	z.SetSize(nz*ny*nx);

	char dummy[MAX_LINE_LENGTH], fname[MAX_STRING_LENGTH];
	int j;

	FILE *fp;

	/* Open the Input file */
	strcpy(fname,meteopathname);
	strcat(fname,".dat");
	/* Open the Input file */
	if((fp=fopen(fname,"r")) == 0) {
		DEBUG("cannot open file %s",fname);
		throw ERROR_INPUT_GRIDSIZE;
	}


	/* Read until the coordinates are found */
	/* Read the dimensions */
	do {
		fscanf(fp," %s ",dummy);;
	} while (strcmp(dummy,"coordinate") != 0);

	/* Read the x- and y-coordinates */
	// fgets(dummy,MAX_LINE_LENGTH,fp);
	for ( j=0; j<nx; j++) {
		fscanf(fp," %16lf%*[\n]",&x[j]);
	}
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for ( j=0; j<ny; j++) {
		fscanf(fp," %16lf%*[\n]",&y[j]);
	}
	for ( j=0; j<6; j++) {
		fgets(dummy,MAX_LINE_LENGTH,fp);
	}

	/* Calculate/read the grid and read the topography (z-values) */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (int iz=0; iz<nz; iz++) {
		for (int iy=ny-1; iy>=0; iy--) {
			for (int ix=0; ix<nx; ix++)
			{
				j = iz*nx*ny + iy*nx + ix;
				fscanf(fp," %16lf%*[\n]",&z[j]);
			}
		}
	}	

	fclose(fp);
}


void LegacyIO::GetGridData(CNodeArray &nodes, char *hour)
{
#ifdef _POPC_ 
	if (hour!=NULL && 0==strcmp(hour,cache_Hour) && cache_WindField.size()) {
		DEBUG("GET WIND FIELD DATA FROM CACHE");
		nodes=cache_WindField;
		cache_WindField.RemoveAll();
		*cache_Hour=0;
		return;
	}
	timer.Start();
#endif
	int nx, ny, nz;
	GetGridSize(nx,ny,nz);

	nodes.SetSize(nx*ny*nz);
	FILE *fp;
	char dummy[MAX_LINE_LENGTH], fname[MAX_STRING_LENGTH];

	strcpy(fname,meteopathname);
	strcat(fname,hour);
	strcat(fname,".dat");

	if ((fp=fopen(fname,"r")) == 0) {
		DEBUG("cannot open file %s",fname);
		throw ERROR_INPUT_GRIDDATA; 
	}

	DEBUG("Reading Wind Data from File %s",fname);

	/* Read the parameters */
	/* Skip header and Read ubar */
	do {
		fscanf(fp," %s ",dummy);;
	} while (strcmp(dummy,"u") != 0);

	int iz,iy,ix,j;

	/* Read u */
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].u);
			}
		}
	}	
	
	/* Read v */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].v);
			}
		}
	}
				
	/* Read w */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].w);
			}
		}
	}

	/* Read Potential Temperature */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].tet);
			}
		}
	}

	/* Read Pressure */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].p);
			}
		}
	}	

	/* Read TKE and see whether all had been read correctly */
	do {
		fscanf(fp," %s ",dummy);;
	} while (strcmp(dummy,"tke") != 0);
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].e);
			}
		}
	}
			
	if (strcmp(dummy,"tke") != 0) {
		DEBUG("Error occurred while reading ARPS data %s",dummy);
		throw ERROR_INPUT_GRIDDATA;
	}




	/* Read Km, the horizontal eddy diffusivity */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].Km);
			}
		}
	}

	/* Read lm, the mixing length but in really here in Arps5 it is now the vertical diffusivity */
	do {
		fscanf(fp," %s ",dummy);;
	} while (strcmp(dummy,"kmv") != 0);
	
	for (iz=0; iz<nz; iz++) {
		for (iy=ny-1; iy>=0; iy--) {
			for (ix=0; ix<nx; ix++) {
				j = iz*nx*ny + iy*nx + ix;  fscanf(fp," %16lf%*[\n]",&nodes[j].lm);
			}
		}
	}

	if (strcmp(dummy,"kmv") != 0){
		DEBUG("Error occurred while reading ARPS data %s",dummy);
		throw ERROR_INPUT_GRIDDATA;
	}

	fclose(fp);

#ifdef _POPC_ 
	timer.Stop();
	DEBUG("Read wind field for hour %s OK",hour);
#endif
}

void LegacyIO::PrepareNextWindField(char *hour)
{
#ifdef _POPC_
	cache_WindField.RemoveAll();
	try {
		DEBUG("Read grid data (hour=%s) for caching",hour);
		GetGridData(cache_WindField,hour);
		strcpy(cache_Hour,hour);
	}
	catch (...) {
		DEBUG("Can not read next wind field for caching");
	}
#else
	(void)hour;
#endif
}
