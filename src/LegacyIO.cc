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
		for (int j=0; j<8; j++) fgets(dummy,MAX_STRING_LENGTH,fp);
		if (fscanf(fp,"%d %d %d\n",&nx,&ny,&nz)!=3) {
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

	x.resize(nx);
	y.resize(ny);
	z.resize(nz*ny*nx);

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
		cache_WindField.clear();
		*cache_Hour=0;
		return;
	}
	timer.Start();
#endif
	int nx, ny, nz;
	GetGridSize(nx,ny,nz);

	nodes.resize(nx*ny*nz);
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
		fscanf(fp," %s ",dummy);
	}
	while (strcmp(dummy,"ubar") != 0);
	
	int iz,iy,ix,j;
	double num;
	
	/* Read ubar */
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].u);}
	
	/* Read vbar */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].v);}
	/* Read wbar */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].w);}
	/* Read tetbar */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].tet);}
	/* Read pbar */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].p);}
	/* Read u */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&num); nodes[j].u += num;
	if (nodes[j].u == 0.) 
		printf("\n Zero u at ix %d, iy %d, iz %d", ix, iy, iz);}
	/* Read v */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&num); nodes[j].v = -nodes[j].v - num;}
	/* Read w */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&num); nodes[j].w += num;}
	/* Read tet */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&num); nodes[j].tet += num;}
	/* Read p */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&num); nodes[j].p += num;}
	/* Read Km, the eddy diffusivity */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].Km);}
	/* Read lm, the mixing length */
	fgets(dummy,MAX_LINE_LENGTH,fp);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix;  fscanf(fp," %16lf%*[\n]",&nodes[j].lm);}
	/* Read TKE and see whether all had been read correctly */
	fscanf(fp," %s",dummy);
	for (iz=0; iz<nz; iz++) for (iy=ny-1; iy>=0; iy--) for (ix=0; ix<nx; ix++) {
		j = iz*nx*ny + iy*nx + ix; fscanf(fp," %16lf%*[\n]",&nodes[j].e);}
	if (strcmp(dummy,"tke") != 0){
		DEBUG("Error occurred while reading ARPS data %s",dummy);
		throw ERROR_INPUT_GRIDDATA; }
	
	fclose(fp);

#ifdef _POPC_ 
	timer.Stop();
	DEBUG("Read wind field for hour %s OK",hour);
#endif
}

void LegacyIO::PrepareNextWindField(char *hour)
{
#ifdef _POPC_
	cache_WindField.clear();
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
