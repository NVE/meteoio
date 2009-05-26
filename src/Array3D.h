#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <errno.h>
#include <memory.h>
#include <stdlib.h>

template<class T> class CArray3D
{
 public:
  CArray3D();
  CArray3D(int anx, int any, int anz);
  ~CArray3D();
  void Create(int anx, int any, int anz);

  void GetSize(int &anx, int &any, int &anz);

  void Destroy();

  inline operator T***();
 protected:
  T ***data;
  int nx;
  int ny;
  int nz;
};

template<class T> 
CArray3D<T>::CArray3D()
{
  nx=ny=nz=0;
  data=0;
}

template<class T> 
CArray3D<T>::CArray3D(int anx, int any, int anz)
{
  nx=ny=nz=0;
  data=0;
  Create(anx,any,anz);
}

template<class T> 
CArray3D<T>::~CArray3D()
{
  Destroy();
}

template<class T> 
void CArray3D<T>::Create(int anx, int any, int anz)
{
  if (anx==nx && any==ny && anz==nz) return;
  
  if (anx<=0 || any<=0 || anz<=0) 
    {
      Destroy();
      return;
    }

  T ***tmpdat;
  //1st dimension
  if ( (tmpdat=new T **[anx])==NULL) throw errno;

  //2nd dimension
  for (int i=0;i<anx;i++)
    {
      if ( (tmpdat[i]=new T * [any])==NULL)
	{
	  for (int k=0;k<i;k++) delete tmpdat[k];
	  delete tmpdat;
	  throw errno;
	}
    }

  //3rd dimension....
  for (int i=0;i<anx;i++)
    for (int j=0;j<any;j++)
      if ( (tmpdat[i][j]=new T[anz])==NULL)
	{
	  //delete all allocated memory and throw errno
	  for (int i1=0;i1<i;i1++)
	    for (int j1=0;j1<ny;j1++) if (tmpdat[i1][j1]!=NULL) delete tmpdat[i1][j1];
	  for (int j1=0;j1<j;j1++) if (tmpdat[i][j1]!=NULL) delete tmpdat[i][j1];

	  for (int i1=0;i1<anx;i1++) delete tmpdat[i1];
	  delete tmpdat;
	  throw errno;
	}

  Destroy();
  data=tmpdat;
  nx=anx;
  ny=any;
  nz=anz;
}

template<class T> 
void CArray3D<T>::GetSize(int &anx, int &any, int &anz)
{
  anx=nx;
  any=ny;
  anz=nz;
}

template<class T> 
void CArray3D<T>::Destroy()
{
  if (data==NULL) return;
  for (int i=0;i<nx;i++) 
    if (data[i]!=NULL)
      {
	for (int j=0;j<ny;j++) if (data[i][j]!=NULL) delete data[i][j];
	delete data[i];
      }
  delete data;
  data=NULL;
  nx=ny=nz=0;
}

template<class T> 
inline CArray3D<T>::operator T***()
{
  return data;
}

#endif
