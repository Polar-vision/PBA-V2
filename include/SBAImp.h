

#pragma once
#include "BAExporter.h"
#include "IBA.h"
#include <stdio.h>
#include <map>
#include <vector>
#include <set>
using namespace std;
//using namespace Eigen;

class SBA : public IBA
{
public:
	SBA(void);
	~SBA(void);
	virtual bool ba_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 );

	virtual	bool ba_run( int argc, char** argv );

	virtual bool ba_initialize( char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL );

	virtual bool ba_motstr_levmar( );

	virtual bool ba_motstr_gn(  FixType ft = BA_FixDefault  );

private:
	//compute reprojection error
	void	sba_cost(double *p, double *hx, int* archor );

	//compute Jacobian, not save Jp, Jc etc
	void	sba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );
	void	sba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );

	void	sba_readAndInitialize( char *camsfname, char *ptsfname, int *ncams, int *n3Dpts, int *n2Dprojs,double **motstruct, 
		double **imgpts, int **archor, char **vmask, char **umask,int **nphoto, int** nfeature, int** archorSort );

	void	sba_readProjectionAndInitilizeFeature(	FILE *fp, double *projs, int ncams, 
		int *archor,char* umask,int* nphoto, int* nfeature);
	void sba_saveXYZ(const char* camera, const char* sz3Dpt, double* p);
};

