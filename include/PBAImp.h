#pragma once
#include "BAExporter.h"
#include "IBA.h"
#include <map>
#include <vector>
#include <set>
using namespace std;
//using namespace Eigen;

class PBA : public IBA
{
public:
	PBA(void);
	~PBA(void);
	virtual bool ba_run(bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL,
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6);

	virtual	bool ba_run( int argc, char** argv );

	virtual bool ba_initialize( char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL );

	virtual bool ba_motstr_levmar( );

	virtual bool ba_motstr_gn(  FixType ft = BA_FixDefault  );

private:
	//initialize feature. You can provide xyz or system also provide them by itself;
	bool    pba_initializeMainArchor( double* imgpts, double* camera,double* K,double* feature, int nP, int FID, double* KR );
	bool    pba_initializeAssoArchor( double* imgpts, int* photo, double* camera,double* K,double* feature,int nMI, int nAI, int FID, bool bLast );
	bool	pba_initializeOtheArchors( double* imgpts, int* photo, double* camera,double* K,double* feature,int* archorSort,int nfeacout, int nOI, int FID );

	//compute reprojection error
	void	pba_cost(double *p, double *hx, int* archor );

	//compute Jacobian, not save Jp, Jc etc
	void	pba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );
	void	pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );

	void	pba_readAndInitialize( char *camsfname, char *ptsfname, int *ncams, int *n3Dpts, int *n2Dprojs,double **motstruct, 
		double **imgpts, int **archor, char **vmask, char **umask,int **nphoto, int** nfeature, int** archorSort );

	void	pba_readProjectionAndInitilizeFeature(	FILE *fp, double *params, double *projs, char *vmask, int ncams, 
		int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort );


	//Add By Zuo
	bool pba_initializeOtheArchors_Mindw(double* imgpts, int* photo, double* camera, double* K, double* feature, int* archorSort, int nfeacout, int nOI, int FID);
	void pba_saveInitialXYZ(const char* sz3Dpt, double* p);
	void pba_saveInitialParallax(const char* sz3Dpt, double* p);
	double* pba_angle2xyz(double* p);

	//transform angle into XYZ
	int		pba_angle2xytGN(double* p);
	int		pba_angle2xytLM(double* p);
	void	pba_saveXYZ(const char* camera, const char* sz3Dpt, double* p, bool gn = true);
	int zu = 0;
};

