
#ifndef _SBA_H
#define _SBA_H

#ifdef SBA_EXPORTS
	#define SBAapi  __declspec(dllexport)
#else
	#define SBAapi __declspec(dllimport)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/StdVector>
#include <eigen3/Eigen/Cholesky>

#include "cholmod/cholmod.h"

#define PI  3.1415926535898 
#define MAXARCHOR 0.5

typedef enum SBA_FixType
{	
	SBA_FixDefault = -1, //fix the longest axis  
	SBA_FixX = 0 ,	    // fix X axis	
	SBA_FixY = 1 ,	    // fix Y axis
	SBA_FixZ = 2 ,	    // fix Z axis
} FixType ; 
	
class SBAapi ISBA
{
public:
	virtual bool sba_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 ) = 0;

	virtual	bool sba_run( int argc, char** argv ) = 0;

	virtual bool sba_initialize( char* szCamera, char* szFeature, char* szCalib =  NULL, char* szXYZ = NULL ) = 0;

	virtual bool sba_motstr_levmar( ) = 0;

	virtual bool sba_motstr_gn( FixType ft = SBA_FixDefault ) = 0;
};

typedef ISBA ISBAPtr;
SBAapi ISBA*	newSBA();
SBAapi void    freeSBA( ISBA* ptr );

#endif