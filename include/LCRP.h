
#ifndef _LCRP_H
#define _LCRP_H

#ifdef LCRP_EXPORTS
	#define LCRPapi  __declspec(dllexport)
#else
	#define LCRPapi __declspec(dllimport)
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

typedef enum LCRP_FixType
{	
	LCRP_FixDefault = -1, //fix the longest axis  
	LCRP_FixX = 0 ,	    // fix X axis	
	LCRP_FixY = 1 ,	    // fix Y axis
	LCRP_FixZ = 2 ,	    // fix Z axis
} FixType ; 
	
class LCRPapi ILCRP
{
public:
	virtual bool lcrp_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 ) = 0;

	virtual	bool lcrp_run( int argc, char** argv ) = 0;

	virtual bool lcrp_initialize( char* szCamera, char* szFeature, char* szCalib =  NULL, char* szXYZ = NULL ) = 0;

	virtual bool lcrp_motstr_levmar( ) = 0;

	virtual bool lcrp_motstr_gn( FixType ft = LCRP_FixDefault ) = 0;
};

typedef ILCRP ILCRPPtr;
LCRPapi ILCRP*	newLCRP();
LCRPapi void    freeLCRP( ILCRP* ptr );

#endif