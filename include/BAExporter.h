#ifndef _BA_H
#define _BA_H

#ifdef BA_EXPORTS
	#define BAapi  __declspec(dllexport)
#else
	#define BAapi __declspec(dllimport)
#endif

typedef enum BA_FixType
{	
	BA_FixDefault = -1, //fix the longest axis  
	BA_FixX = 0 ,	    // fix X axis	
	BA_FixY = 1 ,	    // fix Y axis
	BA_FixZ = 2 ,	    // fix Z axis
} FixType ; 

typedef enum BA_Type
{
	pba = 1,
	sba = 2,
}BAType;

class IBA;//PBA和SBA的基类
class BAapi BAExporter
{
public:
	bool ba_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 );

	bool ba_run( int argc, char** argv );

	bool ba_initialize( char* szCamera, char* szFeature, char* szCalib =  NULL, char* szXYZ = NULL );

	bool ba_motstr_levmar( );

	bool ba_motstr_gn( FixType ft = BA_FixDefault );

	BAExporter(BAType);
	~BAExporter(void);

private:
	IBA* ptr;
};

#endif