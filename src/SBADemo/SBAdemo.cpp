#include "stdafx.h"
#include "SBA.h"

#include <vector>
#include <map>
using namespace std;

#pragma  comment( lib, "../../bin/libamd.lib")
#pragma  comment( lib, "../../bin/libcamd.lib")	
#pragma  comment( lib, "../../bin/libccolamd.lib")
#pragma  comment( lib, "../../bin/libcholmod.lib")
#pragma  comment( lib, "../../bin/libcolamd.lib")
#pragma  comment( lib, "../../bin/libmetis_CHOLMOD.lib")
#pragma  comment( lib, "../../bin/libgoto2.lib")
#pragma  comment( lib, "../../bin/SBA.lib")

int _tmain(int argc, char* argv[] )
{
	ISBA* sBA = newSBA();
	//Toronto dataset
	char* szCam = "C:/txtdata/Toronto/Cam13.txt";
	char* szFea = "C:/txtdata/Toronto/Feature13.txt";
	char* szCalib = "C:/txtdata/Toronto/cal13.txt";
	char* szXYZ = "C:/txtdata/Toronto/xyz_Maxdw.txt";
	//char* szXYZ = NULL;
	char* szReport = "C:/txtdata/Toronto/report.txt";
	char* szPose = "C:/txtdata/Toronto/FinalPose.txt";
	char* sz3D = "C:/txtdata/Toronto/Final3D.txt";
	bool bLM = false;   //true is Levenberg-Marquardt
	sBA->sba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

	//WuLong dataset
	//char* szCam = "C:/txtdata/WuLong/Cam42.txt";
	//char* szFea = "C:/txtdata/WuLong/Feature42.txt";
	//char* szCalib = "C:/txtdata/WuLong/cal42.txt";
	//char* szXYZ = "C:/txtdata/WuLong/xyz_Maxdw.txt";
	////char* szXYZ = NULL;
	//char* szReport = "C:/txtdata/WuLong/report.txt";
	//char* szPose = "C:/txtdata/WuLong/FinalPose.txt";
	//char* sz3D = "C:/txtdata/WuLong/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//sBA->sba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

	//Taian dataset
	//char* szCam = "C:/txtdata/Taian/Cam737.txt";
	//char* szFea = "C:/txtdata/Taian/Feature737.txt";
	//char* szCalib = "C:/txtdata/Taian/cal737.txt";
	//char* szXYZ = "C:/txtdata/Taian/xyz_Maxdw.txt";
	////char* szXYZ = NULL;
	//char* szReport = "C:/txtdata/Taian/report.txt";
	//char* szPose = "C:/txtdata/Taian/FinalPose.txt";
	//char* sz3D = "C:/txtdata/Taian/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//sBA->sba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

	//DunHuan dataset
	//char* szCam = "C:/txtdata/DunHuan/Cam63.txt";
	//char* szFea = "C:/txtdata/DunHuan/Feature63.txt";
	//char* szCalib = "C:/txtdata/DunHuan/cal63.txt";
	//char* szXYZ = "C:/txtdata/DunHuan/xyz_Maxdw.txt";
	////char* szXYZ = NULL;
	//char* szReport = "C:/txtdata/DunHuan/report.txt";
	//char* szPose = "C:/txtdata/DunHuan/FinalPose.txt";
	//char* sz3D = "C:/txtdata/DunHuan/Final3D.txt";
	//bool bLM = false;   //true is Levenberg-Marquardt
	//sBA->sba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

	//Vaihingen dataset
	//char* szCam = "C:/txtdata/Vaihingen/Cam20.txt";
	//char* szFea = "C:/txtdata/Vaihingen/Feature20.txt";
	//char* szCalib = "C:/txtdata/Vaihingen/cal20.txt";
	//char* szXYZ = "C:/txtdata/Vaihingen/xyz_Maxdw.txt";
	////char* szXYZ = NULL;
	//char* szPose = "C:/txtdata/Vaihingen/FinalPose.txt";
	//char* sz3D = "C:/txtdata/Vaihingen/Final3D.txt";
	//char* szReport = "C:/txtdata/Vaihingen/report.txt";
	//bool bLM = false;   //true is Levenberg-Marquardt
	//sBA->sba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

	//College dataset
	//char* szCam = "C:/txtdata/College/Cam468.txt";
	//char* szFea = "C:/txtdata/College/Feature468.txt";
	//char* szCalib = "C:/txtdata/College/cal468.txt";
	////char* szXYZ = "C:/txtdata/College/xyz_Maxdw.txt";
	//char* szXYZ = NULL;
	//char* szReport = "C:/txtdata/College/report.txt";
	//char* szPose = "C:/txtdata/College/FinalPose.txt";
	//char* sz3D = "C:/txtdata/College/Final3D.txt";
	//bool bLM = false;   //true is Levenberg-Marquardt
	//sBA->sba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

	//Village dataset
	//char* szCam = "C:/txtdata/Village/Cam90.txt";
	//char* szFea = "C:/txtdata/Village/Feature90.txt";
	//char* szCalib = "C:/txtdata/Village/cal90.txt";
	//char* szXYZ = "C:/txtdata/Village/xyz_Maxdw.txt";
	////char* szXYZ = NULL;
	//char* szReport = "C:/txtdata/Village/report.txt";	
	//char* szPose = "C:/txtdata/Village/FinalPose.txt";
	//char* sz3D = "C:/txtdata/Village/Final3D.txt";
	//bool bLM     = false;   //true is Levenberg-Marquardt
	//sBA->sba_run( false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D );

	//NewCollege dataset
	//char* szCam = "C:/txtdata/NewCollege/Cam3500.txt";
	//char* szFea = "C:/txtdata/NewCollege/Feature3500.txt";
	//char* szCalib = "C:/txtdata/NewCollege/cal3500.txt";
	//char* szXYZ = "C:/txtdata/NewCollege/xyz_Maxdw.txt";
	////char* szXYZ = NULL;
	//char* szReport = "C:/txtdata/NewCollege/report.txt";
	//char* szPose = "C:/txtdata/NewCollege/FinalPose.txt";
	//char* sz3D = "C:/txtdata/NewCollege/Final3D.txt";	
	//bool bLM     = true;   //true is Levenberg-Marquardt
	//sBA->sba_run( false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D );

	freeSBA(sBA);
	system( "pause" );

	return 0;
}

