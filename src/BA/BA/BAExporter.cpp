#include "stdafx.h"
#include "BAExporter.h"
#include "PBAImp.h"
#include "SBAImp.h"

BAExporter::BAExporter(BAType ba)
{
	if (ba == pba)
		ptr = new PBA;
	else if (ba == sba)
		ptr = new SBA;
}

BAExporter::~BAExporter(void)
{
	delete ptr;
}

bool BAExporter::ba_run(bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ, char* szCalib, char* szReport,
	char* szPose, char* sz3D, double Tau)
{
	return ptr->ba_run(bRobust, bLM, nMaxIter, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D, Tau);
}
bool BAExporter::ba_run(int argc, char** argv)
{
	return ptr->ba_run(argc, argv);
}
bool BAExporter::ba_initialize(char* szCamera, char* szFeature, char* szCalib, char* szXYZ)
{
	return	ptr->ba_initialize(szCamera, szFeature, szCalib, szXYZ);
}
bool BAExporter::ba_motstr_levmar()
{
	return ptr->ba_motstr_levmar();
}
bool BAExporter::ba_motstr_gn(FixType ft)
{
	return ptr->ba_motstr_gn(ft);
}