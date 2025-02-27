//==============================================================================================================
// University of Technology, Sydney, Australia
// 
// Authors:  Liang Zhao         -- Liang.Zhao-1@uts.edu.au 
// 		  Shoudong Huang        -- Shoudong.Huang@uts.edu.au
//        Yanbiao Sun           -- syb51@pku.edu.cn
// 		  Gamini Dissanayake    -- Gamini.Dissanayake@uts.edu.au
// 
// 		  Centre for Autonomous Systems
// 
// 		  Faculty of Engineering and Information Technology
// 
// 		  University of Technology, Sydney
// 
// 		  NSW 2007, Australia
// 
// 		  License
// 
// 		  ParallaxBA by Liang Zhao, Shoudong Huang, Yanbiao Sun, Gamini Dissanayake is licensed under a 
// 
// 		  Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// 
// 		  Please contact Yanbiao Sun {syb51@pku.edu.cn} if you have any questions/comments about the code.
//==============================================================================================================

#include "stdafx.h"
#include "ParallaxBA.h"

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
#pragma  comment( lib, "../../bin/ParallaxBA.lib")
//不需要在 Additional Dependencies 中手动添加 .lib 文件
//也不需要指定 Additional Library Directories，因为编译器会根据 #pragma comment(lib, ...) 自动找到并链接指定的库
//常见方法包括将 DLL 放在与可执行文件相同的目录下，找到dll
int _tmain(int argc, char* argv[] )
{
	IParallaxBA* pBA = newParallaxBA();
	//Toronto dataset
	char* szCam = "C:/txtdata/Toronto/Cam13.txt";
	char* szFea = "C:/txtdata/Toronto/Feature13.txt";
	char* szCalib = "C:/txtdata/Toronto/cal13.txt";
	char* szXYZ = "C:/txtdata/Toronto/xyz_Maxdw.txt";
	//char* szXYZ = NULL;
	char* szReport = "C:/txtdata/Toronto/report.txt";
	char* szPose = "C:/txtdata/Toronto/FinalPose.txt";
	char* sz3D = "C:/txtdata/Toronto/Final3D.txt";
	bool bLM = true;   //true is Levenberg-Marquardt
	pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

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
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

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
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

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
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);
	
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
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

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
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

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
    //pBA->pba_run( false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D );

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
	//pBA->pba_run( false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D );

	freeParallaxBA(pBA);
	system("pause");

	return 0;

/*****************************************************************************************************************************************************/
	//zan shi yong bu dao
	//========================================================================================
	//College dataset
 	  //char*  szCam = "../../../data/College/Cam468.txt";
 	  //char*  szFea = "../../../data/College/Feature468.txt";
 	  //char*  szCalib = "../../../data/College/cal468.txt";
 	  //char*  szReport = "../../../data/College/report.txt";
 	  //bool bLM     = false;
    //pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//NewCollege dataset
 	 // char*  szCam = "../../../data/NewCollege/Cam3500.txt";
 	 // char*  szFea = "../../../data/NewCollege/Feature3500.txt";
 	 // char*  szCalib = "../../../data/NewCollege/cal3500.txt";
 	 // //char*  szXYZ =  "../../../data/NewCollege/XYZ3500.txt";
	  //char* szXYZ = "";
 	 // char*  szReport = "../../../data/NewCollege/report.txt";
 	 // bool bLM     = false;
   //   pBA->pba_run( false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//========================================================================================
	//Venice dataset
    //char*  szCam = "../../../data/Venice/Cam871.txt";
    //char*  szFea = "../../../data/Venice/Feature871.txt";
    //char*  szXYZ =  "../../../data/Venice/XYZ871.txt";
    //char*  szReport = "../../../data/Venice/report.txt";
    //
    //bool bLM     = true;
    //pBA->pba_run( false, bLM, 100, szCam, szFea, szXYZ, NULL, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Toronto dataset
 	  //char*  szCam = "../../../data/Toronto/Cam13.txt";
 	  //char*  szFea = "../../../data/Toronto/Feature13.txt";
 	  //char*  szCalib = "../../../data/Toronto/cal13.txt";
 	  //char*  szReport = "../../../data/Toronto/report.txt";
 	  //bool bLM     = false;   //true is Levenberg-Marquardt
    //pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Vaihingen dataset
    //char*  szCam = "../../../data/Vaihingen/Cam20.txt";
    //char*  szFea = "../../../data/Vaihingen/Feature20.txt";
    //char*  szCalib = "../../../data/Vaihingen/cal20.txt";
    //char*  szReport = "../../../data/Vaihingen/report.txt";
    //
    //bool bLM     = false;   //true is Levenberg-Marquardt
    //pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//DunHuan dataset
//    char*  szCam = "../../../data/DunHuan/Cam63.txt";
//    char*  szFea = "../../../data/DunHuan/Feature63.txt";
//    char*  szCalib = "../../../data/DunHuan/cal63.txt";
//    char*  szReport = "../../../data/DunHuan/report.txt";
//    
//    bool bLM     = false;   //true is Levenberg-Marquardt
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//test1-College
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/College/Cam468.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/College/Feature468.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/College/cal468.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/College/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/College/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/College/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test1-College-LM");
	//bool bLM = true;   //true is Levenberg-Marquardt
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test1-College-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test1-College-GN");
	//bLM = false;
	//t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test1-College-GN-包含初始化内存的平差时间：", t_diff);

	//test2-DunHuan
	//char* szCam = "C:/Users/lwyan/Desktop/ParallaxBA_2012/data_test20241114/DunHuan/Cam63.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/ParallaxBA_2012/data_test20241114/DunHuan/Feature63.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/ParallaxBA_2012/data_test20241114/DunHuan/cal63.txt";
	//char* szXYZ = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/DunHuan/output_ssba_4_no3D_v1_pba.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/ParallaxBA_2012/data_test20241114/DunHuan/report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/ParallaxBA_2012/data_test20241114/DunHuan/OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/ParallaxBA_2012/data_test20241114/DunHuan/OutputOptimal3DPts.txt";
	//printf("%s\n", "test2-DunHuan-LM");
	//bool bLM = true;
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 1000, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test2-DunHuan-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test2-DunHuan-GN");
	//bLM = false;
	//t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test2-DunHuan-GN-包含初始化内存的平差时间：", t_diff);

	//test3-Jinan
	//char* szCam =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Jinan/Cam76.txt";
	//char* szFea =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Jinan/Feature76.txt";
	//char* szCalib =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Jinan/cal76.txt";
	//char* szReport =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Jinan/LM_report.txt";
	//char* szOutputOptimalCamera =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Jinan/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Jinan/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test3-Jinan-LM");
	//bool bLM = true;
	//double t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test3-Jinan-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test3-Jinan-GN");
	//bLM = false;
	//t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test3-Jinan-GN-包含初始化内存的平差时间：", t_diff);

	//test4-Malaga
	//char* szCam =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Malaga/Cam170.txt";
	//char* szFea =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Malaga/Feature170.txt";
	//char* szCalib =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Malaga/cal170.txt";
	//char* szReport =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Malaga/LM_report.txt";
	//char* szOutputOptimalCamera =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Malaga/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Malaga/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test4-Malaga-LM");
	//bool bLM = true;
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test4-Malaga-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test4-Malaga-GN");
	//bLM = false;
	//t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test4-Malaga-GN-包含初始化内存的平差时间：", t_diff);

	//test5-NewCollege
	//char* szCam = "C:/NewCollege_temp/Cam3500.txt";
	//char* szFea = "C:/NewCollege_temp/Feature3500.txt";
	//char* szCalib = "C:/NewCollege_temp/cal3500.txt";
	//char* szReport = "C:/NewCollege_temp/LM_report.txt";
	//char* szXYZ = "C:/NewCollege_temp/Triangulatedxyz.txt";xyz_Maxdw.txt
	//char* szXYZ = "C:/NewCollege_temp/xyz_Maxdw.txt"; 
	//char* szXYZ = NULL;
	//char* szOutputOptimalCamera = "C:/NewCollege_temp/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/NewCollege_temp/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test5-NewCollege-LM");
	//bool bLM = true;
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 300, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test5-NewCollege-LM-包含初始化内存的平差时间：", t_diff);
	//
	//printf("%s\n", "test5-NewCollege-GN");
	//bLM = false;
	//t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test5-NewCollege-GN-包含初始化内存的平差时间：", t_diff);




	//test6-Taian
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Taian/Camera737.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Taian/Feature737.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Taian/cal737.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Taian/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Taian/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Taian/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test6-Taian-LM");
	//bool bLM = false;
	//double t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL);
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test6-Taian-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test6-Taian-GN");
	//bLM = true;
	//t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test6-Taian-GN-包含初始化内存的平差时间：", t_diff);

	//test7-Toronto
	//char* szCam =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Toronto/Cam13.txt";
	//char* szFea =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Toronto/Feature13.txt";
	//char* szCalib =		"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Toronto/cal13.txt";
	//char* szReport =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Toronto/LM_report.txt";
	//char* szOutputOptimalCamera =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Toronto/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts =	"C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Toronto/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test7-Toronto-LM");
	//bool bLM = true;
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test7-Toronto-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test7-Toronto-GN");
	//bLM = false;
	//t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test7-Toronto-GN-包含初始化内存的平差时间：", t_diff);

	//test8-Vaihingen
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Vaihingen/Cam20.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Vaihingen/Feature20.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Vaihingen/cal20.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Vaihingen/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Vaihingen/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Vaihingen/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test8-Vaihingen-LM");
	//bool bLM = true;
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test8-Vaihingen-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test8-Vaihingen-GN");
	//bLM = false;
	//t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test8-Vaihingen-GN-包含初始化内存的平差时间：", t_diff);

	//test10-Village
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Village/Cam90.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Village/Feature90.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Village/cal90.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Village/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Village/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Village/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "test10-Village-LM");
	//bool bLM = true;				//true is Levenberg-Marquardt
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test10-Village-LM-包含初始化内存的平差时间：", t_diff);
	//printf("%s\n", "test10-Village-GN");
	//bLM = false;
	//t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//t2 = clock();
	//t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "test10-Village-GN-包含初始化内存的平差时间：", t_diff);

	//test9-Venice
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Venice_incomplete/Cam871.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Venice_incomplete/Feature871.txt";
	////char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Venice_incomplete/cal.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Venice_incomplete/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Venice_incomplete/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/Venice_incomplete/LM_OutputOptimal3DPts.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//double t1 = clock();//56
	//char* szCalib = NULL;
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);//56 maximum number of iterations
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "包含初始化内存的平差时间：", t_diff);
	//========================================================================================
	//Venice dataset
	//char*  szCam = "../../../data/Venice/Cam871.txt";
	//char*  szFea = "../../../data/Venice/Feature871.txt";
	//char*  szXYZ =  "../../../data/Venice/XYZ871.txt";
	//char*  szReport = "../../../data/Venice/report.txt";
	//
	//bool bLM     = true;
	//pBA->pba_run( false, bLM, 100, szCam, szFea, szXYZ, NULL, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//test-5-1
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/tuniu_8_98_sba/cam.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/tuniu_8_98_sba/projs.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/tuniu_8_98_sba/cal.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/tuniu_8_98_sba/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/tuniu_8_98_sba/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/tuniu_8_98_sba/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "tuniu_8_98_sba-LM");
	//bool bLM = true;				//true is Levenberg-Marquardt
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "tuniu_8_98_sba-LM-包含初始化内存的平差时间：", t_diff);
//test-5-2
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/cug_8_192_sba/cam.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/cug_8_192_sba/projs.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/cug_8_192_sba/cal.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/cug_8_192_sba/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/cug_8_192_sba/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/cug_8_192_sba/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "cug_8_192_sba-LM");
	//bool bLM = true;				//true is Levenberg-Marquardt
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "cug_8_192_sba-LM-包含初始化内存的平差时间：", t_diff);
//test-5-3
	//char* szCam = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/beijin-sba/cam.txt";
	//char* szFea = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/beijin-sba/projs.txt";
	//char* szCalib = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/beijin-sba/cal.txt";
	//char* szReport = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/beijin-sba/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/beijin-sba/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Users/lwyan/Desktop/Todo/ParallaxBA_all_results/data_test20241114/5_datasets_20241222_test/beijin-sba/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "beijin-sba-LM");
	//bool bLM = true;				//true is Levenberg-Marquardt
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "beijin-sba-LM-包含初始化内存的平差时间：", t_diff);

//beijin数据集
	//char* szCam = "C:/beijin_temp/cam.txt";
	//char* szFea = "C:/beijin_temp/projs.txt";
	//char* szCalib = "C:/beijin_temp/cal.txt";
	////char* szXYZ = "C:/beijin_temp/Initial3DPts.txt";
	//char* szXYZ = "C:/beijin_temp/LM_OutputOptimal3DPts.txt";
	//char* szReport = "C:/beijin_temp/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/beijin_temp/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/beijin_temp/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "beijin-pba-LM");
	//bool bLM = true;				//true is Levenberg-Marquardt
	//double t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "beijin-pba-LM-包含初始化内存的平差时间：", t_diff);

	////char* szCam = "G:/testWulong42/pbassba/Cam42.txt";
	////char* szFea = "G:/testWulong42/pbassba/Feature42.txt";
	////char* szCalib = "G:/testWulong42/pbassba/cal42.txt";
	////char* szReport = "G:/testWulong42/pbassba/LM_report.txt";
	////char* szXYZ = "G:/testWulong42/pbassba/LM_OutputInitial3DPts_maxdw.txt";
	//////char* szXYZ = "C:/NewCollege_temp/xyz_Maxdw.txt"; 
	//////char* szXYZ = NULL;
	////char* szOutputOptimalCamera = "G:/testWulong42/pbassba/LM_OutputOptimalCamera.txt";
	////char* szOutputOptimal3DPts = "G:/testWulong42/pbassba/LM_OutputOptimal3DPts.txt";
	////printf("%s\n", "test5-Wulong-LM");
	////bool bLM = true;
	////double t1 = clock();
	////pBA->pba_run(false, bLM, 300, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	////double t2 = clock();
	////double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	////printf("%s %f\n", "test5-Wulong-LM-包含初始化内存的平差时间：", t_diff);

//Taian数据集
	//char* szCam = "C:/Taian_temp/Cam737.txt";
	//char* szFea = "C:/Taian_temp/Feature737.txt";
	//char* szCalib = "C:/Taian_temp/cal737.txt";
	////char* szXYZ = "C:/Taian_temp/Initial3DPts.txt";
	////char* szXYZ = "C:/Taian_temp/LM_OutputOptimal3DPts.txt";
	//char* szXYZ = NULL;
	//char* szReport = "C:/Taian_temp/LM_report.txt";
	//char* szOutputOptimalCamera = "C:/Taian_temp/LM_OutputOptimalCamera.txt";
	//char* szOutputOptimal3DPts = "C:/Taian_temp/LM_OutputOptimal3DPts.txt";
	//printf("%s\n", "Taian-pba-LM");
	//bool bLM = false;				//true is Levenberg-Marquardt
	//double t1 = clock();
	////pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//pBA->pba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "Taian-pba-LM-包含初始化内存的平差时间：", t_diff);
	//========================================================================================
	//Jinan dataset
//    char*  szCam = "../../../data/Jinan/Cam76.txt";
//    char*  szFea = "../../../data/Jinan/Feature76.txt";
//    char*  szCalib = "../../../data/Jinan/cal76.txt";
//    char*  szReport = "../../../data/Jinan/report.txt";
//     
//    bool bLM     = true;   //true is Levenberg-Marquardt
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL, 1E-3 );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//========================================================================================
	//Dataset4 Xifengshan
	//char*  szCam = "../../../data/Dataset4/low_accuracy/Cam100.txt";
	//char*  szFea = "../../../data/Dataset4/low_accuracy/Feature100.txt";
	//char*  szCalib = "../../../data/Dataset4/low_accuracy/cal.txt";
	//char*  szReport = "../../../data/Dataset4/low_accuracy/report.txt";
	//char*  szOutputOptimalCamera = "../../../data/Dataset4/OutputOptimalCamera.txt";
	//char*  szOutputOptimal3DPts = "../../../data/Dataset4/OutputOptimal3DPts.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//double t1 = clock();
	//pBA->pba_run(false, bLM, 30, szCam, szFea, NULL, szCalib, szReport, szOutputOptimalCamera, szOutputOptimal3DPts);//30为最大迭代次数
	//double t2 = clock();
	//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
	//printf("%s %f\n", "包含初始化内存的平差时间：", t_diff);
}

