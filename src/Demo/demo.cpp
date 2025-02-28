#include "stdafx.h"
#include "BAExporter.h"

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
#pragma  comment( lib, "../../bin/BA.lib")

int _tmain(int argc, char* argv[] )
{
	BAExporter BA(sba);
	//BeiJing dataset
	//char* szCam = "../../data/BeiJing/Cam40.txt";
	//char* szFea = "../../data/BeiJing/Feature40.txt";
	//char* szCalib = "../../data/BeiJing/cal40.txt";
	//char* szXYZ = "../../data/BeiJing/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/BeiJing/report.txt";
	//char* szPose = "../../data/BeiJing/FinalPose.txt";
	//char* sz3D = "../../data/BeiJing/Final3D.txt";
	//bool bLM = false;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)       PBA(GN)          SBA(LM)       PBA(LM)
#Iteration          4             4                11             8
Initial MSE      232.408        232.408          232.408       232.408
Final MSE        0.210473       0.210473         0.210473      0.210473
Runtime          0.111          0.135            0.244         0.208
Reason              2             2                 2             2
*/

	//Toronto dataset
	//char* szCam = "../../data/Toronto/Cam13.txt";
	//char* szFea = "../../data/Toronto/Feature13.txt";
	//char* szCalib = "../../data/Toronto/cal13.txt";
	//char* szXYZ = "../../data/Toronto/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/Toronto/report.txt";
	//char* szPose = "../../data/Toronto/FinalPose.txt";
	//char* sz3D = "../../data/Toronto/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 500, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
				   SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            9                      9                    487                   18
Initial MSE    3326138.80221379       3326138.80221379      3326138.80221379      3326138.80221379
Final MSE         0.041856               0.041856              50.083779             0.041856
Runtime           0.717                  0.652                 37.665                1.178
Reason                2                      2                     2                     2
*/

	//WuLong dataset
	char* szCam = "../../data/WuLong/Cam42.txt";
	char* szFea = "../../data/WuLong/Feature42.txt";
	char* szCalib = "../../data/WuLong/cal42.txt";
	char* szXYZ = "../../data/WuLong/XYZ.txt";
	//char* szXYZ = NULL;
	char* szReport = "../../data/WuLong/report.txt";
	char* szPose = "../../data/WuLong/FinalPose.txt";
	char* sz3D = "../../data/WuLong/Final3D.txt";
	bool bLM = true;   //true is Levenberg-Marquardt
	BA.ba_run(false, bLM, 500, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       5                    321                   13
Initial MSE       16.3765                  16.3765               16.3765              16.3765
Final MSE             -                    0.277196              0.277196             0.277196
Runtime               -                    0.139                 6.796                0.295
Reason                -                       2                     2                     2
*/

	//Taian dataset
	//char* szCam = "../../data/Taian/Cam737.txt";
	//char* szFea = "../../data/Taian/Feature737.txt";
	//char* szCalib = "../../data/Taian/cal737.txt";
	//char* szXYZ = "../../data/Taian/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/Taian/report.txt";
	//char* szPose = "../../data/Taian/FinalPose.txt";
	//char* sz3D = "../../data/Taian/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 500, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       11                    427                  61
Initial MSE    747896.41362525        747896.41362525        747896.41362525       747896.41362525
Final MSE             -                    0.038748              5.895712             0.038748
Runtime               -                    20.394                 960.195              123.876
Reason                -                       2                     2                     2
*/

	//DunHuang dataset
	//char* szCam = "../../data/DunHuang/Cam63.txt";
	//char* szFea = "../../data/DunHuang/Feature63.txt";
	//char* szCalib = "../../data/DunHuang/cal63.txt";
	//char* szXYZ = "../../data/DunHuang/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/DunHuang/report.txt";
	//char* szPose = "../../data/DunHuang/FinalPose.txt";
	//char* sz3D = "../../data/DunHuang/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 2000, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            8                       7                    651                  11
Initial MSE    261763.35960011        261763.35960011        261763.35960011       261763.35960011
Final MSE          0.168251                0.168251              0.793898             0.168251
Runtime            1.631                   1.411                 128.398              1.939
Reason                2                       2                     2                     2
*/

	//Vaihingen dataset
	//char* szCam = "../../data/Vaihingen/Cam20.txt";
	//char* szFea = "../../data/Vaihingen/Feature20.txt";
	//char* szCalib = "../../data/Vaihingen/cal20.txt";
	//char* szXYZ = "../../data/Vaihingen/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szPose = "../../data/Vaihingen/FinalPose.txt";
	//char* sz3D = "../../data/Vaihingen/Final3D.txt";
	//char* szReport = "../../data/Vaihingen/report.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 500, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            6                       6                    36                    13
Initial MSE    142089.50626664        142089.50626664        142089.50626664       142089.50626664
Final MSE          0.119365                0.119365              0.119365             0.119365
Runtime            2.560                   2.383                 14.319               4.475
Reason                2                       2                     2                     2
*/

	//College dataset
	//char* szCam = "../../data/College/Cam468.txt";
	//char* szFea = "../../data/College/Feature468.txt";
	//char* szCalib = "../../data/College/cal468.txt";
	//char* szXYZ = "../../data/College/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/College/report.txt";
	//char* szPose = "../../data/College/FinalPose.txt";
	//char* sz3D = "../../data/College/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 1000, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       12                   1001                  17
Initial MSE    202329.44856826         202329.44856826        202329.44856826       202329.44856826
Final MSE             -                   0.734738               16.598062            0.734738
Runtime               -                   12.240                 1765.791             16.395
Reason                -                       2                     3                     2
*/

	//Village dataset
	//char* szCam = "../../data/Village/Cam90.txt";
	//char* szFea = "../../data/Village/Feature90.txt";
	//char* szCalib = "../../data/Village/cal90.txt";
 //   char* szXYZ = "../../data/Village/XYZ.txt";
 //   //char* szXYZ = NULL;
	//char* szReport = "../../data/Village/report.txt";	
 //   char* szPose = "../../data/Village/FinalPose.txt";
 //   char* sz3D = "../../data/Village/Final3D.txt";
	//bool bLM     = false;   //true is Levenberg-Marquardt
 //   BA.ba_run( false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D );

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            6                       6                    15                    11
Initial MSE    28170.98109754           28170.98109754          28170.98109754       28170.98109754
Final MSE          0.083716                 0.083716              0.083716            0.083716
Runtime            1.653                    1.653                 3.751               2.606
Reason                2                       2                     2                     2
*/

	//NewCollege dataset
	//char* szCam = "../../data/NewCollege/Cam3500.txt";
	//char* szFea = "../../data/NewCollege/Feature3500.txt";
	//char* szCalib = "../../data/NewCollege/cal3500.txt";
	//char* szXYZ = "../../data/NewCollege/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/NewCollege/report.txt";
	//char* szPose = "../../data/NewCollege/FinalPose.txt";
	//char* sz3D = "../../data/NewCollege/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 1000, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       -                    1001                  33
Initial MSE      295.79746736           295.79746736           295.79746736          295.79746736
Final MSE             -                       -                  0.991214              1.090911
Runtime               -                       -                  3692.488              92.964
Reason                -                       -                     3                     1
*/


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

