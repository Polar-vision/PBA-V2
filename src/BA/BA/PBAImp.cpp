#include "stdafx.h"
#include "PBAImp.h"

//compute reprojection image coordinates of each image points
static void pba_reprojectEachPts( double *KR, double* pa, double* pb, int nM, int nN, int nP, double n[2] )
{
	double *posei, *posek;
	double pti2k[3];
	double ptXUnit[3];
	double dDot, dDisi2k, dW2;
	double ptXk[3];
	double *posel;
	double pti2l[3];
	double ptXj[3];
	double *pKR;

	if ( nP == nM )	//Main archor
	{
		ptXj[0] = sin( pb[0] ) * cos( pb[1] );
		ptXj[1] = sin( pb[1] );
		ptXj[2] = cos( pb[0] ) * cos( pb[1] );

		pKR = KR + nP*9;
		n[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
			(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

		n[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
			(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);
	}else 
		if ( nP == nN )	//Associate archor
		{
			posei = pa+6*nM+3;
			posek = pa+6*nN+3;

			//主锚点到副锚点的平移向量
			pti2k[0] = posek[0]-posei[0];	pti2k[1] = posek[1]-posei[1];	pti2k[2] = posek[2]-posei[2];	

			//主锚点到特征点的单位向量
			ptXUnit[0] = sin( pb[0] ) * cos( pb[1] );
			ptXUnit[1] = sin( pb[1] );
			ptXUnit[2] = cos( pb[0] ) * cos( pb[1] );

			//compute angle w2
			dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
			dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
			if (dDot/dDisi2k > 1)
				dW2 = 0;
			if ( dDot/dDisi2k < -1)
				dW2 = PI;
			else
				dW2  = acos( dDot/dDisi2k );

			//compute Xk vector according sin theory
			ptXk[0] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[0] - sin( pb[2] ) * pti2k[0];
			ptXk[1] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[1] - sin( pb[2] ) * pti2k[1];
			ptXk[2] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[2] - sin( pb[2] ) * pti2k[2];

			pKR = KR + nN*9;		
			n[0] = ( pKR[0]*ptXk[0]	+ pKR[1]*ptXk[1] + pKR[2]*ptXk[2])/
				( pKR[6]*ptXk[0] + pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);

			n[1] = ( pKR[3]*ptXk[0]	+ pKR[4]*ptXk[1] + pKR[5]*ptXk[2])/
				(pKR[6]*ptXk[0]	+ pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);
		}
		else
		{
			posel = pa + nP*6 + 3;
			posek = pa + nN *6 + 3;
			posei = pa + nM *6 + 3;

			//主锚点到副锚点的平移向量
			pti2k[0] = posek[0] - posei[0];		pti2k[1] = posek[1] - posei[1];		pti2k[2] = posek[2] - posei[2];
			//主锚点到当且锚点的平移向量
			pti2l[0] = posel[0] - posei[0];		pti2l[1] = posel[1] - posei[1];		pti2l[2] = posel[2] - posei[2];
			
			//XUnit 主锚点到特征点的单位向量
			ptXUnit[0] = sin( pb[0] ) * cos( pb[1] );
			ptXUnit[1] = sin( pb[1] );
			ptXUnit[2] = cos( pb[0] ) * cos( pb[1] );

			//compute angle w2
			dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
			dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
			//dW2  = acos( dDot/dDisi2k );
			if (dDot/dDisi2k > 1)
				dW2 = 0;
			if ( dDot/dDisi2k < -1)
				dW2 = PI;
			else
				dW2  = acos( dDot/dDisi2k );

			//compute Xl vector according sin theory
			ptXk[0] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[0] - sin( pb[2] ) * pti2l[0];
			ptXk[1] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[1] - sin( pb[2] ) * pti2l[1];
			ptXk[2] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[2] - sin( pb[2] ) * pti2l[2];

			pKR = KR + nP*9;			
			n[0] = (pKR[0]*ptXk[0] + pKR[1]*ptXk[1] + pKR[2]*ptXk[2])/
				( pKR[6]*ptXk[0] + pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);

			n[1] = (pKR[3]*ptXk[0] + pKR[4]*ptXk[1] + pKR[5]*ptXk[2])/
				( pKR[6]*ptXk[0] + pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);
		}
}
//#pragma optimize("", off) // 关闭优化
//compute Jacobian for each image point
//pPA pPB,uv对相机参数和三维点的一阶导
//pAM uv对主锚点相机位移参数的一阶导
//pAA uv对副锚点相机位移参数的一阶导
static void pba_jacobianEachPts(double* KR, double* KdA, double *KdB, double* KdG, double* pa, double* ppt, int nM, int nN, int nP,double* pAM, double* pAA, double* pPA, double* pPB)
{
	double matXj[3];
	double matxyt[3];
	double matDuvDxyt[6];
	double matDxytDRA[3], matDxytDRB[3], matDxytDRG[3];
	double matDuvDRA[2], matDuvDRB[2], matDuvDRG[2];
	double matDXjDDH[6];
	double matDxytDDH[6];
	double matDuvDDH[4];
	double *ptAngle = NULL;

	double *posei, *posek, *posel;
	double matPosei2k[3], matPosei2l[3];
	double	dDisTi2k, dDot, dArcCosIn, dW;

	double matXk[3],matXl[3];
	double matDDotDDH[2];
	double matDArcCosInDDH[2];
	double dDWDArcCosIn;
	double matDsinwWDDH[2];
	double matDXkDpa[3];
	double matDuvDpa[2];
	double tmp1[6];
	double tmp2[6];
	
	double matDdisDTi[3];
	double matDDotDTi[3];
	double matDArcCosInDTi[3];
	double matDsinwWDTi[3];
	double matDTi2kDTi[9], matDTi2kDTk[9];
	double matDXkDTi[9], matDXkDTk[9], matDXlDTl[9];
	double matDuvDTi[6], matDuvDTk[6], matDuvDTl[6];
	double matDdisDTk[3];
	double matDDotDTk[3], matDArcCosInDTk[3], matDsinwWDTk[3];
	double *pKR, *pKdA, *pKdB, *pKdG;

	if ( nP == nM )	//main archor
	{
		matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
		matXj[1] = sin( ppt[1] );
		matXj[2] = cos( ppt[0] ) * cos( ppt[1] );

		pKR = KR + nP*9;
		matxyt[0] = pKR[0]*matXj[0] +pKR[1]*matXj[1] +pKR[2]*matXj[2];
		matxyt[1] = pKR[3]*matXj[0] +pKR[4]*matXj[1] +pKR[5]*matXj[2];
		matxyt[2] = pKR[6]*matXj[0] +pKR[7]*matXj[1] +pKR[8]*matXj[2];

		matDuvDxyt[0] = 1/matxyt[2];	
		matDuvDxyt[1] = 0;
		matDuvDxyt[2] = -matxyt[0]/(matxyt[2]*matxyt[2]);
		matDuvDxyt[3] = 0;	
		matDuvDxyt[4] = 1/matxyt[2];		
		matDuvDxyt[5] = -matxyt[1]/(matxyt[2]*matxyt[2]);
		
		//pKdG:pKR矩阵关于omega的一阶导
		//pKdB:pKR矩阵关于phi的一阶导
		//pKdA:pKR矩阵关于kappa的一阶导
		//camera angles
		pKdG = KdG + nP*9;
		matDxytDRG[0] =  pKdG[0]*matXj[0] + pKdG[1]*matXj[1] + pKdG[2]*matXj[2];//xyt关于omega的一阶导
		matDxytDRG[1] =  pKdG[3]*matXj[0] + pKdG[4]*matXj[1] + pKdG[5]*matXj[2];
		matDxytDRG[2] =  pKdG[6]*matXj[0] + pKdG[7]*matXj[1] + pKdG[8]*matXj[2];

		pKdB = KdB + nP*9;
		matDxytDRB[0] =  pKdB[0]*matXj[0] + pKdB[1]*matXj[1] + pKdB[2]*matXj[2];//xyt关于phi的一阶导
		matDxytDRB[1] =  pKdB[3]*matXj[0] + pKdB[4]*matXj[1] + pKdB[5]*matXj[2];
		matDxytDRB[2] =  pKdB[6]*matXj[0] + pKdB[7]*matXj[1] + pKdB[8]*matXj[2];

		pKdA = KdA + nP*9;
		matDxytDRA[0] =  pKdA[0]*matXj[0] + pKdA[1]*matXj[1] + pKdA[2]*matXj[2];//xyt关于kappa的一阶导
		matDxytDRA[1] =  pKdA[3]*matXj[0] + pKdA[4]*matXj[1] + pKdA[5]*matXj[2];
		matDxytDRA[2] =  pKdA[6]*matXj[0] + pKdA[7]*matXj[1] + pKdA[8]*matXj[2];

		matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv关于kappa的一阶导
		matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];

		matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv关于phi的一阶导
		matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];

		matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv关于omega的一阶导
		matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];

		pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv对主锚点旋转参数的一阶导
		pPA[3] = 0;						pPA[4] = 0;						pPA[5] = 0;
		pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
		pPA[9] = 0;						pPA[10] = 0;					pPA[11] = 0;
		
		//matXj[0] = sin(ppt[0]) * cos(ppt[1]);
		//matXj[1] = sin(ppt[1]);
		//matXj[2] = cos(ppt[0]) * cos(ppt[1]);
		//matxyt[0] = pKR[0] * matXj[0] + pKR[1] * matXj[1] + pKR[2] * matXj[2];
		//matxyt[1] = pKR[3] * matXj[0] + pKR[4] * matXj[1] + pKR[5] * matXj[2];
		//matxyt[2] = pKR[6] * matXj[0] + pKR[7] * matXj[1] + pKR[8] * matXj[2];
		//azimuth and elevation angle
		matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
		matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
		matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);

		matDxytDDH[0] = pKR[0]*matDXjDDH[0]	+ pKR[1]*matDXjDDH[2] + pKR[2]*matDXjDDH[4];//x关于azimuth的导数
		matDxytDDH[1] = pKR[0]*matDXjDDH[1] + pKR[1]*matDXjDDH[3] + pKR[2]*matDXjDDH[5];//x关于elevation的导数

		matDxytDDH[2] = pKR[3]*matDXjDDH[0]	+ pKR[4]*matDXjDDH[2] + pKR[5]*matDXjDDH[4];//y关于azimuth的导数
		matDxytDDH[3] = pKR[3]*matDXjDDH[1] + pKR[4]*matDXjDDH[3] + pKR[5]*matDXjDDH[5];//y关于elevation的导数

		matDxytDDH[4] = pKR[6]*matDXjDDH[0]	+ pKR[7]*matDXjDDH[2] + pKR[8]*matDXjDDH[4];//t关于azimuth的导数
		matDxytDDH[5] = pKR[6]*matDXjDDH[1] + pKR[7]*matDXjDDH[3] + pKR[8]*matDXjDDH[5];//t关于elevation的导数

		matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//uv关于azimuth的导数
		matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];
		matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//uv关于elevation的导数
		matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];

		pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = 0;//uv对azimuth和elevation的一阶导
		pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = 0;
		
	}else 
		if ( nP == nN )	//associate archor
		{			
			ptAngle = pa + nN*6;
			matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
			matXj[1] = sin( ppt[1] );
			matXj[2] = cos( ppt[0] ) * cos( ppt[1] );
			
			pKR = KR + nN*9;
			pKdA= KdA+ nN*9;
			pKdB= KdB+ nN*9;
			pKdG= KdG+ nN*9;
			
			posei = pa + nM*6 + 3;//主锚点Xc
			posek = pa + nN*6 + 3;//副锚点Xc
			matPosei2k[0] = posek[0]-posei[0];		matPosei2k[1] = posek[1]-posei[1];		matPosei2k[2] = posek[2]-posei[2];//主锚点到副锚点的向量

			dDisTi2k = sqrt( matPosei2k[0]*matPosei2k[0] + matPosei2k[1]*matPosei2k[1] + matPosei2k[2]*matPosei2k[2] );//主锚点到副锚点的向量的模长
			dDot     = matXj[0]*matPosei2k[0] + matXj[1]*matPosei2k[1] + matXj[2]*matPosei2k[2];
			dArcCosIn= dDot / dDisTi2k;
			dW       = acos( dArcCosIn );

			matXk[0] = dDisTi2k * sin( dW + ppt[2] ) * matXj[0] - sin(ppt[2])*matPosei2k[0];
			matXk[1] = dDisTi2k * sin( dW + ppt[2] ) * matXj[1] - sin(ppt[2])*matPosei2k[1];
			matXk[2] = dDisTi2k * sin( dW + ppt[2] ) * matXj[2] - sin(ppt[2])*matPosei2k[2];

			matxyt[0] = pKR[0]*matXk[0] + pKR[1]*matXk[1] + pKR[2]*matXk[2];
			matxyt[1] = pKR[3]*matXk[0] + pKR[4]*matXk[1] + pKR[5]*matXk[2];
			matxyt[2] = pKR[6]*matXk[0] + pKR[7]*matXk[1] + pKR[8]*matXk[2];

			matDuvDxyt[0] = 1/matxyt[2];	
			matDuvDxyt[1] = 0;
			matDuvDxyt[2] = -matxyt[0]/(matxyt[2]*matxyt[2]);
			matDuvDxyt[3] = 0;	
			matDuvDxyt[4] = 1/matxyt[2];		
			matDuvDxyt[5] = -matxyt[1]/(matxyt[2]*matxyt[2]);
			
//pKdG:pKR矩阵关于omega的一阶导
//pKdB:pKR矩阵关于phi的一阶导
//pKdA:pKR矩阵关于kappa的一阶导
			//camera angles
			matDxytDRG[0] =  pKdG[0]*matXk[0] + pKdG[1]*matXk[1] + pKdG[2]*matXk[2];
			matDxytDRG[1] =  pKdG[3]*matXk[0] + pKdG[4]*matXk[1] + pKdG[5]*matXk[2];
			matDxytDRG[2] =  pKdG[6]*matXk[0] + pKdG[7]*matXk[1] + pKdG[8]*matXk[2];

			
			matDxytDRB[0] =  pKdB[0]*matXk[0] + pKdB[1]*matXk[1] + pKdB[2]*matXk[2];
			matDxytDRB[1] =  pKdB[3]*matXk[0] + pKdB[4]*matXk[1] + pKdB[5]*matXk[2];
			matDxytDRB[2] =  pKdB[6]*matXk[0] + pKdB[7]*matXk[1] + pKdB[8]*matXk[2];

			matDxytDRA[0] =  pKdA[0]*matXk[0] + pKdA[1]*matXk[1] + pKdA[2]*matXk[2];
			matDxytDRA[1] =  pKdA[3]*matXk[0] + pKdA[4]*matXk[1] + pKdA[5]*matXk[2];
			matDxytDRA[2] =  pKdA[6]*matXk[0] + pKdA[7]*matXk[1] + pKdA[8]*matXk[2];

			matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv关于kappa的一阶导
			matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];

			matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv关于phi的一阶导
			matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];

			matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv关于omega的一阶导
			matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];
			
			//Xj关于azimuth和elevation的一阶导
			//azimuth and elevation angle
			matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
			matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
			matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);

//matPosei2k[0] = posek[0] - posei[0];		matPosei2k[1] = posek[1] - posei[1];		matPosei2k[2] = posek[2] - posei[2];//主锚点到副锚点的向量

//dDisTi2k = sqrt(matPosei2k[0] * matPosei2k[0] + matPosei2k[1] * matPosei2k[1] + matPosei2k[2] * matPosei2k[2]);//主锚点到副锚点的向量的模长
//dDot = matXj[0] * matPosei2k[0] + matXj[1] * matPosei2k[1] + matXj[2] * matPosei2k[2];
//dArcCosIn = dDot / dDisTi2k;
//dW = acos(dArcCosIn);

			matDDotDDH[0] = cos(ppt[0])*cos(ppt[1])*matPosei2k[0] - sin(ppt[0])*cos(ppt[1])*matPosei2k[2];//dDot关于azimuth的一阶导
			matDDotDDH[1] = -sin(ppt[0])*sin(ppt[1])*matPosei2k[0] + cos(ppt[1])*matPosei2k[1]-cos(ppt[0])*sin(ppt[1])*matPosei2k[2];//dDot关于elevation的一阶导

			matDArcCosInDDH[0] = matDDotDDH[0]/dDisTi2k;//ArcCosIn关于azimuth的一阶导
			matDArcCosInDDH[1] = matDDotDDH[1]/dDisTi2k;//ArcCosIn关于elevation的一阶导

			dDWDArcCosIn = -1/sqrt(1-dArcCosIn*dArcCosIn);//dW关于ArccosIn的一阶导

//matXk[0] = dDisTi2k * sin(dW + ppt[2]) * matXj[0] - sin(ppt[2]) * matPosei2k[0];
//matXk[1] = dDisTi2k * sin(dW + ppt[2]) * matXj[1] - sin(ppt[2]) * matPosei2k[1];
//matXk[2] = dDisTi2k * sin(dW + ppt[2]) * matXj[2] - sin(ppt[2]) * matPosei2k[2];

			matDsinwWDDH[0] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[0];//sin(dW + ppt[2])关于azimuth的一阶导
			matDsinwWDDH[1] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[1];//sin(dW + ppt[2])关于elevation的一阶导

			tmp2[0] =  dDisTi2k*matXj[0]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[0];//matXk[0]关于azimuth的一阶导
			tmp2[1] =  dDisTi2k*matXj[0]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[1];//matXk[0]关于elevation的一阶导
			tmp2[2] =  dDisTi2k*matXj[1]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[2];//matXk[1]关于azimuth的一阶导
			tmp2[3] =  dDisTi2k*matXj[1]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[3];//matXk[1]关于elevation的一阶导
			tmp2[4] =  dDisTi2k*matXj[2]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[4];//matXk[2]关于azimuth的一阶导
			tmp2[5] =  dDisTi2k*matXj[2]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[5];//matXk[2]关于elevation的一阶导

			matDxytDDH[0] = pKR[0]*tmp2[0] + pKR[1]*tmp2[2] + pKR[2]*tmp2[4];//x关于azimuth的一阶导
			matDxytDDH[1] = pKR[0]*tmp2[1] + pKR[1]*tmp2[3] + pKR[2]*tmp2[5];//x关于elevation的一阶导
			matDxytDDH[2] = pKR[3]*tmp2[0] + pKR[4]*tmp2[2] + pKR[5]*tmp2[4];//y关于azimuth的一阶导
			matDxytDDH[3] = pKR[3]*tmp2[1] + pKR[4]*tmp2[3] + pKR[5]*tmp2[5];//y关于elevation的一阶导
			matDxytDDH[4] = pKR[6]*tmp2[0] + pKR[7]*tmp2[2] + pKR[8]*tmp2[4];//t关于azimuth的一阶导
			matDxytDDH[5] = pKR[6]*tmp2[1] + pKR[7]*tmp2[3] + pKR[8]*tmp2[5];//t关于elevation的一阶导

			matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//u关于azimuth的一阶导
			matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];//u关于elevation的一阶导
			matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//v关于azimuth的一阶导
			matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];//v关于elevation的一阶导
			
			//parallax angle
			matDXkDpa[0] = dDisTi2k*cos(dW+ppt[2])*matXj[0] - cos(ppt[2])*matPosei2k[0];//matXk[0]关于parallax的一阶导
			matDXkDpa[1] = dDisTi2k*cos(dW+ppt[2])*matXj[1] - cos(ppt[2])*matPosei2k[1];//matXk[1]关于parallax的一阶导
			matDXkDpa[2] = dDisTi2k*cos(dW+ppt[2])*matXj[2] - cos(ppt[2])*matPosei2k[2];//matXk[2]关于parallax的一阶导

			tmp1[0] = pKR[0]*matDuvDxyt[0] + pKR[3]*matDuvDxyt[1] + pKR[6]*matDuvDxyt[2];
			tmp1[1] = pKR[1]*matDuvDxyt[0] + pKR[4]*matDuvDxyt[1] + pKR[7]*matDuvDxyt[2];
			tmp1[2] = pKR[2]*matDuvDxyt[0] + pKR[5]*matDuvDxyt[1] + pKR[8]*matDuvDxyt[2];
			tmp1[3] = pKR[0]*matDuvDxyt[3] + pKR[3]*matDuvDxyt[4] + pKR[6]*matDuvDxyt[5];
			tmp1[4] = pKR[1]*matDuvDxyt[3] + pKR[4]*matDuvDxyt[4] + pKR[7]*matDuvDxyt[5];
			tmp1[5] = pKR[2]*matDuvDxyt[3] + pKR[5]*matDuvDxyt[4] + pKR[8]*matDuvDxyt[5];

			matDuvDpa[0] = tmp1[0]*matDXkDpa[0] + tmp1[1]*matDXkDpa[1] + tmp1[2]*matDXkDpa[2];//u关于parallax的一阶导
			matDuvDpa[1] = tmp1[3]*matDXkDpa[0] + tmp1[4]*matDXkDpa[1] + tmp1[5]*matDXkDpa[2];//v关于parallax的一阶导

			//uv关于azimuth、elevation、parallax的一阶导
			pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = matDuvDpa[0];
			pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = matDuvDpa[1];

//posei = pa + nM * 6 + 3;//主锚点Xc
//posek = pa + nN * 6 + 3;//副锚点Xc
//matPosei2k[0] = posek[0] - posei[0];		matPosei2k[1] = posek[1] - posei[1];		matPosei2k[2] = posek[2] - posei[2];//主锚点到副锚点的向量
//dDisTi2k = sqrt(matPosei2k[0] * matPosei2k[0] + matPosei2k[1] * matPosei2k[1] + matPosei2k[2] * matPosei2k[2]);//主锚点到副锚点的向量的模长
//dDot = matXj[0] * matPosei2k[0] + matXj[1] * matPosei2k[1] + matXj[2] * matPosei2k[2];
//dArcCosIn = dDot / dDisTi2k;
//dW = acos(dArcCosIn);
//matXk[0] = dDisTi2k * sin(dW + ppt[2]) * matXj[0] - sin(ppt[2]) * matPosei2k[0];
//matXk[1] = dDisTi2k * sin(dW + ppt[2]) * matXj[1] - sin(ppt[2]) * matPosei2k[1];
//matXk[2] = dDisTi2k * sin(dW + ppt[2]) * matXj[2] - sin(ppt[2]) * matPosei2k[2];
// 
			//Ti
			matDdisDTi[0] = -matPosei2k[0]/dDisTi2k;//dDisTi2k对主锚点坐标Ti的一阶导
			matDdisDTi[1] = -matPosei2k[1]/dDisTi2k;
			matDdisDTi[2] = -matPosei2k[2]/dDisTi2k;


			matDDotDTi[0] = -sin(ppt[0])*cos(ppt[1]);//dDot对主锚点坐标Ti的一阶导
			matDDotDTi[1] = -sin(ppt[1]);
			matDDotDTi[2] = -cos(ppt[0])*cos(ppt[1]);

			
			matDArcCosInDTi[0] = (dDisTi2k*matDDotDTi[0]-dDot*matDdisDTi[0])/(dDisTi2k*dDisTi2k);//dArcCosIn对主锚点坐标Ti的一阶导
			matDArcCosInDTi[1] = (dDisTi2k*matDDotDTi[1]-dDot*matDdisDTi[1])/(dDisTi2k*dDisTi2k);
			matDArcCosInDTi[2] = (dDisTi2k*matDDotDTi[2]-dDot*matDdisDTi[2])/(dDisTi2k*dDisTi2k);

			matDsinwWDTi[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[0];//sin(dW + ppt[2])对主锚点坐标Ti的一阶导
			matDsinwWDTi[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[1];
			matDsinwWDTi[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[2];

			matDTi2kDTi[0] = matDTi2kDTi[4] = matDTi2kDTi[8] = -1;//matPosei2k对主锚点坐标Ti的一阶导
			matDTi2kDTi[1] = matDTi2kDTi[2] = matDTi2kDTi[3] = 0;
			matDTi2kDTi[5] = matDTi2kDTi[6] = matDTi2kDTi[7] = 0;

			matDXkDTi[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[0] + dDisTi2k*matXj[0]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[0];//matXk对主锚点坐标Ti的一阶导
			matDXkDTi[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[1] + dDisTi2k*matXj[0]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[1];
			matDXkDTi[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[2] + dDisTi2k*matXj[0]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[2];
			matDXkDTi[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[0] + dDisTi2k*matXj[1]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[3];
			matDXkDTi[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[1] + dDisTi2k*matXj[1]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[4];
			matDXkDTi[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[2] + dDisTi2k*matXj[1]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[5];
			matDXkDTi[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[0] + dDisTi2k*matXj[2]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[6];
			matDXkDTi[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[1] + dDisTi2k*matXj[2]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[7];
			matDXkDTi[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[2] + dDisTi2k*matXj[2]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[8];

			matDuvDTi[0] = tmp1[0]*matDXkDTi[0] + tmp1[1]*matDXkDTi[3] + tmp1[2]*matDXkDTi[6];//uv对主锚点坐标Ti的一阶导
			matDuvDTi[1] = tmp1[0]*matDXkDTi[1] + tmp1[1]*matDXkDTi[4] + tmp1[2]*matDXkDTi[7];
			matDuvDTi[2] = tmp1[0]*matDXkDTi[2] + tmp1[1]*matDXkDTi[5] + tmp1[2]*matDXkDTi[8];
			matDuvDTi[3] = tmp1[3]*matDXkDTi[0] + tmp1[4]*matDXkDTi[3] + tmp1[5]*matDXkDTi[6];
			matDuvDTi[4] = tmp1[3]*matDXkDTi[1] + tmp1[4]*matDXkDTi[4] + tmp1[5]*matDXkDTi[7];
			matDuvDTi[5] = tmp1[3]*matDXkDTi[2] + tmp1[4]*matDXkDTi[5] + tmp1[5]*matDXkDTi[8];

			pAM[0] = matDuvDTi[0];		pAM[1] = matDuvDTi[1];		pAM[2] = matDuvDTi[2];
			pAM[3] = matDuvDTi[3];		pAM[4] = matDuvDTi[4];		pAM[5] = matDuvDTi[5];				 
			
			//Tk
			matDdisDTk[0] = matPosei2k[0]/dDisTi2k;
			matDdisDTk[1] = matPosei2k[1]/dDisTi2k;
			matDdisDTk[2] = matPosei2k[2]/dDisTi2k;

			matDDotDTk[0] = sin(ppt[0])*cos(ppt[1]);
			matDDotDTk[1] = sin(ppt[1]);
			matDDotDTk[2] = cos(ppt[0])*cos(ppt[1]);

			matDArcCosInDTk[0] = (dDisTi2k*matDDotDTk[0] - dDot*matDdisDTk[0])/(dDisTi2k*dDisTi2k);
			matDArcCosInDTk[1] = (dDisTi2k*matDDotDTk[1] - dDot*matDdisDTk[1])/(dDisTi2k*dDisTi2k);
			matDArcCosInDTk[2] = (dDisTi2k*matDDotDTk[2] - dDot*matDdisDTk[2])/(dDisTi2k*dDisTi2k);

			matDsinwWDTk[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[0];
			matDsinwWDTk[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[1];
			matDsinwWDTk[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[2];

			matDTi2kDTk[0] = matDTi2kDTk[4] = matDTi2kDTk[8] = 1;
			matDTi2kDTk[1] = matDTi2kDTk[2] = matDTi2kDTk[3] = 0;
			matDTi2kDTk[5] = matDTi2kDTk[6] = matDTi2kDTk[7] = 0;


			matDXkDTk[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[0] + dDisTi2k*matXj[0]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[0];
			matDXkDTk[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[1] + dDisTi2k*matXj[0]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[1];
			matDXkDTk[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[2] + dDisTi2k*matXj[0]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[2];
			matDXkDTk[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[0] + dDisTi2k*matXj[1]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[3];
			matDXkDTk[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[1] + dDisTi2k*matXj[1]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[4];
			matDXkDTk[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[2] + dDisTi2k*matXj[1]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[5];
			matDXkDTk[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[0] + dDisTi2k*matXj[2]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[6];
			matDXkDTk[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[1] + dDisTi2k*matXj[2]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[7];
			matDXkDTk[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[2] + dDisTi2k*matXj[2]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[8];

			matDuvDTk[0] = tmp1[0]*matDXkDTk[0] + tmp1[1]*matDXkDTk[3] + tmp1[2]*matDXkDTk[6];//uv对副锚点坐标Tk的一阶导
			matDuvDTk[1] = tmp1[0]*matDXkDTk[1] + tmp1[1]*matDXkDTk[4] + tmp1[2]*matDXkDTk[7];
			matDuvDTk[2] = tmp1[0]*matDXkDTk[2] + tmp1[1]*matDXkDTk[5] + tmp1[2]*matDXkDTk[8];
			matDuvDTk[3] = tmp1[3]*matDXkDTk[0] + tmp1[4]*matDXkDTk[3] + tmp1[5]*matDXkDTk[6];
			matDuvDTk[4] = tmp1[3]*matDXkDTk[1] + tmp1[4]*matDXkDTk[4] + tmp1[5]*matDXkDTk[7];
			matDuvDTk[5] = tmp1[3]*matDXkDTk[2] + tmp1[4]*matDXkDTk[5] + tmp1[5]*matDXkDTk[8];
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv对副锚点参数的一阶导
			pPA[3] = matDuvDTk[0];			pPA[4] = matDuvDTk[1];			pPA[5] = matDuvDTk[2];
			pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
			pPA[9] = matDuvDTk[3];			pPA[10] = matDuvDTk[4];			pPA[11] = matDuvDTk[5];
		}
		else
		{
			ptAngle = pa + nP*6;
			matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
			matXj[1] = sin( ppt[1] );
			matXj[2] = cos( ppt[0] ) * cos( ppt[1] );

			pKR = KR + nP*9;
			pKdA= KdA+ nP*9;
			pKdB= KdB+ nP*9;
			pKdG= KdG+ nP*9;
			
			posei = pa + nM*6 + 3;//主锚点
			posek = pa + nN*6 + 3;//副锚点
			posel = pa + nP*6  + 3;//其它锚点

			matPosei2k[0] = posek[0]-posei[0];		matPosei2k[1] = posek[1]-posei[1];		matPosei2k[2] = posek[2]-posei[2];//主锚点到副锚点
			matPosei2l[0] = posel[0]-posei[0];		matPosei2l[1] = posel[1]-posei[1];		matPosei2l[2] = posel[2]-posei[2];//主锚点到其它锚点

			dDisTi2k = sqrt( matPosei2k[0]*matPosei2k[0] + matPosei2k[1]*matPosei2k[1] + matPosei2k[2]*matPosei2k[2] );
			dDot     = matXj[0]*matPosei2k[0] + matXj[1]*matPosei2k[1] + matXj[2]*matPosei2k[2];
			dArcCosIn= dDot / dDisTi2k;
			dW       = acos( dArcCosIn );

			matXl[0] = dDisTi2k * sin( dW + ppt[2] ) * matXj[0] - sin(ppt[2])*matPosei2l[0];
			matXl[1] = dDisTi2k * sin( dW + ppt[2] ) * matXj[1] - sin(ppt[2])*matPosei2l[1];
			matXl[2] = dDisTi2k * sin( dW + ppt[2] ) * matXj[2] - sin(ppt[2])*matPosei2l[2];

			matxyt[0] = pKR[0]*matXl[0] + pKR[1]*matXl[1] + pKR[2]*matXl[2];
			matxyt[1] = pKR[3]*matXl[0] + pKR[4]*matXl[1] + pKR[5]*matXl[2];
			matxyt[2] = pKR[6]*matXl[0] + pKR[7]*matXl[1] + pKR[8]*matXl[2];

			matDuvDxyt[0] = 1/matxyt[2];	
			matDuvDxyt[1] = 0;
			matDuvDxyt[2] = -matxyt[0]/(matxyt[2]*matxyt[2]);
			matDuvDxyt[3] = 0;	
			matDuvDxyt[4] = 1/matxyt[2];		
			matDuvDxyt[5] = -matxyt[1]/(matxyt[2]*matxyt[2]);	
		
			//camera angle
			matDxytDRG[0] =  pKdG[0]*matXl[0] + pKdG[1]*matXl[1] + pKdG[2]*matXl[2];
			matDxytDRG[1] =  pKdG[3]*matXl[0] + pKdG[4]*matXl[1] + pKdG[5]*matXl[2];
			matDxytDRG[2] =  pKdG[6]*matXl[0] + pKdG[7]*matXl[1] + pKdG[8]*matXl[2];

			matDxytDRB[0] =  pKdB[0]*matXl[0] + pKdB[1]*matXl[1] + pKdB[2]*matXl[2];
			matDxytDRB[1] =  pKdB[3]*matXl[0] + pKdB[4]*matXl[1] + pKdB[5]*matXl[2];
			matDxytDRB[2] =  pKdB[6]*matXl[0] + pKdB[7]*matXl[1] + pKdB[8]*matXl[2];

			matDxytDRA[0] =  pKdA[0]*matXl[0] + pKdA[1]*matXl[1] + pKdA[2]*matXl[2];
			matDxytDRA[1] =  pKdA[3]*matXl[0] + pKdA[4]*matXl[1] + pKdA[5]*matXl[2];
			matDxytDRA[2] =  pKdA[6]*matXl[0] + pKdA[7]*matXl[1] + pKdA[8]*matXl[2];			

			matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv对其它锚点kappa参数的一阶导
			matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];

			matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv对其它锚点phi参数的一阶导
			matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];

			matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv对其它锚点omega参数的一阶导
			matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];
		
			//azimuth and elevation angle
			matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
			matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
			matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);

			matDDotDDH[0] = cos(ppt[0])*cos(ppt[1])*matPosei2k[0] - sin(ppt[0])*cos(ppt[1])*matPosei2k[2];
			matDDotDDH[1] = -sin(ppt[0])*sin(ppt[1])*matPosei2k[0] + cos(ppt[1])*matPosei2k[1]-cos(ppt[0])*sin(ppt[1])*matPosei2k[2];

			matDArcCosInDDH[0] = matDDotDDH[0]/dDisTi2k;
			matDArcCosInDDH[1] = matDDotDDH[1]/dDisTi2k;

			dDWDArcCosIn = -1/sqrt(1-dArcCosIn*dArcCosIn);

			matDsinwWDDH[0] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[0];
			matDsinwWDDH[1] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[1];

			tmp2[0] =  dDisTi2k*matXj[0]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[0];
			tmp2[1] =  dDisTi2k*matXj[0]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[1];
			tmp2[2] =  dDisTi2k*matXj[1]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[2];
			tmp2[3] =  dDisTi2k*matXj[1]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[3];
			tmp2[4] =  dDisTi2k*matXj[2]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[4];
			tmp2[5] =  dDisTi2k*matXj[2]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[5];

			matDxytDDH[0] = pKR[0]*tmp2[0] + pKR[1]*tmp2[2] + pKR[2]*tmp2[4];
			matDxytDDH[1] = pKR[0]*tmp2[1] + pKR[1]*tmp2[3] + pKR[2]*tmp2[5];
			matDxytDDH[2] = pKR[3]*tmp2[0] + pKR[4]*tmp2[2] + pKR[5]*tmp2[4];
			matDxytDDH[3] = pKR[3]*tmp2[1] + pKR[4]*tmp2[3] + pKR[5]*tmp2[5];
			matDxytDDH[4] = pKR[6]*tmp2[0] + pKR[7]*tmp2[2] + pKR[8]*tmp2[4];
			matDxytDDH[5] = pKR[6]*tmp2[1] + pKR[7]*tmp2[3] + pKR[8]*tmp2[5];

			matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//u对azimuth的一阶导
			matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];//u对elevation的一阶导
			matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//v对azimuth的一阶导
			matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];//v对elevation的一阶导
		
			//parallax angle
			matDXkDpa[0] = dDisTi2k*cos(dW+ppt[2])*matXj[0] - cos(ppt[2])*matPosei2l[0];
			matDXkDpa[1] = dDisTi2k*cos(dW+ppt[2])*matXj[1] - cos(ppt[2])*matPosei2l[1];
			matDXkDpa[2] = dDisTi2k*cos(dW+ppt[2])*matXj[2] - cos(ppt[2])*matPosei2l[2];	

			tmp1[0] = pKR[0]*matDuvDxyt[0] + pKR[3]*matDuvDxyt[1] + pKR[6]*matDuvDxyt[2];
			tmp1[1] = pKR[1]*matDuvDxyt[0] + pKR[4]*matDuvDxyt[1] + pKR[7]*matDuvDxyt[2];
			tmp1[2] = pKR[2]*matDuvDxyt[0] + pKR[5]*matDuvDxyt[1] + pKR[8]*matDuvDxyt[2];
			tmp1[3] = pKR[0]*matDuvDxyt[3] + pKR[3]*matDuvDxyt[4] + pKR[6]*matDuvDxyt[5];
			tmp1[4] = pKR[1]*matDuvDxyt[3] + pKR[4]*matDuvDxyt[4] + pKR[7]*matDuvDxyt[5];
			tmp1[5] = pKR[2]*matDuvDxyt[3] + pKR[5]*matDuvDxyt[4] + pKR[8]*matDuvDxyt[5];

			matDuvDpa[0] = tmp1[0]*matDXkDpa[0] + tmp1[1]*matDXkDpa[1] + tmp1[2]*matDXkDpa[2];//u对parallax的一阶导
			matDuvDpa[1] = tmp1[3]*matDXkDpa[0] + tmp1[4]*matDXkDpa[1] + tmp1[5]*matDXkDpa[2];//v对parallax的一阶导

			pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = matDuvDpa[0];
			pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = matDuvDpa[1];
		
			//Ti
			matDdisDTi[0] = -matPosei2k[0]/dDisTi2k;
			matDdisDTi[1] = -matPosei2k[1]/dDisTi2k;
			matDdisDTi[2] = -matPosei2k[2]/dDisTi2k;

			matDDotDTi[0] = -sin(ppt[0])*cos(ppt[1]);
			matDDotDTi[1] = -sin(ppt[1]);
			matDDotDTi[2] = -cos(ppt[0])*cos(ppt[1]);

			matDArcCosInDTi[0] = (dDisTi2k*matDDotDTi[0]-dDot*matDdisDTi[0])/(dDisTi2k*dDisTi2k);
			matDArcCosInDTi[1] = (dDisTi2k*matDDotDTi[1]-dDot*matDdisDTi[1])/(dDisTi2k*dDisTi2k);
			matDArcCosInDTi[2] = (dDisTi2k*matDDotDTi[2]-dDot*matDdisDTi[2])/(dDisTi2k*dDisTi2k);

			matDsinwWDTi[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[0];
			matDsinwWDTi[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[1];
			matDsinwWDTi[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[2];

			matDTi2kDTi[0] = matDTi2kDTi[4] = matDTi2kDTi[8] = -1;
			matDTi2kDTi[1] = matDTi2kDTi[2] = matDTi2kDTi[3] = 0;
			matDTi2kDTi[5] = matDTi2kDTi[6] = matDTi2kDTi[7] = 0;

			matDXkDTi[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[0] + dDisTi2k*matXj[0]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[0];
			matDXkDTi[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[1] + dDisTi2k*matXj[0]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[1];
			matDXkDTi[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[2] + dDisTi2k*matXj[0]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[2];
			matDXkDTi[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[0] + dDisTi2k*matXj[1]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[3];
			matDXkDTi[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[1] + dDisTi2k*matXj[1]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[4];
			matDXkDTi[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[2] + dDisTi2k*matXj[1]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[5];
			matDXkDTi[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[0] + dDisTi2k*matXj[2]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[6];
			matDXkDTi[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[1] + dDisTi2k*matXj[2]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[7];
			matDXkDTi[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[2] + dDisTi2k*matXj[2]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[8];

			matDuvDTi[0] = tmp1[0]*matDXkDTi[0] + tmp1[1]*matDXkDTi[3] + tmp1[2]*matDXkDTi[6];
			matDuvDTi[1] = tmp1[0]*matDXkDTi[1] + tmp1[1]*matDXkDTi[4] + tmp1[2]*matDXkDTi[7];
			matDuvDTi[2] = tmp1[0]*matDXkDTi[2] + tmp1[1]*matDXkDTi[5] + tmp1[2]*matDXkDTi[8];
			matDuvDTi[3] = tmp1[3]*matDXkDTi[0] + tmp1[4]*matDXkDTi[3] + tmp1[5]*matDXkDTi[6];
			matDuvDTi[4] = tmp1[3]*matDXkDTi[1] + tmp1[4]*matDXkDTi[4] + tmp1[5]*matDXkDTi[7];
			matDuvDTi[5] = tmp1[3]*matDXkDTi[2] + tmp1[4]*matDXkDTi[5] + tmp1[5]*matDXkDTi[8];

			pAM[0] = matDuvDTi[0];		pAM[1] = matDuvDTi[1];		pAM[2] = matDuvDTi[2];//uv对主锚点平移向量的一阶导
			pAM[3] = matDuvDTi[3];		pAM[4] = matDuvDTi[4];		pAM[5] = matDuvDTi[5];	
		
			//Tk
			matDdisDTk[0] = matPosei2k[0]/dDisTi2k;
			matDdisDTk[1] = matPosei2k[1]/dDisTi2k;
			matDdisDTk[2] = matPosei2k[2]/dDisTi2k;

			matDDotDTk[0] = sin(ppt[0])*cos(ppt[1]);
			matDDotDTk[1] = sin(ppt[1]);
			matDDotDTk[2] = cos(ppt[0])*cos(ppt[1]);

			matDArcCosInDTk[0] = (dDisTi2k*matDDotDTk[0] - dDot*matDdisDTk[0])/(dDisTi2k*dDisTi2k);
			matDArcCosInDTk[1] = (dDisTi2k*matDDotDTk[1] - dDot*matDdisDTk[1])/(dDisTi2k*dDisTi2k);
			matDArcCosInDTk[2] = (dDisTi2k*matDDotDTk[2] - dDot*matDdisDTk[2])/(dDisTi2k*dDisTi2k);

			matDsinwWDTk[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[0];
			matDsinwWDTk[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[1];
			matDsinwWDTk[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[2];

			matDTi2kDTk[0] = matDTi2kDTk[4] = matDTi2kDTk[8] = 1;
			matDTi2kDTk[1] = matDTi2kDTk[2] = matDTi2kDTk[3] = 0;
			matDTi2kDTk[5] = matDTi2kDTk[6] = matDTi2kDTk[7] = 0;

			matDXkDTk[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[0] + dDisTi2k*matXj[0]*matDsinwWDTk[0];
			matDXkDTk[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[1] + dDisTi2k*matXj[0]*matDsinwWDTk[1];
			matDXkDTk[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[2] + dDisTi2k*matXj[0]*matDsinwWDTk[2];
			matDXkDTk[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[0] + dDisTi2k*matXj[1]*matDsinwWDTk[0];
			matDXkDTk[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[1] + dDisTi2k*matXj[1]*matDsinwWDTk[1];
			matDXkDTk[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[2] + dDisTi2k*matXj[1]*matDsinwWDTk[2];
			matDXkDTk[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[0] + dDisTi2k*matXj[2]*matDsinwWDTk[0];
			matDXkDTk[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[1] + dDisTi2k*matXj[2]*matDsinwWDTk[1];
			matDXkDTk[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[2] + dDisTi2k*matXj[2]*matDsinwWDTk[2];

			matDuvDTk[0] = tmp1[0]*matDXkDTk[0] + tmp1[1]*matDXkDTk[3] + tmp1[2]*matDXkDTk[6];
			matDuvDTk[1] = tmp1[0]*matDXkDTk[1] + tmp1[1]*matDXkDTk[4] + tmp1[2]*matDXkDTk[7];
			matDuvDTk[2] = tmp1[0]*matDXkDTk[2] + tmp1[1]*matDXkDTk[5] + tmp1[2]*matDXkDTk[8];
			matDuvDTk[3] = tmp1[3]*matDXkDTk[0] + tmp1[4]*matDXkDTk[3] + tmp1[5]*matDXkDTk[6];
			matDuvDTk[4] = tmp1[3]*matDXkDTk[1] + tmp1[4]*matDXkDTk[4] + tmp1[5]*matDXkDTk[7];
			matDuvDTk[5] = tmp1[3]*matDXkDTk[2] + tmp1[4]*matDXkDTk[5] + tmp1[5]*matDXkDTk[8];
			
			pAA[0] = matDuvDTk[0];		pAA[1] = matDuvDTk[1];		pAA[2] = matDuvDTk[2];//uv对副锚点平移向量的一阶导
			pAA[3] = matDuvDTk[3];		pAA[4] = matDuvDTk[4];		pAA[5] = matDuvDTk[5];
		
			//Tl
			matDXlDTl[0] = matDXlDTl[4] = matDXlDTl[8] = -sin(ppt[2]);
			matDXlDTl[1] = matDXlDTl[2] = matDXlDTl[3] = 0;
			matDXlDTl[5] = matDXlDTl[6] = matDXlDTl[7] = 0;	

			matDuvDTl[0] = tmp1[0]*matDXlDTl[0] + tmp1[1]*matDXlDTl[3] + tmp1[2]*matDXlDTl[6];
			matDuvDTl[1] = tmp1[0]*matDXlDTl[1] + tmp1[1]*matDXlDTl[4] + tmp1[2]*matDXlDTl[7];
			matDuvDTl[2] = tmp1[0]*matDXlDTl[2] + tmp1[1]*matDXlDTl[5] + tmp1[2]*matDXlDTl[8];
			matDuvDTl[3] = tmp1[3]*matDXlDTl[0] + tmp1[4]*matDXlDTl[3] + tmp1[5]*matDXlDTl[6];
			matDuvDTl[4] = tmp1[3]*matDXlDTl[1] + tmp1[4]*matDXlDTl[4] + tmp1[5]*matDXlDTl[7];
			matDuvDTl[5] = tmp1[3]*matDXlDTl[2] + tmp1[4]*matDXlDTl[5] + tmp1[5]*matDXlDTl[8];
		
			pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv对其它锚点参数的一阶导
			pPA[3] = matDuvDTl[0];			pPA[4] = matDuvDTl[1];			pPA[5] = matDuvDTl[2];
			pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
			pPA[9] = matDuvDTl[3];			pPA[10] = matDuvDTl[4];			pPA[11] = matDuvDTl[5];
		}
}

PBA::PBA(void)
{
	m_e1 = m_e2 = m_e3 = 1E-10;
	m_e4 = 0;	
	m_bProvideXYZ = false;
	m_bFocal = false;
	m_szCameraInit = m_szFeatures = m_szCalibration = m_szXYZ = m_sz3Dpts = m_szCamePose = m_szReport = NULL;

	m_nMaxIter = 100;
	m_Tau      = 1E-6;	
	m_bRobustKernel = false;
	m_nRobustType = 1;
	m_bsolverLM = true;
	m_bsolverGN = false;
	m_delt      = 1;
}


PBA::~PBA(void)
{
	if ( m_szCameraInit!= NULL)	free(m_szCameraInit);
	if ( m_szFeatures!= NULL)	free(m_szFeatures);
	if ( m_szCalibration!= NULL)	free(m_szCalibration);
	if ( m_szXYZ!= NULL)		free(m_szXYZ);
	if ( m_szCamePose!= NULL)	free(m_szCamePose);
	if ( m_sz3Dpts!= NULL)		free(m_sz3Dpts);
	if ( m_szReport!= NULL)		free(m_szReport);
}

bool PBA::ba_run(bool bRobust,
	bool bLM,
	int nMaxIter,
	char* szCam,
	char* szFea,
	char* szXYZ,
	char* szCalib,
	char* szReport,
	char* szPose,
	char* sz3D,
	double Tau)
{
	m_szCameraInit = szCam;                                                       //相机外方位元素
	m_szFeatures   = szFea;                                                       //匹配点
	m_szCalibration= szCalib;                                                     //相机内方位元素
	m_szXYZ        = szXYZ;                                                       //
	m_bRobustKernel= bRobust;                                                     //
	m_bsolverGN    = !bLM;                                                        //
	m_nMaxIter     = nMaxIter;                                                    //
	m_szCamePose   = szPose;                                                      //
	m_sz3Dpts      = sz3D;                                                        //
	m_szReport     = szReport;                                                    //
	m_Tau          = Tau;                                                         //

	BAType ba = pba;
	ba_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );     //给定初值

	if (m_bsolverGN)
		ba_motstr_gn( );
	else
		ba_motstr_levmar();
	return true;
}

bool PBA::ba_run( int argc, char** argv )
{
	bool bTrue = ba_parseArgs( argc, argv );
	if (!bTrue)
	{
		fprintf( stderr, "ParallaxBA: Input wrong commands, please check them!\n");
		return false;
	}
	
	ba_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );

	if (m_bsolverGN)		
		ba_motstr_gn( );
	else					
	{
		ba_motstr_levmar();
	}
	return true;
}


bool PBA::ba_initialize( char* szCamera, char* szFeature,  char* szCalib, char* szXYZ )
{
	printf("ParallaxBA: Parallax Angle Bundle Adjustment Version 1.0\n");
	FILE* fp;

	//must input initial initial camera pose file and projection image points file
	fp = fopen( szCamera, "r" );
	if ( fp == NULL )
	{
		fprintf( stderr, "ParallaxBA: Missing initial camera poses file! \n");
		exit(1);
	}
	else
		fclose(fp);

	fp = fopen( szFeature, "r" );
	if ( fp == NULL )
	{	
		fprintf( stderr, "ParallaxBA: Missing feature projection points file! \n");
		exit(1); 
	}
	else
		fclose(fp);

	if ( szCalib != NULL )
	{
		m_bFocal = false;																											
		m_K = (double*)malloc(9*sizeof(double));                           
		ba_readCameraPoseration(szCalib, m_K);                            
	}	

	if ( szXYZ != NULL )
		m_bProvideXYZ = true;	
	zu = 0;
	savePara = 0;
	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	pba_readAndInitialize( szCamera, szFeature, &m_ncams, &m_n3Dpts, &m_n2Dprojs,&m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
				 &m_imgpts, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort );
	
	std::string originalPath(szCamera); 
	size_t pos = originalPath.find_last_of("/\\");
	std::string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	std::string newPath1 = parentPath + "/" + "InitialParallax.txt";
	std::string newPath2 = parentPath + "/" + "Triangulatedxyz.txt";
	std::string newPath3 = parentPath + "/" + "XYZ.txt";

	if(savePara == 1)
		pba_saveInitialParallax(newPath1.c_str(), m_motstruct);
	if (zu == 1)
		ba_saveTriangulatedxyz(newPath2.c_str(), m_motstruct);
	else if (zu == 2)
		pba_saveInitialXYZ(newPath3.c_str(), pba_angle2xyz(m_motstruct));//转换为xyz
	else
		;//Other codes

	printf( "Number of cameras: %d\n", m_ncams );
	printf( "Number of points: %d\n", m_n3Dpts );
	printf( "Number of projections: %d\n", m_n2Dprojs );
	return true;
}
double *PBA::pba_angle2xyz(double* p)
{
	double *n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));
	static int i, j;
	double* pAngle;
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	int cnp = 6, pnp = 3;
	for (i = 0; i < m_n3Dpts; i++)
	{
		pAngle = p + m_ncams * 6 + i * 3;//azimuth angle
		w = *(p + m_ncams * 6 + i * 3 + 2);//parallax angle

		xj[0] = sin(*(pAngle)) * cos(*(pAngle + 1));//x
		xj[1] = sin(*(pAngle + 1));//y
		xj[2] = cos(*(pAngle)) * cos(*(pAngle + 1));//z

		nM = *(m_archor + i * 3 + 1);
		nA = *(m_archor + i * 3 + 2);

		Tik[0] = -*(p + nM * 6 + 3) + *(p + nA * 6 + 3);
		Tik[1] = -*(p + nM * 6 + 4) + *(p + nA * 6 + 4);
		Tik[2] = -*(p + nM * 6 + 5) + *(p + nA * 6 + 5);

		Dik = sqrt(Tik[0] * Tik[0] + Tik[1] * Tik[1] + Tik[2] * Tik[2]);

		w2 = acos((xj[0] * Tik[0] + xj[1] * Tik[1] + xj[2] * Tik[2]) / Dik);

		//特征点相对于主锚点的xyz坐标
		xk[0] = (Dik * sin(w2 + w) * xj[0]) / sin(w);
		xk[1] = (Dik * sin(w2 + w) * xj[1]) / sin(w);
		xk[2] = (Dik * sin(w2 + w) * xj[2]) / sin(w);

		//特征点的全局xyz坐标
		*(n3DptsT + i * 3) = *(p + nM * 6 + 3) + xk[0];
		*(n3DptsT + i * 3 + 1) = *(p + nM * 6 + 4) + xk[1];
		*(n3DptsT + i * 3 + 2) = *(p + nM * 6 + 5) + xk[2];
		if (i ==2)
		{
			printf("%f %f %f\n", *(n3DptsT + i * 3), *(n3DptsT + i * 3 + 1), *(n3DptsT + i * 3 + 2));
		}
	}

	return n3DptsT;
}
void PBA::pba_saveInitialXYZ(const char* sz3Dpt, double* p)
{
	static int i = 0;
	double dx, dy, dz;
	FILE* fp = NULL, * fpc = NULL;
	int nM, nA;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		fp = fopen(sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			dx = *(p + i * 3);
			dy = *(p + i * 3 + 1);
			dz = *(p + i * 3 + 2);
			nM = *(m_archor + i * 3 + 1);
			nA = *(m_archor + i * 3 + 2);
			//fprintf(fp, "%0.5lf %0.5lf %0.5lf %d %d\n", dx, dy, dz,nM, nA);
			//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1,dx, dy, dz);
			fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", dx, dy, dz);
			if (i == 2)
			{
				printf("%0.5lf %0.5lf %0.5lf\n", dx, dy, dz);
				printf("%d\n", m_n3Dpts);
			}
				
		}
		printf("%d\n", i);
		fclose(fp);
	}
	free(p);
}
void PBA::pba_saveInitialParallax(const char* sz3Dpt, double* p)
{
	//double* n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));

	static int i = 0;
	double azimuth, elevation, parallax;
	FILE* fp = NULL, * fpc = NULL;
	double* pAngle;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		fp = fopen(sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			pAngle = p + m_ncams * 6 + i * 3;

			azimuth = *pAngle;
			elevation = *(pAngle + 1);
			parallax = *(pAngle + 2);
			fprintf(fp, "%0.5lf     %0.5lf     %0.5lf\n", azimuth, elevation, parallax);
			//fprintf(fp, "%0.5lf\n", parallax*180/PI);
		}
		fclose(fp);
	}
	free(p);
}

bool PBA::ba_motstr_levmar( )
{
	if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
	{	
		fprintf( stderr, "ParallaxBA: Missing related source file, please input them first\n");	
		return false;
	}
	
	FILE *fpRe = NULL;
	if (m_szReport!= NULL)
	{
		fpRe = fopen(m_szReport, "w ");
		fprintf(fpRe, "%d  poses, %d 3D features, %d projection\n", m_ncams, m_n3Dpts, m_n2Dprojs );
		fprintf(fpRe, "Levenberg-Marquardt is used\n" );
	}

	int i, ii;
	int n = m_n3Dpts, m = m_ncams, cnp = 6, pnp = 3, mnp = 2;
	int nvis, nuis, itno, issolved, nobs, nvars, nMaxS, nu=2, stop=0;

	double *p = m_motstruct, *x = m_imgpts;	//p pointer refers to unknown parameters, x pointer refers to image coordinate	
	double *U, *V, *W, *e, *eab, *E, *S, *dp, *IV; // pointers into U V W E S IV
	double *pa = NULL, *pb = NULL, *ea, *eb, *dpa, *dpb, *hx, *Ex, *rx, *pdp; // pointers into p, jac, eab and dp respectively 	double *Ex, *rx; 
	double initialerror = 0, error = 0;
	int    nIter = 0, nLinear = 0;
	int Usz, Vsz, Wsz, easz, esz, ebsz, Sblsz, Sdim; 
	static double mu, tmp, p_eL2, eab_inf, pdp_eL2, x2, delt2, init_p_eL2 = 0; 
	double p_L2, dp_L2=DBL_MAX, dF, dL;
	double tau=fabs(m_Tau), eps1=fabs(m_e1),eps3_sq=m_e3*m_e3;
	bool init = false, ordering = true;
	double tStart, tEnd, tTimeUse, t0, t1, t11, t2, t3, t4, t5;
	struct sba_crsm Sidxij, Uidxij;		// S mask and U mask
	double tCons = 0, tSolve = 0, tCost = 0, tTimeIter = 0, tTimeCal = 0;
	int nLmIterations = 0;

	Usz=cnp*cnp;												//6*6
	Vsz=pnp * pnp;												//3*3
	Wsz=cnp * pnp;												//6*3
	esz=mnp; easz=cnp; ebsz=pnp;								//2 、 6 、 3
	Sblsz=cnp*cnp;												//6*6
	Sdim=m * cnp;												//m_ncams（相机数目）*6
	nvis = m_n2Dprojs;											//三维点对应的2D投影点的序号
	mu=eab_inf=0.0;												//0.0

	nuis = ba_ConstructSmask( Sidxij, Uidxij );
	nobs=nvis*mnp;
	nvars=m*cnp + n*pnp;
			
	S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));
	W	=	(double *)malloc(nvis*Wsz*sizeof(double));
	U	=	(double *)malloc(nuis*Usz*sizeof(double));
	V	=	(double *)malloc(n*Vsz*sizeof(double));
	IV	=	(double *)malloc(n*Vsz*sizeof(double));
	e	=	(double *)malloc(nobs*sizeof(double));
	eab	=	(double *)malloc(nvars*sizeof(double));
	E	=	(double *)malloc(m*cnp*sizeof(double));		
	dp	=	(double *)malloc(nvars*sizeof(double));
	hx	=	(double *)malloc(nobs*sizeof(double));	
	pdp	=	(double *)malloc(nvars*sizeof(double));	
	pa=p; pb=p+m*cnp; ea=eab; eb=eab+m*cnp;	dpa=dp; dpb=dp+m*cnp;	

	cholmod_start (&m_cS) ;  	
	//m_cS.print_function = NULL;	
	int *Ap  = (int*)malloc((m+1)*sizeof(int));
	int * Aii = (int*)malloc(m_nS*sizeof(int));
	ba_constructAuxCSSLM( Ap, Aii );

	m_cholSparseE = cholmod_zeros( cnp*m, 1, CHOLMOD_REAL, &m_cS);
	Ex = (double*)m_cholSparseE->x;
	nMaxS = (m_nS-m_ncams)*36+m_ncams*21;	//maximum non-zero element in S matrix 

	m_cholSparseS = cholmod_allocate_sparse(m_ncams*6,m_ncams*6,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_cholSparseS->p;		//column pointer
	Si = (int*)m_cholSparseS->i;		//row pointer
	
	//Compute initial error and initial reprojection error 
	tStart = clock();
	pba_cost(p, hx, m_archor ); 

	//for (int i = 0; i < nobs; i++)
	//{
	//	if(x[i]== 2027.86||x[i] == 456.817)
	//		printf("%f  %f\n", x[i], hx[i]);
	//}
	

	p_eL2 = nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */

	initialerror = p_eL2/nvis;
	printf("Initial Error %0.1lf [%0.8lf]\n", p_eL2, initialerror);

	if( m_szReport != NULL )
		fprintf( fpRe, "Initial Error  %lf\n", initialerror );

	init_p_eL2 = p_eL2;

	//Iteration 
	for(itno=0; itno<m_nMaxIter && !stop; ++itno)
	{
		//Setup S matrix, include two step
		memset( U, 0, nuis*Usz*sizeof(double) );
		memset( ea, 0,m*easz*sizeof(double) );
		memset( V, 0, n*Vsz*sizeof(double));
		memset( IV, 0, n*Vsz*sizeof(double));
		memset( eb, 0, n*ebsz*sizeof(double));
		memset( W, 0, nvis*Wsz*sizeof(double));	
						
		//Step one; compute W V U directly, don't save each projection image Jacobian 
		t0 = clock();
		if ( m_bRobustKernel)
			pba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			pba_jacobian(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		t1 = clock();
		if( itno == 0)
			mu = ba_computeInitialmu( U, V, Uidxij, tau, nvars );

		nLmIterations = 0;
		while(1) //determine increment using adaptive damping
		{
			nLmIterations++;
			t11 = clock();

			ba_inverseVLM( V, IV, Uidxij, mu ); //compute inverse matrix of V
			
			//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
			memset( E, 0, m*easz*sizeof(double));
			memset( S, 0, m_nS*Sblsz*sizeof(double) );
			ba_constructSLM( S, E, U, IV, W, ea, eb, Sidxij, mu );
			t2 = clock();
			
			//Solve linear equation
			//set CSS format using S matrix
			ba_constructCSSLM( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init ); 
			for ( ii = 0; ii < cnp*m; ii++  )
				Ex[ii] = E[ii];	
			
			ba_solveCholmodLM( Ap, Aii, init, ordering);
			nLinear++;
			
			init = true;
			rx = (double*)m_cholSparseR->x;

			if (m_cS.status != CHOLMOD_NOT_POSDEF )
			{
				for ( ii = 0; ii < cnp*m; ii++ )
					dpa[ii] = rx[ii];
				issolved = 1;
			}
			else
				issolved = 1;
			
			t3 = clock();
						
			if(issolved)
			{
				//Solve features
				ba_solveFeatures( W, IV, ea, eb, dpa, dpb );
			
				// Compute ||J^T e||_inf and ||p||^2 
				for(i=0, p_L2=eab_inf=0.0; i<nvars; ++i)
				{
					if(eab_inf < (tmp=fabs(eab[i])))
						eab_inf=tmp;
					
					p_L2+=p[i]*p[i];
				}								
	
				//update
				for(i=0, dp_L2=0.0; i<nvars; ++i)// compute p's new estimate and ||dp||^2 
				{
					pdp[i]=p[i] + (tmp=dp[i]);
					dp_L2+=tmp*tmp;
				}

				//m_e1 = 0;
				if (sqrt(dp_L2)<=m_e1*(sqrt(p_L2)+m_e1))
				{	stop = 1;	break;	}

				ba_updateKR( m_KR, m_KdA, m_KdB, m_KdG, m_K, pdp );
				t4 = clock();
				
				pba_cost(pdp, hx, m_archor );
				pdp_eL2=nrmL2xmy(hx, x, hx, nobs); 	
				error = pdp_eL2/nvis;
				t5 = clock();
								
				if ( m_bRobustKernel )
				{
					pdp_eL2 = 0;
					delt2 = m_delt*m_delt;
					if ( m_nRobustType==1)						//Cauchy Kernel Function
					{
						for ( i = 0; i < m_n2Dprojs; i++ )
						{					
							x2 = hx[i*2]*hx[i*2]+hx[i*2+1]*hx[i*2+1];
							x2 = delt2*log( x2/delt2 + 1 );
							pdp_eL2 += x2;
						}
					}
					else										//Huber Kernel Function
					{
						for ( i = 0; i < m_n2Dprojs; i++ )
						{					
							x2 = hx[i*2]*hx[i*2]+hx[i*2+1]*hx[i*2+1];
							
							if (x2 <= delt2)  // inlier
								x2 = x2;
							else  // outliers
								x2 = 2*sqrt(x2)*m_delt - delt2;
								
							pdp_eL2 += x2;
						}
					}
					error = pdp_eL2/nvis;
				}				

				for(i=0, dL=0.0; i<nvars; ++i)
					dL+=dp[i]*(mu*dp[i]+eab[i]);  //low  
				dF=p_eL2-pdp_eL2;	

				if((dF/dL)>0.0)
				{ 
					if((sqrt(p_eL2)-sqrt(pdp_eL2))<m_e3*sqrt(p_eL2)) 
					{	stop=2;		break;	}

					for(i=0; i<nvars; ++i) 
						p[i]=pdp[i];
					for(i=0; i<nobs; ++i) 
						e[i]=hx[i];		
			

					p_eL2=pdp_eL2;
					if((eab_inf <= eps1))
					{	dp_L2=0.0; 		stop=4;		break;	}

					tmp=(2.0*dF/dL-1.0);
					tmp=1.0-tmp*tmp*tmp;
					mu=mu*( (tmp>=1.0/3.0)? tmp : 1.0/3.0 );
					nu=2;

					tTimeIter = (t5 - t0)*0.001;
					tTimeCal += tTimeIter;

					printf( "Iteration=%d  MSE=%0.8lf   LmIters=%d  Pertime=%0.2lf TotalTime=%0.2lf\n", itno, pdp_eL2/nvis, nLmIterations, tTimeIter, tTimeCal ); 					
					if( m_szReport!= NULL )
						fprintf( fpRe, "Iteration %d  Error  %0.8lf\n", itno, pdp_eL2/nvis );
					nIter++;

					break;
				}
				else
				{
					mu*=nu;
					nu *= 2;
				}
			} 			
		}
		
		if(p_eL2<=eps3_sq) stop=5; 
	}
	if(itno>=m_nMaxIter)
		stop=3;	
	
	//clear memory and print
iterstop:

	sba_crsm_free(&Uidxij);
	sba_crsm_free(&Sidxij);

	cholmod_free_factor(&m_cholFactorS, &m_cS) ;              
	cholmod_l_free_dense(&m_cholSparseE, &m_cS);
	cholmod_l_free_dense(&m_cholSparseR, &m_cS);
	cholmod_finish (&m_cS) ;  
	free(Ap);
	free(Aii);
	cholmod_free_sparse(&m_cholSparseS, &m_cS) ;
	
	tEnd  = clock();
	tTimeUse = tEnd - tStart;

	//save optimal camera pose and feature
	pba_saveXYZ( m_szCamePose, m_sz3Dpts, p, false );

	free(S);	free(W);	free(U);	free(V);	free(IV);
	free(e);	free(eab);	free(E);   	free(dp);	free(hx);	free(pdp);	
	free(m_KR); free(m_KdA);free(m_KdB);free(m_KdG);	

	printf( "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
		nvars, m_n2Dprojs*2, stop, p_eL2/nvis, initialerror, nIter+1, nLinear, tTimeUse / CLOCKS_PER_SEC );
	printf( "ParallaxBA reasons is listed as following:\n" );
	printf( "reason 1: relative change of state vector is small\n" );
	printf( "reason 2: relative change of projection error is small\n" );
	printf( "reason 3: maximum iteration\n");
	printf( "reason 4: maximum value of b is small\n");
	printf( "reason 5: total reprojection error is small\n");

	if( m_szReport!=NULL )
	{
		fprintf( fpRe, "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, error, initialerror, nIter+1, nLinear, tTimeUse*0.001 );

		fprintf( fpRe, "ParallaxBA reasons is listed as following:\n" );
		fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
		fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
		fprintf( fpRe, "reason 3: maximum iteration\n");
		fprintf( fpRe, "reason 4: maximum value of b is small\n");
		fprintf( fpRe, "reason 5: total reprojection error is small\n");

		fclose(fpRe);
	}	
	
	return true;
}

bool PBA::ba_motstr_gn( FixType ft )
{	
	if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
	{	
		printf( "ParallaxBA: Missing related file, please input them first\n");	
		return false;
	}

	FILE *fpRe = NULL;
	if( m_szReport != NULL )
	{
		fpRe = fopen(m_szReport, "w ");
		fprintf(fpRe, "%d  poses, %d 3D features, %d projection\n", m_ncams, m_n3Dpts, m_n2Dprojs );
		fprintf(fpRe, "Gauss-Newton is used\n" );
	}

	int i, ii, nft;
	int n = m_n3Dpts, m = m_ncams, cnp = 6, pnp = 3, mnp = 2;
	int nvis, nuis, itno, nobs, nvars, nMaxS;

	double *p = m_motstruct, *x = m_imgpts;	//p pointer refers to unknown parameters, x pointer refers to image coordinate	
	double *U, *V, *W, *e, *eab, *E, *S, *dp; // pointers into U V W E S IV
	double *pa, *pb, *ea, *eb, *dpa, *dpb, *hx, *Ex, *rx; // pointers into p, jac, eab and dp respectively 	double *Ex, *rx; 
	double initialerror = 0, error = 0, lasterror = 0;
	int Usz, Vsz, Wsz, esz, easz, ebsz, Sblsz, Sdim; 

	double tmp, tmp1, init_p_eL2, p_eL2, eab_inf, pdp_eL2, x2, delt2; 
	double dp_L2=DBL_MAX;
	int nu=2, stop=0;
	eab_inf=0.0;		

	bool init = false;
	bool ordering = true;
	double tStart, tEnd, tTimeUse, tPerTime, t0, t1, t2, t3, t4, t5, tTotalTime = 0;
	struct sba_crsm Sidxij, Uidxij;		// S mask and U mask          
	nuis = ba_ConstructSmask( Sidxij, Uidxij );   

	Usz=cnp*cnp; Vsz=pnp * pnp; Wsz=cnp * pnp; 	
	esz=mnp; easz=cnp; ebsz=pnp; Sblsz=cnp*cnp;	Sdim=m * cnp;	
	nvis = m_n2Dprojs;                                            
	nobs=nvis*mnp;                                               
	nvars=m*cnp + n*pnp;                                         

	S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));         
	W	=	(double *)malloc(nvis*Wsz*sizeof(double));            
	U	=	(double *)malloc(nuis*Usz*sizeof(double));            
	V	=	(double *)malloc(n*Vsz*sizeof(double));               
	e	=	(double *)malloc(nobs*sizeof(double));                 
	eab	=	(double *)malloc(nvars*sizeof(double));              
	E	=	(double *)malloc(m*cnp*sizeof(double));	            
	dp	=	(double *)malloc(nvars*sizeof(double));              
	hx	=	(double *)malloc(nobs*sizeof(double));	              
	pa=p; 
	pb=p+m*cnp; 

	ea=eab;
	eb=eab+m*cnp;

	dpa=dp; 
	dpb=dp+m*cnp;
	
	//Select fix axis for Gauss-Newton
	if ( ft == BA_FixDefault )
	{
		tmp = abs(*(m_motstruct+9)-*(m_motstruct+3));
		ft = BA_FixX;									
		tmp1 = abs(*(m_motstruct+10)-*(m_motstruct+4));
		if ( tmp1 > tmp )
		{ft= BA_FixY; tmp = tmp1;	}
		tmp1 = abs(*(m_motstruct+11)-*(m_motstruct+5));
		if ( tmp1 > tmp )
		{ft= BA_FixZ; }
	}
	nft = ft;
	
	cholmod_start (&m_cS) ;  
	int *Ap  = (int*)malloc((m)*sizeof(int));
	int * Aii = (int*)malloc(m_nS*sizeof(int));
	ba_constructAuxCSSGN( Ap, Aii );

	m_cholSparseE = cholmod_zeros( cnp*m-7, 1, CHOLMOD_REAL, &m_cS);
	Ex = (double*)m_cholSparseE->x;
	nMaxS = (m_nS-m_ncams)*36+m_ncams*21;	//maximum non-zero element in S matrix 
	m_cholSparseS = cholmod_allocate_sparse(m_ncams*6-7,m_ncams*6-7,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_cholSparseS->p;		//column pointer
	Si = (int*)m_cholSparseS->i;		//row pointer
	
	//Compute initial error and initial reprojection error 
	tStart = clock();
	//dpa[0] = dpa[1] = dpa[2] = dpa[3] = dpa[4] = dpa[5] = 0;//未用到
	//dpa[9+nft] = 0;
	pba_cost(p, hx, m_archor );
	p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */

	initialerror = p_eL2/nvis;
	lasterror = initialerror;//Add by zuo
	printf("Initial Error %0.1lf [%0.8lf]\n", p_eL2, initialerror);

	if( m_szReport!= NULL )
		fprintf( fpRe, "Initial Error  %0.8lf\n", initialerror );
	init_p_eL2=p_eL2;	
	
	//Iteration 
	for(itno=0; itno<m_nMaxIter && !stop; ++itno)
	{
		//Setup S matrix, include two step
		memset( U, 0, nuis*Usz*sizeof(double) );
		memset( ea, 0,m*easz*sizeof(double) );
		memset( V, 0, n*Vsz*sizeof(double));
		memset( eb, 0, n*ebsz*sizeof(double));
		memset( W, 0, nvis*Wsz*sizeof(double));

		//Step one; compute W V U directly, don't save each projection image Jacobian 
		t0 = clock();
		if ( m_bRobustKernel)
			pba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			pba_jacobian(p, m_archor, &Uidxij, e, U, ea, V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature);

		t1 = clock();

		ba_inverseVGN( U, V, Uidxij ); //compute inverse matrix of V

		//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
		memset( E, 0, m*easz*sizeof(double));
		memset( S, 0, m_nS*Sblsz*sizeof(double) );
		ba_constructSGN( S, E, U, V, W, ea, eb, Sidxij );//S为方程左侧，E为方程右侧
		t2 = clock();
	
		//Solve equation
		ba_constructCSSGN( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init, nft ); //set CSS format using S matrix//break
		for ( ii = 0; ii < 3+nft; ii++ )
			Ex[ii] = E[6+ii];//第2帧的姿态角

		for ( ii = 0; ii < 6*m-(10+nft); ii++  )
			Ex[3+nft+ii] = E[10+nft+ii];//跳过第2帧的Xc，只存储Yc和Zc

		//Sx为左侧，Ex为右侧

		ba_solveCholmodGN( Ap, Aii, init, ordering);//break
		init = true;
		rx = (double*)m_cholSparseR->x;//解算的相机参数

		if (m_cS.status == CHOLMOD_NOT_POSDEF)
		{
			printf ( "Cholesky failure, writing debug.txt (Hessian loadable by Octave)" );
			goto iterstop;
		}
		else
		{
			for ( ii = 0; ii < 6*m-(10+nft); ii++ )
				dpa[6+ii] = rx[ii];
			for ( ii = 0; ii < 6*m-(10+nft); ii++ )
				dpa[ii+10+nft] = rx[3+nft+ii];//存储解算的相机参数
		}					
		t3 = clock();
		
		//Solve features
		ba_solveFeatures( W, V, ea, eb, dpa, dpb );//break
		
		//update
		for(i=0, dp_L2=0.0; i<nvars; ++i)//nvars：未知参数个数
		{
			p[i]=p[i] + (tmp=dp[i]);//dpa指向dp中的外参改正量，dpb指向dp中的三维点改正量，p指向未知参数
			dp_L2+=tmp*tmp;//存储改正量的平方和
		}

		//dp_L2=sqrt(dp_L2);
		if (dp_L2<1E-8*1E-8)//break
		{	stop = 1;	goto iterstop;	}	//参数的变化已经很小						

		ba_updateKR( m_KR, m_KdA, m_KdB, m_KdG, m_K, p );
		t4 = clock();
		
		pba_cost(p, hx, m_archor ); 
		pdp_eL2=nrmL2xmy(e, x, hx, nobs); 				
		error = pdp_eL2/nvis;//break 计算mse
		t5 = clock();
				
		if ( m_bRobustKernel )
		{
			pdp_eL2 = 0;
			delt2 = m_delt*m_delt;
			if ( m_nRobustType==1)						//Cauchy Kernel Function
			{
				for ( i = 0; i < m_n2Dprojs; i++ )
				{					
					x2 = e[i*2]*e[i*2]+e[i*2+1]*e[i*2+1];
					x2 = delt2*log( x2/delt2 + 1 );
					pdp_eL2 += x2;
				}
			}
			else										//Huber Kernel Function
			{
				for ( i = 0; i < m_n2Dprojs; i++ )
				{					
					x2 = e[i*2]*e[i*2]+e[i*2+1]*e[i*2+1];

					if (x2 <= delt2)  // inlier
						x2 = x2;
					else  // outlier
						x2 = 2*sqrt(x2)*m_delt - delt2;

					pdp_eL2 += x2;
				}
			}
			error = pdp_eL2/nvis;
		}
		
		if (abs(error-lasterror) < 1E-9)//change by zuo
		{	stop = 2;	goto iterstop;		}//与上次mse相比，是否显著改善
		
		tPerTime = (t5-t0)*0.001;//0.001将返回的计时值转换成秒
		tTotalTime += tPerTime;//单次迭代的秒

		printf( "Iteration=%d  MSE=%0.8lf Pertime=%0.2lf TotalTime=%0.2lf\n", itno, pdp_eL2/nvis, tPerTime, tTotalTime );
		//printf("setup: %0.2lf  solve: %0.2lf  update: %0.2lf cost: %0.2lf totoal: %0.2lf\n",
		//(t2-t0)/1000.0,(t3-t2)/1000.0, (t4-t3)/1000.0, (t5-t4)/1000.0,(t5-t0)/1000.0);
		//printf("\n");

		if( m_szReport != NULL )
			fprintf( fpRe, "Iteration %d  Error  %0.8lf\n", itno, pdp_eL2/nvis );

		lasterror = error;		//change by zuo
	}
	if(itno>=m_nMaxIter) stop=3;
	
	//clear memory and print
iterstop:
	cholmod_finish (&m_cS) ;  
	sba_crsm_free(&Uidxij);
	sba_crsm_free(&Sidxij);

	cholmod_free_factor(&m_cholFactorS, &m_cS) ;              
	cholmod_l_free_dense(&m_cholSparseE, &m_cS);
	cholmod_l_free_dense(&m_cholSparseR, &m_cS);
	free(Ap);
	free(Aii);
	cholmod_free_sparse(&m_cholSparseS, &m_cS) ;

	tEnd  = clock();
	tTimeUse = tEnd - tStart;
	

	//save optimal camera pose and feature
	pba_saveXYZ( m_szCamePose, m_sz3Dpts, p );

	free(S);	free(W);	free(U);	free(V);	
	free(e);	free(eab);	free(E);   	free(dp);	free(hx);	
	free(m_KR); free(m_KdA);free(m_KdB);free(m_KdG);

	printf( "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
		nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
	printf( "ParallaxBA reasons is listed as following:\n " );
	printf( "reason 1: relative change of state vector is small\n" );
	printf( "reason 2: relative change of projection error is small\n " );
	printf( "reason 3: maximum iteration\n");

	if( m_szReport != NULL )
	{
		fprintf( fpRe, "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
		fprintf( fpRe, "ParallaxBA reasons is listed as following:\n" );
		fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
		fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
		fprintf( fpRe, "reason 3: maximum iteration\n");
		fclose(fpRe);
	}
	
	return true;
}


void PBA::pba_readProjectionAndInitilizeFeature(FILE *fp,
	double *params, double *projs, char *vmask, int ncams, 
	int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort )
{
	int n;
	int nframes;
	int ptno = 0, cur;

	int nproj2D = 0;
	
	int count = 0;

	int frameno;
	int feastart = 0;

	int nP, nP2;
	
	double* ptr1 = projs;//记录每个投影点的xy而不是序号的（2个2个这样记录和存储）

	int i, j; 
	int  sum, cnp = 6;

	int nFlag;
	
	int *ptr2;
	bool bAdjust;	



	m_smask = (char*)malloc(m_ncams*m_ncams*sizeof(char));
	memset( m_smask, 0, m_ncams*m_ncams*sizeof(char) );

	int* archorEx = new int[m_n3Dpts*2];

	bool bM, bN;
	int max_nframes = -1;
	//read all projection point, initialize three feature angle at the same time
	//读取所有投影点，同时初始化三个特征角***重点***
	while(!feof(fp))//一行一行读Match点
	{
		nFlag = 0;
		n = readNInts( fp, &nframes, 1 );  //读取Match-FeaturePoint文件每一行的第一列，即同名点的个数为nframes/3D点的2D投影个数
		if( n!=1 )
			break;//一般都是1，表示成功
		
		archor[ptno*3] = nframes;//从archor[0]开始，记录3D点的2D投影个数nframes
		cur = 0;
		//test
		//if (nframes > 10)
		//{
		//	int gg = 0;
		//}

		//Add By Zuo
		//double* A = (double*)malloc(nframes * 2 * 4 * sizeof(double));
		//Eigen::MatrixXd A(2 * nframes, 4);

		//test
		bM = bN = false;//这是类中定义的，为实例化的对象服务，重点：
		for( i=0, sum = 0; i<nframes; ++i )//一个投影点一个投影点的读取，按总数nframes
		{
			n = readNInts( fp, &frameno, 1 ); //第二次：读取第一个投影点所在-image index

			nphoto[nproj2D] = frameno;//这是一个投影点的id，存在nphoto，nproj2D从0开始的，
			//后一直累计到全部2D投影，而nphoto刚好就是存了全部2D投影
			nfeature[nproj2D] = ptno;//ptno也是从0开始的，也就是nfeature[0]=0?
			//nfeature目前还不清楚
			nproj2D++;//正常+1

			if(frameno>=ncams)//一般不会发生
			{
				fprintf(stderr, "ParallaxBA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles( fp, ptr1, 2 ); //第二次：读取第一个投影点的像点坐标 ptr1[0] ptr1[1]
			ptr1+=2;//一次记录两个所以+2
			//也就是1个id+2个坐标
			if(n!=3)//一般不会发生
			{
				fprintf(stderr, "ParallaxBA:reading image projections wrong!\n");
				return;
			}

			//第三次：执行完if ( bM && !bN )就执行这里
			if ( bM && bN )//这是其他的
			{
				ptr2 = archorSort+ptno*2;
				//bAdjust = pba_initializeOtheArchors_Mindw( //_Mindw
				bAdjust = pba_initializeOtheArchors( //Maxdw
					projs+feastart*2,
					nphoto+feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams*cnp + ptno * 3, 
					ptr2,
					sum,
					i,
					ptno );
				if ( bAdjust )
				{
					archor[ptno*3+1] = *(nphoto+feastart+ptr2[0]);
					archor[ptno*3+2] = *(nphoto+feastart+ptr2[1]);

					archorEx[ptno*2] = ptr2[0];
					archorEx[ptno*2+1] = ptr2[1];
				}
				sum++;
			}

			//第二次：执行完if ( !bM )就执行这里
			if ( bM && !bN )
			{	//bLast是临时变量
				bool bLast = (i == nframes-1);//第一次比如3，i=0，所以bLast也就是false，不是最后一个
				bool bT = pba_initializeAssoArchor( 
					projs+feastart*2,	//记录每个投影点的xy，2个这样记录&存储, feastart一开始是0，间隔为2
					nphoto+feastart,	//记录每个点的序号，feastart一开始是0，间隔为1
					m_motstruct,
					m_K,
					m_motstruct+m_ncams*cnp+ptno*3,
					0,
					1,
					ptno,
					bLast );

				if (bT)
				{
					archorSort[ptno*2+1] = i;//为锚点排序，比如archorSort[0]=2,则第3个投影点所在的像片为主锚点，archorSort[1]=0表示第1个投影点为副锚点
					//archorSort[2]=1表示第2个投影点为其它锚点
					archor[ptno*3+2] = nphoto[count];
					sum++;

					archorEx[ptno*2+1] = i;

					bN = true;
				}
			}


			//第一次：前面的bM = bN = false;的话，就是执行这里的
			if ( !bM )
			{//bLast是临时变量,如果nframes是2的话那么这个就是主锚点了
				bool bLast = (i == nframes-2);
				bool bT = pba_initializeMainArchor( 
					projs+feastart*2,	//记录每个投影点的xy，2个这样记录&存储, feastart一开始是0，间隔为2
					m_motstruct,		//外参六个参数所有 + 三维点三个所有
					m_K,				//内参（都知道）
					m_motstruct+m_ncams*cnp+ptno*3,//m_ncams是多少个相机/像片,cnp从6开始也就是外参的间隔，ptno从0开始，*3为PBA特征参数的间隔
					nphoto[count],		//nphoto存了全部2D投影点的序号
					ptno,				//ptno存了3D点的当前序号
					m_KR );				//KR共线方程变换矩阵
				/*重点注意！：
				* 数组*1为序号存储
				* 数组*2为2D点存储
				* 数组*3为3D点存储
				* 其他的数组为外参，内参等等存储
				* 非数组为的int/float或者double为当前序号
				*/

				///*调试用*/
				//if (bT == true)
				//	printf("%s\n", "true");
				//else
				//	printf("%s\n", "false");
				//printf("%f %f %f\n", (m_motstruct + m_ncams*cnp + ptno * 3)[0], (m_motstruct + m_ncams*cnp + ptno * 3)[1], (m_motstruct + m_ncams*cnp + ptno * 3)[2]);
				///*调试用*/
				archorSort[ptno*2] = i;//第i个投影点即第i个anchor
				archor[ptno*3+1] = nphoto[count];
				sum++;

				archorEx[ptno*2] = i;
				bM = true;//第二次：执行他处
			}	
			count++;//这里和i一样把？每循环一次就—+1，可能只是想多一个变量？		
		}
		
		//set masks for U and S matrix              umask存储U矩阵，m_smask存储S矩阵
		for( i = 0; i < nframes; i++ )
		{
			nP = nphoto[feastart+i];                         //第i个观测的视图号
			int nM_ = archor[ptno * 3 + 1];
			int nA_ = archor[ptno * 3 + 2];
			int tmp3 = archor[ptno * 3 + 3];
			//umask 代表的是 U 矩阵的稀疏结构，只在 主锚点、副锚点、当前观测相机 之间记录 1，因为这些相机参数直接影响该点的观测方程。
			if (nM_<nP)                         
				umask[nM_*(ncams)+nP] = 1;//(主锚点，当前锚点)
			else                                             
				umask[nP*(ncams)+nM_] = 1;//(当前锚点，主锚点)

			
			if (nA_<nP)                         
				umask[nA_*(ncams)+nP] = 1;//（副锚点，当前锚点）
			else                                             
				umask[nP*(ncams)+nA_] = 1;//（当前锚点，副锚点）

			umask[nP*ncams+nP] = 1;//（当前锚点，当前锚点）

			for ( j = i; j < nframes; j++  )
			{
				nP2 = nphoto[feastart+j];//第j个观测的视图号                     

				if ( nP == nP2 )                              //
					m_smask[nP*m_ncams+nP2] = 1;//（第i个观测的视图号，第j个观测的视图号）
				else if ( nP < nP2 )
					m_smask[nP*m_ncams+nP2] = 1;//（第i个观测的视图号，第j个观测的视图号）
				else
				{
					m_smask[nP2*m_ncams + nP] = 1;//（第j个观测的视图号，第i个观测的视图号）未执行过
				}
					
			}
		}					
		feastart += nframes;
		ptno++;//3D特征点的序号，point NO. OR photo NO.
	}
	//count number of non-zero element in S matrix
	m_nS = 0;
	for ( i = 0; i < m_ncams; i++ ) 
	{
		for (j = 0; j < m_ncams; j++)
		{
			if (m_smask[i*m_ncams + j] == 1)
			{
				m_nS++;
			}
		}
	}
}

void PBA::pba_readAndInitialize(char *camsfname, char *ptsfname, int *ncams,
	int *n3Dpts, int *n2Dprojs,
	double **motstruct, double **imgpts,
	int **archor, char **vmask,
	char **umask, int **nphoto,
	int** nfeature, int** archorSort)
{
	FILE *fpc, *fpp, *fpXYZ;
	int i, tmp1, tmp2;
	double ptMain[3], ptA[3];
	double dW1, dW2;	

	//calculate number of cameras, 3D points and projection points
	fpc		=	fopen( camsfname, "r" );//打开外方位元素文件
	*ncams	=	findNcameras( fpc );//读取像片的张数
	m_ncams =	*ncams;

	fpp		=	fopen( ptsfname, "r" );//打开匹配的同名点文件
	readNpointsAndNprojections( fpp, n3Dpts, 3, n2Dprojs, 2 );//读取3D特征点的个数，投影点个数

	*motstruct	=	(double *)malloc( (*ncams*6 + *n3Dpts*3)*sizeof(double) );//用于存储所有的外方位元素 和 所有特征点的3D坐标
	if(	*motstruct==NULL )
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'motstruct' failed \n");
		exit(1);
	}

	*imgpts	=	(double *)malloc(*n2Dprojs*2*sizeof(double));//用于存储所有投影点的2D坐标
	if(	*imgpts==NULL )
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'imgpts' failed\n");
		exit(1);
	}
	//如果要把文件从头读出，需把指针移动到文件头，利用该函数
	rewind(fpc);//每读取一个字符，文件内部位置指针向后移动一个字节，读取完毕，该指针已指向文件末尾，
	rewind(fpp);

	//allocate indicator of U
	*umask = (char*)malloc(*ncams * *ncams );
	memset(*umask, 0, *ncams * *ncams * sizeof(char));//第1个参数一定要是一个已知的，已经被分配内存的地址，第3个参数一定要使用sizeof操作符。
	//对一块已经分配地址的内存进行初始化，并且通常初始化0或者字符'\0'

	//allocate main and associate anchors
	*archor = (int*)malloc(*n3Dpts*3*sizeof(int));//
	memset( *archor, -1, *n3Dpts*3*sizeof(int) ); //将所有主锚点的位置初始化为-1

	*nphoto		= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*nfeature	= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*archorSort = (int*)malloc(*n3Dpts*3*sizeof(int));

	ba_readCameraPose(fpc, *motstruct);
	///*-------------------------------调试用-----------------------------------------*/
	//for (int i = 0;i < *ncams * 6;i++)
	//{
	//	printf("%f ", motstruct[0][i]);
	//	if ((i + 1) % 6 == 0)
	//		printf("\n");
	//}
	///*-------------------------------调试用-----------------------------------------*/
	fclose(fpc);//已关闭
	
	//Update KR
	m_KR  = (double*)malloc(m_ncams*9*sizeof(double));//相机内参矩阵 与 旋转矩阵的乘积
	m_KdA = (double*)malloc(m_ncams*9*sizeof(double));//旋转矩阵M_x (ω) M_y (φ) M_z (κ)的一阶导
	m_KdB = (double*)malloc(m_ncams*9*sizeof(double));//旋转矩阵M_x (ω) M_y (φ) M_z (κ)的一阶导
	m_KdG = (double*)malloc(m_ncams*9*sizeof(double));//旋转矩阵M_x (ω) M_y (φ) M_z (κ)的一阶导
	
	if (zu == 1)
	{
		//构建P矩阵
		m_P = (double*)malloc(m_ncams * 12 * sizeof(double));
		ba_constructP(m_P, m_K, *motstruct);
	}
	else
		ba_updateKR(m_KR, m_KdA, m_KdB, m_KdG, m_K, *motstruct);//第5个参数为相机内参，第6个参数为相机外参


	//test
	//for (int i = 0; i < m_ncams; i++)
	//{
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9], m_KR[i * 9 + 1], m_KR[i * 9 + 2], m_P[i * 12], m_P[i * 12 + 1], m_P[i * 12 + 2]);
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9 + 3], m_KR[i * 9 + 4], m_KR[i * 9 + 5], m_P[i * 12 + 4], m_P[i * 12 + 5], m_P[i * 12 + 6]);
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9 + 6], m_KR[i * 9 + 7], m_KR[i * 9 + 8], m_P[i * 12 + 8], m_P[i * 12 + 9], m_P[i * 12 + 10]);
	//}


	//if XYZ are provided, we can use them as feature initialization.
	if (m_bProvideXYZ)
	{
		fpXYZ = fopen( m_szXYZ, "r");
		m_XYZ = (double*)malloc(m_n3Dpts*3*sizeof(double));

		for( i = 0; i < m_n3Dpts; i++)
			fscanf( fpXYZ, "%lf  %lf  %lf", m_XYZ+i*3, m_XYZ+i*3+1, m_XYZ+i*3+2 );
		fclose(fpXYZ);
	}
	//6个EOP，3个FeatureXYZ
	/*
	* 这是个很关键的函数
	*/
	if(zu==1)
		ba_readProjectionAndTriangulateFeature(fpp, *imgpts, *ncams);
	else
	{
		pba_readProjectionAndInitilizeFeature(fpp,
			*motstruct + *ncams * 6,
			*imgpts,
			*vmask,
			*ncams,
			*archor,
			*umask,
			*nphoto,
			*nfeature,
			*archorSort);
	}

	fclose(fpp);//上面结束后关闭

	
	int nCount = 0;
	double pti2k[3];
	int cur = 0;
	if ( m_bProvideXYZ )
	{
		for ( i = 0; i < m_n3Dpts; i++ )
		{
			int nM = m_archor[i*3+1];
			int nN = m_archor[i*3+2];

			ptMain[0] = *(*motstruct + nM*6 + 3);
			ptMain[1] = *(*motstruct + nM*6 + 4);
			ptMain[2] = *(*motstruct + nM*6 + 5);

			ptA[0] = *(*motstruct + nN*6 + 3);
			ptA[1] = *(*motstruct + nN*6 + 4);
			ptA[2] = *(*motstruct + nN*6 + 5);

			pti2k[0] = ptA[0] - ptMain[0];
			pti2k[1] = ptA[1] - ptMain[1];
			pti2k[2] = ptA[2] - ptMain[2];

			double dispti2k;
			dispti2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1]  + pti2k[2]*pti2k[2] );

			ptMain[0] = m_XYZ[i*3] - ptMain[0];
			ptMain[1] = m_XYZ[i*3+1] - ptMain[1];
			ptMain[2] = m_XYZ[i*3+2] - ptMain[2];	

			ptA[0] = m_XYZ[i*3] - ptA[0];
			ptA[1] = m_XYZ[i*3+1] - ptA[1];
			ptA[2] = m_XYZ[i*3+2] - ptA[2];

			dW1 = ptMain[0]*ptMain[0] + ptMain[1]*ptMain[1] + ptMain[2]*ptMain[2];
			dW2 = ptA[0]*ptA[0] + ptA[1]*ptA[1] + ptA[2]*ptA[2];

			double disDot2;
			disDot2 = ptMain[0]*pti2k[0] + ptMain[1]*pti2k[1] + ptMain[2]*pti2k[2]; 
			double dww = disDot2/(dispti2k * sqrt(dW1));


			double* pKR = m_KR + nM*9;
			double n[2], n2[2], ptXj[3];

			ptXj[0] = ptMain[0];	ptXj[1] = ptMain[1];	ptXj[2] = ptMain[2];
			n[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
				(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			n[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
				(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			pKR = m_KR + nN*9;
			ptXj[0] = ptA[0];	ptXj[1] = ptA[1];	ptXj[2] = ptA[2];
			n2[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
				(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			n2[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
				(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			//printf("%d %d %d\n", m_archorSort[i * 2], m_archorSort[i * 2 + 1], m_archor[i * 3]);
			int id1 = cur + m_archorSort[i*2];
			int id2 = cur + m_archorSort[i*2+1];
			double err1 = (m_imgpts[id1*2]-n[0])*(m_imgpts[id1*2]-n[0])+(m_imgpts[id1*2+1]-n[1])*(m_imgpts[id1*2+1]-n[1]);
			double err2 = (m_imgpts[id2*2]-n2[0])*(m_imgpts[id2*2]-n2[0])+(m_imgpts[id2*2+1]-n2[1])*(m_imgpts[id2*2+1]-n2[1]);

			cur += m_archor[i*3];
			if ( (sqrt(dW1)/sqrt(dW2)>30) || (err1>err2)&&(m_archor[3*i]==2) )
			{
				nCount++;
				m_archor[i*3+1] = nN;
				m_archor[i*3+2] = nM;

				tmp1 = m_archorSort[i*2] ;
				tmp2 = m_archorSort[i*2+1] ;

				m_archorSort[i*2] = tmp2;
				m_archorSort[i*2+1] = tmp1;

				double dDAngle = atan2( ptA[0], ptA[2] );
				double dHAngle = atan2( ptA[1], sqrt(ptA[0]*ptA[0]+ ptA[2]*ptA[2]) );

				(*motstruct)[m_ncams*6+i*3] = dDAngle;
				(*motstruct)[m_ncams*6+i*3+1] = dHAngle;

				double dwwDot = ptMain[0]*ptA[0] + ptMain[1]*ptA[1] + ptMain[2]*ptA[2];				
			}	 
		}
	}

	if(m_bProvideXYZ)
		free(m_XYZ);

	fpXYZ = NULL;
}

void PBA::pba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
{
	int nM, nN;	  
	register int i, j, ii, jj, k;
	int cnp, pnp, mnp;
	double *pa, *pb, *ppt;	
	double *ppUpa, *pea, *peb, *pe, *pV, *pW;
	int pos, pos2, nP, nF, numfea, cur;
	double sum;
	register double tmp = 0;

	double pAM[6], pAA[6], pPA[12], pPB[6];

	cnp = 6; pnp = 3; mnp = 2;
	pa=p; pb=p+m*cnp;

	pos2 = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];//三维点的共视图数量
		cur = 0;
		nF = i;
		//printf("%d %d\n", m_archorSort[i * 2], m_archorSort[i * 2 + 1]);
		
		for ( j = 0; j < numfea; j++ )
		{
			
			nP = nphoto[pos2];//共视图编号
			ppt=pb + nF*pnp;//三维点坐标  
			pe= e + pos2*2;		//error

			nM = archor[nF*3+1];//主锚点
			nN = archor[nF*3+2];//副锚点	
			//pa指向相机参数
			//ppt指向三维点参数
			//nP当前视图编号
			//nM主锚点编号
			//nA副锚点编号
			//pPA指向uv对当前视图参数的一阶导
			//pPB指向uv对三维点的一阶导
			//pAM指向uv（非主锚点视图的uv）对主锚点视图位移参数的一阶导
			//pAA指向uv（非主副锚点视图的uv）对副锚点视图位移参数的一阶导

			pba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );
			
			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);//返回U矩阵中nP行nP列的元素在稀疏结构中的列索引
			//printf("%d  ", pos);
			ppUpa = U + pos*(6*6);//按块填充，会不会存在同一块被重复填充的问题？同一块是不同观测偏导的相加过程，不会被重复填充

			//printf("\n\n\n\n");
			for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//由于是对称矩阵，只填上三角
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii]*pPA[k*6+jj];//U^T*U
				ppUpa[ii*6+jj] += sum;
				//printf("%d ", ppUpa[ii * 6 + jj]);
			}

			//ea
			pea = ea + nP*6;
			for ( ii = 0; ii < 6; ii++ )
			{
				for( jj = 0, sum = 0; jj < 2; jj++ )
					sum += pPA[jj*6+ii] * pe[jj];
				pea[ii] += sum;
			}
			//V
			pV = V + nF*3*3;
			for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii]*pPB[k*3+jj];
				pV[ii*3+jj] += sum;

			}
			//eb
			peb = eb + nF*3;
			for ( ii = 0; ii < 3; ii++ )
			{
				for( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii] * pe[k];
				peb[ii] += sum;
			}
			//W
			pW = W + pos2 * 3*6;
			for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii] * pPB[k*3+jj];
				pW[ii*3+jj] += sum;
			}
			if ( nP == nM )
			{
				cur++;
				pos2++;
				continue;
			}
			else
			if( nP == nN )
			{
				//U		  
				pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
				ppUpa = U + pos*(6*6);
				for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pAM[k*3+ii] * pAM[k*3+jj];
					ppUpa[(ii+3)*6+3+jj] += sum;
				}
				if( nM < nP )
				{
					pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPA[k*6+jj];
						ppUpa[(ii+3)*6+jj] += sum;
					}
				}
				else
				{
					pos = sba_crsm_elmidx(Uidxij, nP, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum +=  pPA[k*6+ii]*pAM[k*3+jj];
						ppUpa[ii*6+jj+3] += sum;
					}
				}

				//ea
				pea = ea + nM*6;
				for ( ii = 0; ii < 3; ii++ )
				{
					for( k = 0, sum = 0; k < 2; k++ )
						sum += pAM[k*3+ii] * pe[k];//当下视图（非主锚点）坐标对主锚点位移参数的偏导 * 当下视图坐标的预测残差
					pea[ii+3] += sum;
				}
				//W
				pos = pos2 + (m_archorSort[i*2]-j);//循环尾部pos2++，所以这个地方-j
				pW = W + pos*6*3;
				for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pAM[k*3+ii] * pPB[k*3+jj];
					pW[(ii+3)*3+jj] += sum;
				}
			}
			else
			{
				//U
				pos = sba_crsm_elmidx(Uidxij, nM, nM );		
				ppUpa = U + pos*(6*6);
				for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
					{
						sum += pAM[k*3+ii] * pAM[k*3+jj];
					}
					
					ppUpa[(ii+3)*6+3+jj] += sum;
				}

				pos = sba_crsm_elmidx(Uidxij, nN, nN );		
				ppUpa = U + pos*(6*6);
				for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
					{
						sum += pAA[k*3+ii] * pAA[k*3+jj];
					}
					ppUpa[(ii+3)*6+3+jj] += sum;
				}
				if ( nM < nN )
				{
					pos = sba_crsm_elmidx(Uidxij, nM, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAM[k*3+ii] * pAA[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum;
					}
				}	
				else
				{
					pos = sba_crsm_elmidx(Uidxij, nN, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAA[k*3+ii] * pAM[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum;
					}
				}
				if ( nM < nP )
				{
					pos = sba_crsm_elmidx(Uidxij, nM, nP );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPA[k*6+jj];
						ppUpa[(ii+3)*6+jj] += sum;
					}
				}
				else
				{
					pos = sba_crsm_elmidx(Uidxij, nP, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum +=  pPA[k*6+ii]*pAM[k*3+jj];
						ppUpa[ii*6+jj+3] += sum;
					}
				}
				
				if ( nN < nP )
				{
					pos = sba_crsm_elmidx(Uidxij, nN, nP );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pPA[k*6+jj];
						ppUpa[(ii+3)*6+jj] += sum;
					}
				}	
				else
				{
					pos = sba_crsm_elmidx(Uidxij, nP, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum +=  pPA[k*6+ii]*pAA[k*3+jj];
						ppUpa[ii*6+jj+3] += sum;
					}
				}
				//ea
				pea = ea + nM*6;
				for ( ii = 0; ii < 3; ii++ )
				{
					for( k = 0, sum = 0; k < 2; k++ )
						sum += pAM[k*3+ii] * pe[k];
					pea[ii+3] += sum;
				}
				pea = ea + nN*6;
				for ( ii = 0; ii < 3; ii++ )
				{
					for( k = 0, sum = 0; k < 2; k++ )
						sum += pAA[k*3+ii] * pe[k];
					pea[ii+3] += sum;
				}  
				//W
				pos = pos2 + (m_archorSort[i*2]-j);
				pW = W + pos*6*3;
				for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pAM[k*3+ii] * pPB[k*3+jj];
					pW[(ii+3)*3+jj] += sum;
				}
				pos = pos2 + (m_archorSort[i*2+1]-j);
				pW = W + pos*6*3;
				for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pAA[k*3+ii] * pPB[k*3+jj];
					pW[(ii+3)*3+jj] += sum;
				}
			}
			pos2++;
			cur++;
		}
	}
}

void PBA::pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
{
	int nM, nN;	  
	register int i, j, ii, jj, k;
	int cnp, pnp, mnp;
	double *pa, *pb, *ppt;	
	double *ppUpa, *pea, *peb, *pe, *pV, *pW;
	int pos, pos2, nP, nF, numfea, cur;
	double sum, x2;
	register double tmp = 0;

	double pAM[6], pAA[6], pPA[12], pPB[6];

	cnp = 6; pnp = 3; mnp = 2;
	pa=p; pb=p+m*cnp;

	pos2 = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		cur = 0;
		nF = i;
		for ( j = 0; j < numfea; j++ )
		{
			nP = nphoto[pos2];
			ppt=pb + nF*pnp;      
			pe= e + pos2*2;		//error

			nM = archor[nF*3+1];
			nN = archor[nF*3+2];

			double delt2 = m_delt*m_delt;
			if( m_nRobustType == 1)
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];
				x2 = 1.0/( x2/delt2+1 );
			}
			else
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];

				if (sqrt(x2) < m_delt) //inlier
					x2 = 1;
				else // outliers 
					x2 = m_delt/sqrt(x2);
			}			

			pba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );

			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);		
			ppUpa = U + pos*(6*6);

			for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//diag
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii]*pPA[k*6+jj];
				//ppUpa[ii*6+jj] += sum;
				ppUpa[ii*6+jj] += sum*x2;
			}

			//ea
			pea = ea + nP*6;
			for ( ii = 0; ii < 6; ii++ )
			{
				for( jj = 0, sum = 0; jj < 2; jj++ )
					sum += pPA[jj*6+ii] * pe[jj];
				//pea[ii] += sum;
				pea[ii] += sum*x2;
			}
			//V
			pV = V + nF*3*3;
			for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii]*pPB[k*3+jj];
				//pV[ii*3+jj] += sum;
				pV[ii*3+jj] += sum*x2;

			}
			//eb
			peb = eb + nF*3;
			for ( ii = 0; ii < 3; ii++ )
			{
				for( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii] * pe[k];
				//peb[ii] += sum;
				peb[ii] += sum*x2;
			}
			//W
			pW = W + pos2 * 3*6;
			for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii] * pPB[k*3+jj];
				//pW[ii*3+jj] += sum;
				pW[ii*3+jj] += sum*x2;
			}

			if ( nP == nM )
			{
				cur++;
				pos2++;
				continue;
			}
			else
				if( nP == nN )
				{
					//U		  
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						//ppUpa[(ii+3)*6+3+jj] += sum;
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							//ppUpa[(ii+3)*6+jj] += sum;
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							//ppUpa[ii*6+jj+3] += sum;
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						//pea[ii+3] += sum;
						pea[ii+3] += sum*x2;
					}
					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						//pW[(ii+3)*3+jj] += sum;
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}
				else
				{
					//U
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						}
						//ppUpa[(ii+3)*6+3+jj] += sum;
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					pos = sba_crsm_elmidx(Uidxij, nN, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAA[k*3+ii] * pAA[k*3+jj];
						}
						//ppUpa[(ii+3)*6+3+jj] += sum;
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if ( nM < nN )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAM[k*3+ii] * pAA[k*3+jj];
							}
							//ppUpa[(ii+3)*6+3+jj] += sum;
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAA[k*3+ii] * pAM[k*3+jj];
							}
							//ppUpa[(ii+3)*6+3+jj] += sum;
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}

					if ( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							//ppUpa[(ii+3)*6+jj] += sum;
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							//ppUpa[ii*6+jj+3] += sum;
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					if ( nN < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pPA[k*6+jj];
							//ppUpa[(ii+3)*6+jj] += sum;
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAA[k*3+jj];
							//ppUpa[ii*6+jj+3] += sum;
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}


					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						//pea[ii+3] += sum;
						pea[ii+3] += sum*x2;
					}

					pea = ea + nN*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pe[k];
						//pea[ii+3] += sum;
						pea[ii+3] += sum*x2;
					}  

					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						//pW[(ii+3)*3+jj] += sum;
						pW[(ii+3)*3+jj] += sum*x2;
					}

					pos = pos2 + (m_archorSort[i*2+1]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pPB[k*3+jj];
						//pW[(ii+3)*3+jj] += sum;
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}

				pos2++;
				cur++;

		}
	}
}
/*void PBA::pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
{
	int nM, nN;	  
	register int i, j, ii, jj, k;
	int cnp, pnp, mnp;
	double *pa, *pb, *ppt;	
	double *ppUpa, *pea, *peb, *pe, *pV, *pW;
	int pos, pos2, nP, nF, numfea, cur;
	static double sum, x2, delt2;
	register double tmp = 0;

	double pAM[6], pAA[6], pPA[12], pPB[6];

	cnp = 6; pnp = 3; mnp = 2;
	pa=p; pb=p+m*cnp;
	delt2 = m_delt*m_delt;

	pos2 = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		cur = 0;
		nF = i;
		for ( j = 0; j < numfea; j++ )
		{
			nP = nphoto[pos2];
			ppt=pb + nF*pnp;      
			pe= e + pos2*2;		//error

			nM = archor[nF*3+1];
			nN = archor[nF*3+2];	

			if( m_nRobustType == 1)
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];
				x2 = 1.0/( x2/delt2+1 );
			}
			else
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];

				if (sqrt(x2) < m_delt) //inlier
					x2 = 1;
				else // outliers 
					x2 = m_delt/sqrt(x2);
			}			

			pba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );

			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);		
			ppUpa = U + pos*(6*6);

			for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//diag
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii]*pPA[k*6+jj];
				ppUpa[ii*6+jj] += sum*x2;
			}

			//ea
			pea = ea + nP*6;
			for ( ii = 0; ii < 6; ii++ )
			{
				for( jj = 0, sum = 0; jj < 2; jj++ )
					sum += pPA[jj*6+ii] * pe[jj];
				pea[ii] += sum*x2;
			}
			//V
			pV = V + nF*3*3;
			for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii]*pPB[k*3+jj];
				pV[ii*3+jj] += sum*x2;

			}
			//eb
			peb = eb + nF*3;
			for ( ii = 0; ii < 3; ii++ )
			{
				for( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii] * pe[k];
				peb[ii] += sum*x2;
			}
			//W
			pW = W + pos2 * 3*6;
			for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii] * pPB[k*3+jj];
				pW[ii*3+jj] += sum*x2;
			}

			if ( nP == nM )
			{
				cur++;
				pos2++;
				continue;
			}
			else
				if( nP == nN )
				{
					//U		  
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						pea[ii+3] += sum*x2;
					}
					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}
				else
				{
					//U
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					pos = sba_crsm_elmidx(Uidxij, nN, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAA[k*3+ii] * pAA[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if ( nM < nN )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAM[k*3+ii] * pAA[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAA[k*3+ii] * pAM[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}

					if ( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					if ( nN < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAA[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}


					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						pea[ii+3] += sum*x2;
					}

					pea = ea + nN*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
						pea[ii+3] += sum*x2;
					}  

					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}

					pos = pos2 + (m_archorSort[i*2+1]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}

				pos2++;
				cur++;
		}
	}
}
*/

void PBA::pba_cost(double* p, double* hx, int* archor)
{
	register int i;
	int cnp, pnp, mnp;
	double* pa, * pb, * ppt, * pmeas;
	int nF, nP, nM, nN;
	int m;

	cnp = 6, pnp = 3, mnp = 2;
	m = m_ncams;
	pa = p; pb = p + m * cnp;

	for (i = 0; i < m_n2Dprojs; i++)
	{
		nF = m_feature[i];
		nP = m_photo[i];

		ppt = pb + nF * pnp;
		pmeas = hx + i * mnp; // set pmeas to point to hx_ij 存储重投影像素点

		nM = archor[nF * 3 + 1];
		nN = archor[nF * 3 + 2];

		pba_reprojectEachPts(m_KR, pa, ppt, nM, nN, nP, pmeas);
		//printf("%f %f\n", pmeas[0], pmeas[1]);
	}
}


/*初始化主锚点*/
bool PBA::pba_initializeMainArchor( 
	double* imgpts,	// 投影点的 2D 图像坐标（x, y）。
	double* camera,	// 相机外参（包括位置和旋转）。
	double* K,		// 相机内参矩阵（焦距、主点位置等）。
	double* feature,// 输出的特征点信息（方位角Azimuth、俯仰角Elevation、3D 初始化）。
	int nP,			// 当前锚点对应的图像编号（相机ID）。
	int FID,		// 特征点的编号。
	double* KR )	// 内参与旋转矩阵的乘积（K * R）。
{
	/*
	bool bT = pba_initializeMainArchor( 
					projs+feastart*2,	//记录每个投影点的xy，2个这样记录&存储, feastart一开始是0，间隔为2
					m_motstruct,		//外参六个参数所有 + 三维点三个所有
					m_K,				//内参（都知道）
					m_motstruct+m_ncams*cnp+ptno*3,//m_ncams是多少个相机/像片,cnp从6开始也就是外参的间隔，ptno从0开始，*3为三维点的间隔
					nphoto[count],		//nphoto存了全部2D投影点的序号
					ptno,				//ptno存了3D点的当前序号
					m_KR );				//KR共线方程变换矩阵
	*/
	//solve  KRX = x
	Vector3d x;//特征点相对于/投影于主锚点的XYZ坐标
	//这里三维点肯定是不知道的。
	if (m_bProvideXYZ)//如果特征点的3D坐标已知
	{
		x(0) = m_XYZ[FID*3] - *(camera + nP*6+3);
		x(1) = m_XYZ[FID*3+1] - *(camera + nP*6+4);
		x(2) = m_XYZ[FID*3+2] - *(camera + nP*6+5);
	}
	
	else//如果特征点的3D坐标未知，则估计特征点的3D坐标，注意用三个角度表示特征点的3D坐标
	{
		double *ptr = m_KR + nP*9;	//m_KR：内参矩阵x旋转矩阵；np是当前3D点的序号，*9是间隔
		Matrix3d  A;
		A << ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7], ptr[8];//输入初始化的KR矩阵
		//从 KR 中提取 3x3 矩阵（A）
		double matx[3];				//齐次的像点坐标
		matx[0] = imgpts[0];
		matx[1] = imgpts[1];
		matx[2] = 1;
		//Matrix A是KR矩阵，而行向量（数组）是2D投影坐标
		Vector3d  b(matx);
		//Matrix A是KR矩阵，而Matrix b是齐次坐标
		//这事矩阵运算就不说了。就是2D投影点通过KR矩阵，可以算出3D点坐标，这时候还没有极坐标呢，都是XYZ过程
		x = A.colPivHouseholderQr().solve(b);//求解Ax=b  解算出特征点相对于主锚点的XYZ坐标
		///*------------------------------------调试用*/
		//printf("%f %f %f\n", x(0), x(1), x(2));
		//printf("%f %f %f\n", ptr[6], ptr[7], ptr[8]);
		///*------------------------------------调试用*/
		//x已经是IJRR的公式（15）
	}

	double* pKR = KR + nP*9;//KR表示内参矩阵与旋转矩阵的乘积，和m_KR不一样嘛？可能KR是传进来的，而m_KR是记录的
	double t = pKR[6]*x(0) + pKR[7]*x(1) + pKR[8]*x(2);	//？不懂
	///*------------------------------------调试用----------------------------------*/
	//printf("%f %f %f\n", pKR[6], pKR[7], pKR[8]);
	///*------------------------------------调试用----------------------------------*/
	//compute azimuth and elevation angle
	double dDAngle = atan2( x(0), x(2) );				//\Psi，也就是Azimuth angle-主锚点
	double dHAngle = atan2( x(1), sqrt(x(0)*x(0)+ x(2)*x(2)) );//\theta，也就是Elevation angle
	//feature就是视差角的坐标
	feature[0] = dDAngle;
	feature[1] = dHAngle;
	feature[2] = 0;

	if ( t < 0 )//？
		return true;
	else
		return false;
}

/*初始化辅锚点*/
bool PBA::pba_initializeAssoArchor( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int nMI,
	int nAI,
	int FID,
	bool bLast )
{
	/*
	bool bT = pba_initializeAssoArchor( 
					projs+feastart*2,	//记录每个投影点的xy，2个这样记录&存储, feastart一开始是0，间隔为2
					nphoto+feastart,	//记录每个点的序号，feastart一开始是0，间隔为1
					m_motstruct,
					m_K,
					m_motstruct+m_ncams*cnp+ptno*3,
					0,
					1,
					ptno,
					bLast );
	*/
	int nM = photo[nMI];                              //主锚点号
	int nA = photo[nAI];                              //副锚点号

	Vector3d  xM, xA;

	if (m_bProvideXYZ)
	{
		xM[0] = m_XYZ[FID*3]   - *(camera + nM*6+3);
		xM[1] = m_XYZ[FID*3+1] - *(camera + nM*6+4);
		xM[2] = m_XYZ[FID*3+2] - *(camera + nM*6+5);

		xA[0] = m_XYZ[FID*3]   - *(camera + nA*6+3);
		xA[1] = m_XYZ[FID*3+1] - *(camera + nA*6+4);
		xA[2] = m_XYZ[FID*3+2] - *(camera + nA*6+5);
	}
	else
	{
		//Main anchor ray
		double *ptr1 = m_KR + nM*9;
		Matrix3d  AM;	
		AM << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

		double matxM[3];
		matxM[0] = *(imgpts+2*nMI);
		matxM[1] = *(imgpts+2*nMI+1);
		matxM[2] = 1;

		Vector3d  bM(matxM);
		xM = AM.colPivHouseholderQr().solve(bM);			//特征点相对于主锚点的XYZ坐标

		//Associate archor ray
		double *ptr2 = m_KR + nA*9;
		Matrix3d  AA;	
		AA << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

		double matxA[3];
		matxA[0] = *(imgpts+2*nAI);
		matxA[1] = *(imgpts+2*nAI+1);
		matxA[2] = 1;

		Vector3d  bA(matxA);
		xA = AA.colPivHouseholderQr().solve(bA);			//特征点相对于辅锚点的XYZ坐标
	}
		
	//Parallax Angle
	double dDot = xM(0)*xA(0) + xM(1)*xA(1) + xM(2)*xA(2);			//IJRR的公式（3.13）的分子

	double dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );		//IJRR的公式（3.13）的分母
	double dDisA = sqrt( xA(0)*xA(0)+xA(1)*xA(1)+xA(2)*xA(2) );		//IJRR的公式（3.13）的分母

	if (dDot/(dDisM*dDisA)>1)			//特殊情况
		feature[2] = 0;
	else if (dDot/(dDisM*dDisA)<-1)		//特殊情况
		feature[2] = PI;
	else
	{
		double dw = acos( dDot/(dDisM*dDisA) );
		feature[2] = dw;
	}
	//这个就是/omega了

	
	double pti2k[3];
	//+3是X，也就是副锚点-主锚点，也就是t^m_a，也就是从主锚点m到辅锚点的向量（X方向）XA-XM
	pti2k[0] = *(camera + nA*6+3) - *(camera + nM*6+3);    
	//YA-YM（Y方向）
	pti2k[1] = *(camera + nA*6+4) - *(camera + nM*6+4);    
	//ZA-ZM（Z方向）
	pti2k[2] = *(camera + nA*6+5) - *(camera + nM*6+5);    
	//所以现在pti2k就是向量t_m->a

	//dDot1就是主锚点到特征店F_j（x）即xM和t_m->a的点积
	double dDot1 = xM[0]*pti2k[0] + xM[1]*pti2k[1] + xM[2]*pti2k[2];
	double dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
	double tmp = dDot1/(dDisM*dDisi2k);//其实就是IJRR图2的/varphi的cos
	double dW2;
	if (tmp > 1)						//特殊情况
		dW2 = 0;
	if ( tmp < -1)						//特殊情况
		dW2 = PI;
	else
		dW2  = acos( tmp );//其实就是IJRR图2的/varphi

	return true;//dW2是/varphi_j;dw=feature[2]是/omega_j;但是似乎不需要dW2，或者说是tmp中间值（IJRR的BA公式中）
}

/*初始化其他锚点*/
bool PBA::pba_initializeOtheArchors( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int* archorSort,
	int nfeacout,
	int nOI,
	int FID )
{
	/*
	* bAdjust = pba_initializeOtheArchors( 
					projs+feastart*2,
					nphoto+feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams*cnp + ptno * 3,
					ptr2,
					sum,
					i,
					ptno );
	*/
	static int i = 0;
	double dw = feature[2];                   //视差角
	double dwNew;
	double dmaxw = dw;
	int   nNewI = 0;
	bool bAdjust = false;
	double dDot,dDisM,dDisA;

	if ( dw < MAXARCHOR   )
	{
		//current archor vector 
		int nO = photo[nOI];

		Vector3d  xO;
		if ( m_bProvideXYZ)
		{
			xO(0) = m_XYZ[FID*3] - *(camera + nO*6 + 3);
			xO(1) = m_XYZ[FID*3+1]-*(camera + nO*6 + 4);
			xO(2) = m_XYZ[FID*3+2]-*(camera + nO*6 + 5);
		}
		else
		{
			double *ptr1 = m_KR + nO*9;
			Matrix3d  AO;	
			AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

			double matxO[3];
			matxO[0] = *(imgpts+nOI*2);
			matxO[1] = *(imgpts+nOI*2+1);
			matxO[2] = 1;
			///*-------------------------------调试用--------------------------------*/
			//printf("%f %f %f\n", matxO[0], matxO[1], matxO[2]);
			///*-------------------------------调试用--------------------------------*/
			Vector3d  bO(matxO);
			xO = AO.colPivHouseholderQr().solve(bO);
		}

		double dDAngle = atan2( xO(0), xO(2) );
		double dHAngle = atan2( xO(1), sqrt(xO(0)*xO(0)+ xO(2)*xO(2)) );

		for ( i = 0; i < nfeacout; i++ )
		{
			//Main Archor Vector
			int nM = photo[i];                             //像片号
			Vector3d  xM;

			if (m_bProvideXYZ)
			{
				xM(0) = m_XYZ[FID*3]  -*(camera + nM*6+3);
				xM(1) = m_XYZ[FID*3+1]-*(camera + nM*6+4);
				xM(2) = m_XYZ[FID*3+2]-*(camera + nM*6+5);
			}
			else
			{
				double *ptr2 = m_KR + nM*9;
				Matrix3d  AM;	
				AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

				double matxM[3];
				matxM[0] = *(imgpts+i*2);
				matxM[1] = *(imgpts+i*2+1);
				matxM[2] = 1;

				Vector3d  bM(matxM);
				xM = AM.colPivHouseholderQr().solve(bM);
			}

			//Parallax angle between current archor and main archor
			dDot = xM(0)*xO(0) + xM(1)*xO(1) + xM(2)*xO(2);
			dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );
			dDisA = sqrt( xO(0)*xO(0)+xO(1)*xO(1)+xO(2)*xO(2) );

			if( dDot/(dDisM*dDisA) > 1 )
				dwNew = 0;
			else if(dDot/(dDisM*dDisA)<-1)
				dwNew = PI;
			else
				dwNew = acos( dDot/(dDisM*dDisA) );

			if ( dwNew > dmaxw )
			{
				dmaxw = dwNew;
				archorSort[0] = nOI;   //存放主锚点号
				archorSort[1] = i;     //存放副锚点号
				feature[0] = dDAngle;
				feature[1] = dHAngle;
				feature[2] = dmaxw;
				bAdjust = true;
			}
		}
	}
	return bAdjust;
}


bool PBA::pba_initializeOtheArchors_Mindw(
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int* archorSort,
	int nfeacout,
	int nOI,
	int FID)
{
	/*
	* bAdjust = pba_initializeOtheArchors(
					projs+feastart*2,
					nphoto+feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams*cnp + ptno * 3,
					ptr2,
					sum,
					i,
					ptno );
	*/
	static int i = 0;
	double dw = feature[2];                   //视差角
	double dwNew;
	double dminw = dw;
	int   nNewI = 0;
	bool bAdjust = false;
	double dDot, dDisM, dDisA;

	if (dw < MAXARCHOR)
	{
		//current archor vector 
		int nO = photo[nOI];

		Vector3d  xO;
		if (m_bProvideXYZ)
		{
			xO(0) = m_XYZ[FID * 3] - *(camera + nO * 6 + 3);
			xO(1) = m_XYZ[FID * 3 + 1] - *(camera + nO * 6 + 4);
			xO(2) = m_XYZ[FID * 3 + 2] - *(camera + nO * 6 + 5);
		}
		else
		{
			double* ptr1 = m_KR + nO * 9;
			Matrix3d  AO;
			AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

			double matxO[3];
			matxO[0] = *(imgpts + nOI * 2);
			matxO[1] = *(imgpts + nOI * 2 + 1);
			matxO[2] = 1;
			///*-------------------------------调试用--------------------------------*/
			//printf("%f %f %f\n", matxO[0], matxO[1], matxO[2]);
			///*-------------------------------调试用--------------------------------*/
			Vector3d  bO(matxO);
			xO = AO.colPivHouseholderQr().solve(bO);
		}

		double dDAngle = atan2(xO(0), xO(2));
		double dHAngle = atan2(xO(1), sqrt(xO(0) * xO(0) + xO(2) * xO(2)));

		for (i = 0; i < nfeacout; i++)
		{
			//Main Archor Vector
			int nM = photo[i];                             //像片号
			Vector3d  xM;

			if (m_bProvideXYZ)
			{
				xM(0) = m_XYZ[FID * 3] - *(camera + nM * 6 + 3);
				xM(1) = m_XYZ[FID * 3 + 1] - *(camera + nM * 6 + 4);
				xM(2) = m_XYZ[FID * 3 + 2] - *(camera + nM * 6 + 5);
			}
			else
			{
				double* ptr2 = m_KR + nM * 9;
				Matrix3d  AM;
				AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

				double matxM[3];
				matxM[0] = *(imgpts + i * 2);
				matxM[1] = *(imgpts + i * 2 + 1);
				matxM[2] = 1;

				Vector3d  bM(matxM);
				xM = AM.colPivHouseholderQr().solve(bM);
			}

			//Parallax angle between current archor and main archor
			dDot = xM(0) * xO(0) + xM(1) * xO(1) + xM(2) * xO(2);
			dDisM = sqrt(xM(0) * xM(0) + xM(1) * xM(1) + xM(2) * xM(2));
			dDisA = sqrt(xO(0) * xO(0) + xO(1) * xO(1) + xO(2) * xO(2));

			if (dDot / (dDisM * dDisA) > 1)
				dwNew = 0;
			else if (dDot / (dDisM * dDisA) < -1)
				dwNew = PI;
			else
				dwNew = acos(dDot / (dDisM * dDisA));

			if (dwNew < dminw)
			{
				dminw = dwNew;
				archorSort[0] = nOI;   //存放主锚点号
				archorSort[1] = i;     //存放副锚点号
				feature[0] = dDAngle;
				feature[1] = dHAngle;
				feature[2] = dminw;
				bAdjust = true;
			}
		}
	}
	return bAdjust;
}

//--------------------------------------------------------
int PBA::pba_angle2xytGN( double *p )
{
	static int i, j;
	double* pAngle;
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	static double aa;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		pAngle = p + m_ncams*6 + i*3;
		w = *(p + m_ncams*6 + i*3 + 2);

		if ( i == 408741 )
		{
			aa = w;
		}

		xj[0] = sin(*(pAngle)) * cos(*(pAngle+1));
		xj[1] = sin(*(pAngle+1));
		xj[2] = cos(*(pAngle)) * cos(*(pAngle+1));

		nM = *(m_archor+i*3+1);
		nA = *(m_archor+i*3+2);

		Tik[0] = -*(p+nM*6+3)+*(p+nA*6+3);
		Tik[1] = -*(p+nM*6+4)+*(p+nA*6+4);
		Tik[2] = -*(p+nM*6+5)+*(p+nA*6+5);
		
		Dik = sqrt( Tik[0]*Tik[0] + Tik[1]*Tik[1] + Tik[2]*Tik[2] );

		w2 = acos( (xj[0]*Tik[0]+xj[1]*Tik[1]+xj[2]*Tik[2])/Dik );

		xk[0] = ( Dik * sin(w2+w) * xj[0] )/sin(w);
		xk[1] = ( Dik * sin(w2+w) * xj[1] )/sin(w);
		xk[2] = ( Dik * sin(w2+w) * xj[2] )/sin(w);

		*(p+m_ncams*6+i*3)   = *(p+nM*6+3) + xk[0];
		*(p+m_ncams*6+i*3+1) = *(p+nM*6+4) + xk[1];
		*(p+m_ncams*6+i*3+2) = *(p+nM*6+5) + xk[2];
	}

	return 1;
}

void PBA::pba_saveXYZ( const char* camera, const char* sz3Dpt, double *p,bool gn  )
{
	static int i = 0;
	double dx, dy, dz;
	FILE *fp = NULL, *fpc = NULL;
	
	//transform from angle to xyz
	if( gn )
		pba_angle2xytGN(p);
	else
		pba_angle2xytLM(p);
	
	//save camera poss
	if( camera != NULL )
	{
		fpc = fopen( camera, "w" );

		for( i = 0; i < m_ncams; i++ )
		{
			fprintf( fpc, "%0.5lf     %0.5lf      %0.5lf     %0.5lf       %0.5lf     %0.5lf\n",
				*(p+i*6), *(p+i*6+1), *(p+i*6+2), *(p+i*6+3), *(p+i*6+4), *(p+i*6+5) );
		}
		fclose(fpc);
	}	

	//save features xyz
	if ( sz3Dpt!= NULL )
	{
		fp	= fopen( sz3Dpt, "w" );
		for ( i = 0; i < m_n3Dpts; i++ )
		{
			dx = *(p+m_ncams*6+i*3);
			dy = *(p+m_ncams*6+i*3+1);
			dz = *(p+m_ncams*6+i*3+2);
			fprintf( fp, "%0.5lf     %0.5lf     %0.5lf\n", dx, dy, dz );
		}
		fclose(fp);
	}

}

int	PBA::pba_angle2xytLM( double *p )
{
	static int i, j;
	double* pAngle;
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	int cnp = 6, pnp = 3;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		pAngle = p + m_ncams*6 + i*3;
		w = *(p + m_ncams*6 + i*3 + 2);

		xj[0] = sin(*(pAngle)) * cos(*(pAngle+1));//x
		xj[1] = sin(*(pAngle+1));//y
		xj[2] = cos(*(pAngle)) * cos(*(pAngle+1));//z

		nM = *(m_archor+i*3+1);
		nA = *(m_archor+i*3+2);

		Tik[0] = -*(p+nM*6+3)+*(p+nA*6+3);
		Tik[1] = -*(p+nM*6+4)+*(p+nA*6+4);
		Tik[2] = -*(p+nM*6+5)+*(p+nA*6+5);

		Dik = sqrt( Tik[0]*Tik[0] + Tik[1]*Tik[1] + Tik[2]*Tik[2] );

		w2 = acos( (xj[0]*Tik[0]+xj[1]*Tik[1]+xj[2]*Tik[2])/Dik );

		//特征点相对于主锚点的xyz坐标
		xk[0] = ( Dik * sin(w2+w) * xj[0] )/sin(w);
		xk[1] = ( Dik * sin(w2+w) * xj[1] )/sin(w);
		xk[2] = ( Dik * sin(w2+w) * xj[2] )/sin(w);

		//特征点的全局xyz坐标
		*(p+m_ncams*6+i*3)   = *(p+nM*6+3) + xk[0];
		*(p+m_ncams*6+i*3+1) = *(p+nM*6+4) + xk[1];
		*(p+m_ncams*6+i*3+2) = *(p+nM*6+5) + xk[2];
	}

	return 1;
}
