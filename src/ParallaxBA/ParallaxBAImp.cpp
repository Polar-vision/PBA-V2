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
#include "ParallaxBAImp.h"

#define MAXSTRLEN  2048 /* 2K */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
	while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}

//cite from SBA code in order to search sub matrix of S according to (i,j)
static void sba_crsm_alloc(struct sba_crsm *sm, int nr, int nc, int nnz)
{
	int msz;
	sm->nr=nr;
	sm->nc=nc;
	sm->nnz=nnz;
	msz=2*nnz+nr+1;
	sm->val=(int *)malloc(msz*sizeof(int));  /* required memory is allocated in a single step */
	if(!sm->val){
		fprintf(stderr, "memory allocation request failed in sba_crsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
		exit(1);
	}
	sm->colidx=sm->val+nnz;
	sm->rowptr=sm->colidx+nnz;
}

static void sba_crsm_free(struct sba_crsm *sm)
{
	 sm->nr=sm->nc=sm->nnz=-1;
	free(sm->val);
	sm->val=sm->colidx=sm->rowptr=NULL;
}


/* returns the index of the (i, j) element. No bounds checking! */
static int sba_crsm_elmidx(struct sba_crsm *sm, int i, int j)//返回i行j列元素的列索引
{
	register int low, high, mid, diff;

	low=sm->rowptr[i];
	high=sm->rowptr[i+1]-1;

	/* binary search for finding the element at column j */
	while(low<=high)
	{
		mid=(low+high)>>1; //(low+high)/2;
		diff=j-sm->colidx[mid];
		if(diff<0)//在j右侧
			 high=mid-1;
		else if(diff>0)//在j左侧
			low=mid+1;
		else
		return mid;
	}

	return -1; /* not found */
}


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

CParallaxBA::CParallaxBA(void)
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


CParallaxBA::~CParallaxBA(void)
{
	if ( m_szCameraInit!= NULL)	free(m_szCameraInit);
	if ( m_szFeatures!= NULL)	free(m_szFeatures);
	if ( m_szCalibration!= NULL)	free(m_szCalibration);
	if ( m_szXYZ!= NULL)		free(m_szXYZ);
	if ( m_szCamePose!= NULL)	free(m_szCamePose);
	if ( m_sz3Dpts!= NULL)		free(m_sz3Dpts);
	if ( m_szReport!= NULL)		free(m_szReport);
}

bool CParallaxBA::pba_run(bool bRobust,
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

	pba_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );     //给定初值

	if (m_bsolverGN)
		pba_motstr_gn( );
	else
	{
		//double t1 = clock();
		pba_motstr_levmar();
		//double t2 = clock();
		//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
		//printf("%s %f\n", "不包含初始化内存的平差时间：", t_diff);
	}
	return true;
}

bool CParallaxBA::pba_run( int argc, char** argv )
{
	bool bTrue = pba_parseArgs( argc, argv );
	if (!bTrue)
	{
		fprintf( stderr, "ParallaxBA: Input wrong commands, please check them!\n");
		return false;
	}
	
	pba_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );

	if (m_bsolverGN)		
		pba_motstr_gn( );
	else					
	{
		pba_motstr_levmar();
	}
	return true;
}


bool CParallaxBA::pba_initialize( char* szCamera, char* szFeature,  char* szCalib, char* szXYZ )
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
		m_K = (double*)malloc(9*sizeof(double));                           //存放相机内方位元素
		pba_readCameraPoseration(szCalib, m_K);                            //
		///*-------------------------------------------------调试用---------------------------------------------------*/
		//for (int i = 0;i < 9;i++)
		//{
		//	printf("%f ", m_K[i]);
		//	if ((i + 1) % 3 == 0)
		//		printf("\n");
		//}
		///*-------------------------------------------------调试用---------------------------------------------------*/
	}	

	if ( szXYZ != NULL )
		m_bProvideXYZ = true;	
	zu = 0;
	savePara = 0;//Add by Lu
	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	pba_readAndInitialize( szCamera, szFeature, &m_ncams, &m_n3Dpts, &m_n2Dprojs,&m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
				 &m_imgpts, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort );
	
	if(savePara == 1)//Add by Lu
		pba_saveInitialParallax("C:/txtdata/Village/InitialParallax.txt", m_motstruct);//Add by Lu
	//Modify the directory to export point file
	if (zu == 1)
		pba_saveTriangulatedxyz("C:/txtdata/Village/Triangulatedxyz.txt", m_motstruct);
	else if (zu == 2)//Add by Lu
	{
		pba_saveInitialXYZ("C:/txtdata/Taian/xyz_Maxdw.txt", pba_angle2xyz(m_motstruct));//转换为xyz
	}
	else//Add by Lu
		;//Other codes

	printf( "Number of cameras: %d\n", m_ncams );
	printf( "Number of points: %d\n", m_n3Dpts );
	printf( "Number of projections: %d\n", m_n2Dprojs );
	return true;
}
double *CParallaxBA::pba_angle2xyz(double* p)
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
void CParallaxBA::pba_saveInitialXYZ(char* sz3Dpt, double* p)
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
void CParallaxBA::pba_saveInitialParallax(char* sz3Dpt, double* p)
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
void CParallaxBA::pba_saveTriangulatedxyz(char* sz3Dpt, double* p)
{
	static int i = 0;
	double x, y, z;
	FILE* fp = NULL;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		fp = fopen(sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			x = *(p + m_ncams * 6 + i * 3);
			y = *(p + m_ncams * 6 + i * 3 + 1);
			z = *(p + m_ncams * 6 + i * 3 + 2);
			//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1, x, y, z);
			fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", x, y, z);
			//fprintf(fp, "%0.5lf\n", parallax * 180 / PI);
		}
		fclose(fp);
	}
	free(p);
}

bool CParallaxBA::pba_motstr_levmar( )
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

	nuis = pba_ConstructSmask( Sidxij, Uidxij );
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
	pba_constructAuxCSSLM( Ap, Aii );

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
		//第一步.1就是计算雅克比矩阵
		//计算出的雅各布矩阵 U、V、W ：代表误差相对于摄像机参数和 3D 点的部分导数？
		if ( m_bRobustKernel)
			pba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			pba_jacobian(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		t1 = clock();
		//第一步.2就是确定调整阻尼参数，阻尼参数控制步长。
		if( itno == 0)//iter No，第一次就是0就肯定会执行
			mu = pba_computeInitialmu( U, V, Uidxij, tau, nvars );

		nLmIterations = 0;//迭代从0开始
		while(1) //determine increment using adaptive damping ，这就是上面的pba_computeInitialmu()函数
		{//雅各布矩阵表示重投误差相对于参数的导数。
			nLmIterations++;
			t11 = clock();
			//求V的逆矩阵V^-1
			pba_inverseVLM( V, IV, Uidxij, mu ); //compute inverse matrix of V
			
			//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
			memset( E, 0, m*easz*sizeof(double));
			memset( S, 0, m_nS*Sblsz*sizeof(double) );
			pba_constructSLM( S, E, U, IV, W, ea, eb, Sidxij, mu );
			t2 = clock();
			
			//啥玩意？
			//Solve linear equation
			//set CSS format using S matrix
			pba_constructCSSLM( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init ); 
			for ( ii = 0; ii < cnp*m; ii++  )
				Ex[ii] = E[ii];	
			
			pba_solveCholmodLM( Ap, Aii, init, ordering);
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
				pba_solveFeatures( W, IV, ea, eb, dpa, dpb );
			
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

				pba_updateKR( m_KR, m_KdA, m_KdB, m_KdG, m_K, pdp );
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

//#pragma optimize("", off) // 关闭优化
bool CParallaxBA::pba_motstr_gn( FixType ft )
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
	struct sba_crsm Sidxij, Uidxij;		// S mask and U mask          存储稀疏矩阵
	nuis = pba_ConstructSmask( Sidxij, Uidxij );   //返回U矩阵非零元素个数

	Usz=cnp*cnp; Vsz=pnp * pnp; Wsz=cnp * pnp; 	
	esz=mnp; easz=cnp; ebsz=pnp; Sblsz=cnp*cnp;	Sdim=m * cnp;	
	nvis = m_n2Dprojs;                                            //像点个数
	nobs=nvis*mnp;                                                //像点坐标个数，观测方程个数
	nvars=m*cnp + n*pnp;                                          //未知参数个数

	S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));          //S矩阵：非0元素个数*6*6
	W	=	(double *)malloc(nvis*Wsz*sizeof(double));            //W矩阵：像点数*6*3
	U	=	(double *)malloc(nuis*Usz*sizeof(double));            //U矩阵：非0元素个数×6*6
	V	=	(double *)malloc(n*Vsz*sizeof(double));               //V矩阵：三维点数×3*3
	e	=	(double *)malloc(nobs*sizeof(double));                //误差项 
	eab	=	(double *)malloc(nvars*sizeof(double));               //未知参数
	E	=	(double *)malloc(m*cnp*sizeof(double));	              //外方位元素	
	dp	=	(double *)malloc(nvars*sizeof(double));               //未知参数
	hx	=	(double *)malloc(nobs*sizeof(double));	              //观测方程
	pa=p; //pa指向未知参数
	pb=p+m*cnp; //pb指向三维点参数

	ea=eab;//ea指向
	eb=eab+m*cnp;

	dpa=dp; 
	dpb=dp+m*cnp;
	
	//Select fix axis for Gauss-Newton
	//由于 BA 只关心相机和点 相对位置，如果不固定某些参数，优化过程中可能会发生 自由漂移 (Drift)，即：整个场景可以任意平移、旋转、缩放，而不会改变投影误差（导致解不唯一）。固定第一帧和最后一帧的轴？
	if ( ft == PBA_FixDefault )
	{
		//double test1 = *(m_motstruct + 9);
		//double test2 = *(m_motstruct + 3);
		//double test3 = *(m_motstruct + 10);
		//double test4 = *(m_motstruct + 4);
		//double test5 = *(m_motstruct + 11);
		//double test6 = *(m_motstruct + 5);
		tmp = abs(*(m_motstruct+9)-*(m_motstruct+3));//第2个相机的Xc-第1个相机的Xc
		ft = PBA_FixX;									
		tmp1 = abs(*(m_motstruct+10)-*(m_motstruct+4));//第2个相机的Yc-第1个相机的Yc
		if ( tmp1 > tmp )
		{ft= PBA_FixY; tmp = tmp1;	}
		tmp1 = abs(*(m_motstruct+11)-*(m_motstruct+5));//第2个相机的Zc-第1个相机的Zc
		if ( tmp1 > tmp )
		{ft= PBA_FixZ; }
	}
	nft = ft;
	
	// 使用 CHOLMOD 进行稀疏矩阵因子分解（Cholesky 分解）的准备工作
	cholmod_start (&m_cS) ;  	//初始化 CHOLMOD 内部的数据结构，主要是分配工作空间，并设置默认参数
	//m_cS.print_function = NULL;	
	//在解算方程的时候用到
	int *Ap  = (int*)malloc((m)*sizeof(int));
	int * Aii = (int*)malloc(m_nS*sizeof(int));//m_nS，S矩阵中非零元素总数
	pba_constructAuxCSSGN( Ap, Aii );//构造稀疏矩阵的索引结构，从1行1列开始，存储每列的起始索引,存储非零元素的行索引

	//在 CHOLMOD 里，m_cholSparseE 是一个 cholmod_dense 结构体
	//typedef struct cholmod_dense
	//{
	//	size_t nrow;  // 行数
	//	size_t ncol;  // 列数
	//	size_t nzmax; // 最大非零元素数（通常是 nrow * ncol）
	//	void* x;      // 指向数据的指针
	//} cholmod_dense;

	m_cholSparseE = cholmod_zeros( cnp*m-7, 1, CHOLMOD_REAL, &m_cS);//固定第1帧的外参和第2帧的Xc轴，实际未知参数个数是6*m-7
	Ex = (double*)m_cholSparseE->x;
	
	//(m_nS-m_ncams)表示S矩阵中非对角块（相机之间的约束）的数量，其中的每个非对角块占据6*6个非零元素（不对称）
	//m_ncams表示S矩阵中对角块（单个相机自身的约束）的数量，其中每个对角块占据21个非零元素（对称，只存储下三角部分）
	nMaxS = (m_nS-m_ncams)*36+m_ncams*21;	//maximum non-zero element in S matrix 

	//分配稀疏矩阵S，第4个参数为true表示按列索引排序，第5个参数为true代表紧密存储，第6个参数为1代表矩阵的存储类型为实数。
	//紧密存储（Packed Storage）指的是只存储矩阵中非零元素，并且所有数据连续存放在内存中，不预留额外空间。这样可以减少内存占用，提高计算效率。
	m_cholSparseS = cholmod_allocate_sparse(m_ncams*6-7,m_ncams*6-7,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_cholSparseS->p;		//column pointer
	Si = (int*)m_cholSparseS->i;		//row pointer
	
	//Compute initial error and initial reprojection error 
	tStart = clock();
	dpa[0] = dpa[1] = dpa[2] = dpa[3] = dpa[4] = dpa[5] = 0;//未用到
	dpa[9+nft] = 0;
	pba_cost(p, hx, m_archor );

	//double sum = 0;
	//for (int i = 0; i < nvis; i+=10000)
	//{
	//	printf("%f %f    %f %f\n", hx[2 * i], hx[2 * i + 1], x[2 * i], x[2 * i + 1]);
	//	sum += (pow(hx[2 * i] - x[2 * i], 2) + pow(hx[2 * i + 1] - x[2 * i + 1], 2));
	//}
	
	p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
	//printf("%f %f\n", sum, p_eL2);

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
		memset( U, 0, nuis*Usz*sizeof(double) );//nuis是U矩阵的非零元素个数，Usz=6*6，存储U矩阵的非零元素
		memset( ea, 0,m*easz*sizeof(double) );//m是相机数量，easz=6，ea存储相机参数
		memset( V, 0, n*Vsz*sizeof(double));//n是三维点数量，Vsz=3*3，
		memset( eb, 0, n*ebsz*sizeof(double));//n是三维点数量，ebsz=3，eb存储三维点参数
		memset( W, 0, nvis*Wsz*sizeof(double));	//nvis是像点数，Wsz=18，每行的非零元素（每行9个非零元素）

		//Step one; compute W V U directly, don't save each projection image Jacobian 
		t0 = clock();
		if ( m_bRobustKernel)
			pba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			pba_jacobian(p, m_archor, &Uidxij, e, U, ea, V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature);

		t1 = clock();

		pba_inverseVGN( U, V, Uidxij ); //compute inverse matrix of V

		//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
		memset( E, 0, m*easz*sizeof(double));
		memset( S, 0, m_nS*Sblsz*sizeof(double) );
		pba_constructSGN( S, E, U, V, W, ea, eb, Sidxij );//S为方程左侧，E为方程右侧
		t2 = clock();
	
		//Solve equation
		pba_constructCSSGN( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init, nft ); //set CSS format using S matrix//break
		for ( ii = 0; ii < 3+nft; ii++ )
			Ex[ii] = E[6+ii];//第2帧的姿态角

		for ( ii = 0; ii < 6*m-(10+nft); ii++  )
			Ex[3+nft+ii] = E[10+nft+ii];//跳过第2帧的Xc，只存储Yc和Zc

		//Sx为左侧，Ex为右侧

		pba_solveCholmodGN( Ap, Aii, init, ordering);//break
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
		pba_solveFeatures( W, V, ea, eb, dpa, dpb );//break
		
		//update
		for(i=0, dp_L2=0.0; i<nvars; ++i)//nvars：未知参数个数
		{
			p[i]=p[i] + (tmp=dp[i]);//dpa指向dp中的外参改正量，dpb指向dp中的三维点改正量，p指向未知参数
			dp_L2+=tmp*tmp;//存储改正量的平方和
		}

		//dp_L2=sqrt(dp_L2);
		if (dp_L2<1E-8*1E-8)//break
		{	stop = 1;	goto iterstop;	}	//参数的变化已经很小						

		pba_updateKR( m_KR, m_KdA, m_KdB, m_KdG, m_K, p );
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

double CParallaxBA::pba_computeInitialmu(double* U, double* V, sba_crsm& Uidxij, double tau, int nvars )
{
	int i, j;
	int pos, m = m_ncams,n = m_n3Dpts, cnp = 6, pnp = 3, Usz = 36, Vsz = 9;
	double tmp = 0;
	double *ptr1, *ptr2;
	double mu;

	double *diagUV = (double *)malloc(nvars*sizeof(double));

	double *diagU = diagUV; 
	double *diagV = diagUV + m*cnp;

	for(j=0; j<m; ++j)
	{
		pos  = sba_crsm_elmidx(&Uidxij, j, j);
		ptr1 = U + pos*Usz;
		ptr2 = diagU + j*cnp; 
		for(i=0; i<cnp; ++i)
			ptr2[i]=ptr1[i*cnp+i];
	}
	for(i=0; i<n; ++i)
	{
		ptr1=V + i*Vsz; // set ptr1 to point to V_i
		ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
		for(j=0; j<pnp; ++j)
			ptr2[j]=ptr1[j*pnp+j];
	}	

	/* find max diagonal element */
	for(i=0, tmp=DBL_MIN; i<m*cnp; ++i)
		if(diagUV[i]>tmp) 
			tmp=diagUV[i];
	for(i=m*cnp; i<nvars; ++i) /* tmp is not re-initialized! */
		if(diagUV[i]>tmp) 
			tmp=diagUV[i];

	mu=m_Tau*tmp;
	
	free(diagUV);
	diagU = diagV = NULL;

	return mu;
}

void CParallaxBA::pba_inverseVLM( double* V, double* IV, sba_crsm& Uidxij, double mu )
{
	int i, j;
	int m = m_ncams, n = m_n3Dpts;
	int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
	double *ptr1, *ptr2;
	Matrix3d MatInv;
	
	//IV save inverse V matrix, V must unchange for the next step
	memcpy( IV, V, n*Vsz*sizeof(double) );	
	for(i=0; i<n; ++i)
	{
		ptr1=V + i*Vsz; 
		ptr2=IV+ i*Vsz;

		for(j=0; j<pnp; ++j)
			ptr2[j*pnp+j] += mu;

		Eigen::Matrix3d matV(ptr2);
		MatInv = matV.inverse();
		ptr2[0] = MatInv(0,0);
		ptr2[4] = MatInv(1,1);
		ptr2[8] = MatInv(2,2);
		ptr2[1] = ptr2[3] = MatInv(0,1);
		ptr2[2] = ptr2[6] = MatInv(0,2);
		ptr2[5] = ptr2[7] = MatInv(1,2);
	} 
}

void CParallaxBA::pba_inverseVGN( double* U, double* V, sba_crsm& Uidxij )
{
	int i;
	int m = m_ncams, n = m_n3Dpts;
	int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
	double *ptr1;
	Matrix3d matV, MatInv;
	
	//compute V inverse matrix using Eigen that has better performance than Lapack
	for(i=0; i<n; ++i)
	{
		ptr1=V + i*Vsz; // set ptr1 to point to V_i
		matV << ptr1[0],ptr1[1],ptr1[2],ptr1[3],ptr1[4],ptr1[5],ptr1[6],ptr1[7],ptr1[8];		
		MatInv = matV.inverse();
		ptr1[0] = MatInv(0,0);
		ptr1[4] = MatInv(1,1);
		ptr1[8] = MatInv(2,2);
		ptr1[1] = ptr1[3] = MatInv(0,1);
		ptr1[2] = ptr1[6] = MatInv(0,2);
		ptr1[5] = ptr1[7] = MatInv(1,2);
	}	
}

int CParallaxBA::pba_ConstructSmask( sba_crsm& Sidxij, sba_crsm& Uidxij )//分别是S矩阵和U矩阵的稀疏索引
{
	int i, j, k, ii, jj;
	int nuis, m = m_ncams;//nuis是U矩阵的非零元素个数
	//compute total smask
	for ( i = 0; i < m; i++ ) for ( j = 0; j < m; j++ )
	{
		//这段代码未执行
		if ( m_umask[i*m+j] == 1 && m_smask[i*m+j] == 0 )//m_umask[i*m+j] == 1 表示 U 矩阵在该位置是非零的，如果 S 还没标记，就把它标记上。
		{
			m_smask[i*m+j] = 1;//S 矩阵的稀疏掩码 (Sidxij 的结构),m_smask[i*m+j] 是一个 m×m 维的二维数组（用 1D 存储），表示 S 矩阵在 (i,j) 位置是否有非零元素。
			m_nS += 1;
		}
	}	
	//分配 S 矩阵的存储空间,Sidxij是S 矩阵的结构体,m, m为矩阵大小
	sba_crsm_alloc(&Sidxij, m, m, m_nS);//m_nS是S矩阵的非零元素个数
	for(i=k=0; i<m; ++i)
	{
		Sidxij.rowptr[i]=k;// 记录 S 矩阵第 i 行的起始索引
		ii=i*m;
		for(j=0; j<m; ++j)
			if(m_smask[ii+j])// 如果该位置是非零元素
			{
				Sidxij.val[k]=k;// 这里存的值是索引 k
				Sidxij.colidx[k++]=j; // 记录该非零元素的列索引
			}
	}
	Sidxij.rowptr[m]=m_nS;// 最后一行指向末尾

	for(i=nuis=0, jj=m*m; i<jj; ++i)
		nuis+=(m_umask[i]!=0);//U 矩阵的稀疏掩码 (Uidxij 的结构)

	sba_crsm_alloc(&Uidxij, m, m, nuis);
	for(i=k=0; i<m; ++i)
	{
		Uidxij.rowptr[i]=k;
		ii=i*m;
		for(j=0; j<m; ++j)
			if(m_umask[ii+j])
			{
				Uidxij.val[k]=k;
				Uidxij.colidx[k++]=j;
			}
	}
	Uidxij.rowptr[m]=nuis;

	////测试
	//for (i = 0; i < m; ++i)
	//{
	//	for (j = 0; j < m; ++j)
	//	{
	//		printf("%d ", m_umask[i * m + j]);
	//	}
	//	printf("\n");
	//}
	//printf("\n\n\n");
	//for (i = 0; i <= m; ++i)
	//	printf("%d ", Uidxij.rowptr[i]);
	//printf("\n\n\n");
	//for (i = 0; i < nuis; ++i)
	//	printf("%d ", Uidxij.val[i]);
	//printf("\n\n\n");
	//for (i = 0; i < nuis; ++i)
	//	printf("%d ", Uidxij.colidx[i]);
	//printf("\n\n\n");

	////测试
	//for (i = 0; i < m; ++i)
	//{
	//	for (j = 0; j < m; ++j)
	//	{
	//		printf("%d ", m_smask[i * m + j]);
	//	}
	//	printf("\n");
	//}
	//printf("\n\n\n");
	//for (i = 0; i <= m; ++i)
	//	printf("%d ", Sidxij.rowptr[i]);
	//printf("\n\n\n");
	//for (i = 0; i < m_nS; ++i)
	//	printf("%d ", Sidxij.val[i]);
	//printf("\n\n\n");
	//for (i = 0; i < m_nS; ++i)
	//	printf("%d ", Sidxij.colidx[i]);
	//printf("\n\n\n");


	return nuis;

}

void CParallaxBA::pba_constructAuxCSSLM( int *Ap, int *Aii )
{
	int* Cp = Ap;
	int* Ci = Aii;
	int ii, jj;
	int m = m_ncams, nZ = 0;
	for ( ii = 0; ii < m; ii++ ) 
	{
		*Cp = nZ;
		for( jj=0; jj<=ii; jj++ )
		{
			if (m_smask[jj*m+ii]==1)
			{
				*Ci++ = jj;
				nZ++;
			}
		}
		Cp++;
	}
	*Cp=nZ;
}
//构造稀疏矩阵的索引结构，从1行1列开始，存储每列的起始索引,存储非零元素的行索引
void CParallaxBA::pba_constructAuxCSSGN( int *Ap, int *Aii )
{
	//从0行0列开始是因为固定第一帧
	int* Cp = Ap;//存储每列的起始索引
	int* Ci = Aii;//存储非零元素的行索引
	int ii, jj;
	int m = m_ncams, nZ = 0;
	//printf("\n\n\n");
	for ( ii = 1; ii < m; ii++ ) //列
	{
		*Cp = nZ;
		for( jj=1; jj<=ii; jj++ )//行
		{
			if (m_smask[jj*m+ii]==1)//按列遍历S矩阵（对称矩阵）的上三角
			{
				*Ci++ = jj-1;//不能从0行开始
				//printf("%d ", jj-1);
				nZ++;
			}
		}
		Cp++;
	}
	*Cp=nZ;
	
	//测试
	//Cp = Ap;//将指针移到文件头
	//Ci = Aii;
	//printf("\n\n\n");
	//for (int i = 0; i < m; ++i)
	//	printf("%d ", *Cp++);
	//printf("\n\n\n");
	//for (int i = 0; i < m_nS; ++i)
	//	printf("%d ", *Ci++);
	//printf("\n\n\n");
}

void CParallaxBA::pba_constructSLM( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij, double mu )
{
	int i, j, ii, jj, k, l;
	int m = m_ncams, cnp = 6, pnp = 3, Usz = 36, ebsz=3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptrS, *ptrE;
	double WV[6*3], sum;

	//Copy U matrix to S matrix 
	pos = 0;
	for ( i = 0; i < m; i++ ) for( j = 0; j < m; j++ )
	{
		if ( m_umask[i*m+j] == 1)// save upper triangle for diagonal element S
		{
			pos1 = sba_crsm_elmidx( &Sidxij, i, j);
			ptr2 = S + pos1*36;
			if ( i == j )
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=6)
				{
					ptr2[ii]= ptr1[ii*cnp+ii] + mu;
					for(jj=ii+1; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				}
				pos++;
			}
			else
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=6)
					for(jj=0; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				pos++;
			}
		}
	}

	for ( i = 0; i < m*cnp; i++ )
		E[i] = ea[i];

	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		for ( j = 0; j < numfea; j++ )
		{
			nF1 = m_feature[pos];
			nP1 = m_photo[pos];
			memset( WV, 0, sizeof(double)*cnp*3 );

			ptr1 = W + pos*cnp*3;	
			ptr2 = V + nF1*3*3;	
			ptrE = E + nP1*cnp;


			//WV
			for(ii=0; ii<cnp; ++ii)
			{
				ptr3=ptr1+ii*pnp;
				for(jj=0; jj<pnp; ++jj)
				{
					for(k=0, sum=0.0; k<=jj; ++k)
						sum+=ptr3[k]*ptr2[jj*pnp+k]; 
					for( ; k<pnp; ++k)
						sum+=ptr3[k]*ptr2[k*pnp+jj]; 
					for(k=0, sum= 0.0; k<pnp; k++ )
						sum+=ptr3[k]*ptr2[jj*pnp+k];
					WV[ii*pnp+jj]=sum;
				}
			}

			for ( k = j; k < numfea; k++ )
			{
				nP2 = m_photo[pos+(k-j)];

				//W(V^-1)W^T
				ptr3 = W + (pos+(k-j))*cnp*3;
				//ptrS = S + (nP1*m*36) + nP2*cnp;

				if ( nP1 == nP2 )
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=ii; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}else 
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=0; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}
			}
			//-W^tb
			ptr5 = eb + nF1*ebsz;
			for(ii=0; ii<cnp; ++ii)
			{
				ptr4=WV+ii*pnp;
				for(jj=0, sum=0.0; jj<pnp; ++jj)
					sum+=ptr4[jj]*ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
				ptrE[ii]-= sum;
			}
			pos++;					
		}
	}	
}

void CParallaxBA::pba_constructSGN( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij )
{
	int i, j, ii, jj, k, l;
	int m = m_ncams, cnp = 6, pnp = 3, Usz = 36, ebsz=3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptrS, *ptrE;
	double WV[6*3], sum;

	//Copy U matrix to S matrix 
	pos = 0;
	for ( i = 0; i < m; i++ ) for( j = 0; j < m; j++ )
	{
		if ( m_umask[i*m+j] == 1)// save upper triangle for diagonal element S
		{
			pos1 = sba_crsm_elmidx( &Sidxij, i, j);
			//printf("%d ", pos1);
			ptr2 = S + pos1*36;
			if ( i == j )
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=6)
				{
					ptr2[ii]= ptr1[ii*cnp+ii];
					for(jj=ii+1; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				}
				pos++;
			}
			else
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=6)
					for(jj=0; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				pos++;
			}
		}
	}

	for ( i = 0; i < m*cnp; i++ )
		E[i] = ea[i];//外方位元素

	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		for ( j = 0; j < numfea; j++ )
		{
			nF1 = m_feature[pos];
			nP1 = m_photo[pos];
			memset( WV, 0, sizeof(double)*cnp*3 );

			ptr1 = W + pos*cnp*3;	
			ptr2 = V + nF1*3*3;	
			ptrE = E + nP1*cnp;

			//W(V^-1)，V已经变成了V^-1
			for(ii=0; ii<cnp; ++ii)
			{
				ptr3=ptr1+ii*pnp;//W矩阵第ii行
				for(jj=0; jj<pnp; ++jj)
				{
					for(k=0, sum=0.0; k<=jj; ++k)
						sum+=ptr3[k]*ptr2[jj*pnp+k]; 
					for( ; k<pnp; ++k)
						sum+=ptr3[k]*ptr2[k*pnp+jj]; 
					for(k=0, sum= 0.0; k<pnp; k++ )//(V^-1)是对称矩阵
						sum+=ptr3[k]*ptr2[jj*pnp+k];//此处的sum与上面的sum相等
					WV[ii*pnp+jj]=sum;
				}
			}

			for ( k = j; k < numfea; k++ )
			{
				nP2 = m_photo[pos+(k-j)];

				//W(V^-1)W^T
				ptr3 = W + (pos+(k-j))*cnp*3;
				pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
				ptrS = S + pos1*36;
				
				if ( nP1 == nP2 )
				{
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;	//WV矩阵第ii行								
						for(jj=ii; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;//W矩阵第jj行

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; //-W(V^-1)W^T

							ptrS[jj]-=sum; //S-W(V^-1)W^T
						}
					}
				}else 
				{
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=0; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; //方程左侧
						}
					}
				}
			}
			//-W^tb 方程右侧
			ptr5 = eb + nF1*ebsz;//三维点
			for(ii=0; ii<cnp; ++ii)
			{
				ptr4=WV+ii*pnp;//W(V^-1)第ii行
				for(jj=0, sum=0.0; jj<pnp; ++jj)
					sum+=ptr4[jj]*ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];   W(V^-1)eb
				ptrE[ii]-= sum; //ea-W(V^-1)eb
			}
			pos++;					
		}
	}	
}

void CParallaxBA::pba_solveFeatures( double *W, double *IV, double *ea, double *eb, double *dpa, double *dpb)
{
	int i, j, ii, jj, pos, numfea;
	int nP1, cnp = 6, pnp = 3;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5;
	double sum, eb2[6];
	pos = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		ptr1 = eb  + i*3;//指向三维点误差向量
		ptr2 = IV   + i*3*3;//指向V矩阵
		ptr5 = dpb + i*3;//指向三维点未知参数
		memset( eb2, 0, sizeof(double)*cnp );
		numfea = m_archor[i*3];

		for ( j = 0; j < numfea; j++ )
		{
			nP1  = m_photo[pos];//视图编号
			ptr3 = W + pos*cnp*3;//指向W矩阵
			ptr4 = dpa + nP1*cnp;//指向相机未知参数
			//Wta
			for(ii=0; ii<pnp; ++ii)
			{
				for(jj=0, sum=0; jj < cnp; ++jj )
					sum += ptr3[jj*3+ii]*ptr4[jj];
				eb2[ii] += sum;
			}
			pos++;
		}

		//V*(eb-Wta)
		for(ii=0; ii<pnp; ++ii)
		{
			for(jj=0, sum = 0; jj < pnp; jj++ )
				sum += ptr2[ii*3+jj]*(ptr1[jj]-eb2[jj]);
			ptr5[ii] = sum;
		}
	}
}

void CParallaxBA::pba_constructCSSLM( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init)
{
	int ii, jj, jjj, k;
	int pos1, m = m_ncams;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double *ptr5;
	int nZ = 0;

	Sx = (double*)m_cholSparseS->x;
	if ( !init)
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			for ( k = 0; k < 6; k++ )
			{
				*Sp = nZ;
				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (m_smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								*Si++ = jj*6 + jjj;
								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
							{
								*Si++ = jj*6 + jjj;
								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp=nZ;
	}
	else
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			for ( k = 0; k < 6; k++ )
			{
				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (m_smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
								*Sx++ = ptr5[jjj*6+k];
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
								*Sx++ = ptr5[jjj*6+k];
						}
					}
				}
			}
		}
	}
}

void CParallaxBA::pba_constructCSSGN( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init, int nft)
{
	int ii, jj, jjj, k;
	int pos1, m = m_ncams;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double *ptr5;
	int nZ = 0;
	Sx = (double*)m_cholSparseS->x;//存储稀疏矩阵中所有非零元素的实际值
	//printf("\n\n\n\n");
	if ( !init)
	{
		for ( ii = 1; ii < m; ii++ )  //column，第0个视图是参考帧，通常被固定，不作为优化变量
		{
			for ( k = 0; k < 6; k++ )
			{
				*Sp = nZ;//表示稀疏矩阵中每一列的起始位置（非零元素的索引）,0,1,3,6,
				//printf("%d ", *Sp);
				if ((ii*6+k)==(9+nft))//k=3代表Xc，4代表Yc，5代表Zc。nft=0代表固定Xc，1代表固定Yc，2代表固定Zc
					continue;//固定第2帧的Xc作为尺度约束

				for ( jj = 1; jj <= ii; jj++ )	//row
				{
					if ((m_smask[jj*m+ii]==1))
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						//printf("%d ", pos1);
						ptr5 = S + pos1*36;//指向S矩阵对应（jj,ii）的子矩阵数据
						
						if( ii == jj )//对角线块
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								if ( (jj*6+jjj) != (9+nft))//除了jj=1，jjj=3（代表第2个相机的Xc固定）外，
								{
									if (jj * 6 + jjj < 9 + nft)//第2个相机的姿态角
										*Si++ = jj * 6 + jjj - 6;
									else//第2个相机的Yc，Zc
										*Si++ = jj * 6 + jjj - 7;
									
									*Sx++ = ptr5[jjj*6+k];
									//printf("%d ", ptr5[jjj * 6 + k]);
									nZ++;
								}	
							}
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
							{
								if ((jj*6+jjj) != (9+nft) )
								{
									if (jj * 6 + jjj < 9 + nft)
										*Si++ = jj * 6 + jjj - 6;
									else
										*Si++ = jj * 6 + jjj - 7;
																			
									*Sx++ = ptr5[jjj*6+k];
									nZ++;
								}
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp=nZ;
	}
	else
	{
		for ( ii = 1; ii < m; ii++ )  //column
		{
			for ( k = 0; k < 6; k++ )
			{
				if ((ii*6+k)==(9+nft))
					continue;

				for ( jj = 1; jj <= ii; jj++ )	//row
				{
					if ((m_smask[jj*m+ii]==1))
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								if ( (jj*6+jjj) != (9+nft))
									*Sx++ = ptr5[jjj*6+k];
							}
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
							{
								if ((jj*6+jjj) != (9+nft) )
									*Sx++ = ptr5[jjj*6+k];
							}
						}
					}
				}
			}
		}
	}
}

//learning this skill from G2O
bool CParallaxBA::pba_solveCholmodLM( int* Ap, int* Aii, bool init, bool ordering)
{
	int i, j;
	VectorXi scalarPermutation, blockPermutation;

	ordering = true;
	if (!init)
	{
		if (!ordering)
			m_cholFactorS = cholmod_analyze(m_cholSparseS, &m_cS); // symbolic factorization
		else
		{
			// get the ordering for the block matrix
			if (blockPermutation.size() == 0)
				blockPermutation.resize(m_ncams);
			if (blockPermutation.size() < m_ncams) // double space if resizing
				blockPermutation.resize(2*m_ncams);

			// prepare AMD call via CHOLMOD
			cholmod_sparse auxCholmodSparse;
			auxCholmodSparse.nzmax = m_nS;
			auxCholmodSparse.nrow = auxCholmodSparse.ncol = m_ncams;
			auxCholmodSparse.p = Ap;
			auxCholmodSparse.i = Aii;
			auxCholmodSparse.nz = 0;
			auxCholmodSparse.x = 0;
			auxCholmodSparse.z = 0;
			auxCholmodSparse.stype = 1;
			auxCholmodSparse.xtype = CHOLMOD_PATTERN;
			auxCholmodSparse.itype = CHOLMOD_INT;
			auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
			auxCholmodSparse.sorted = 1;
			auxCholmodSparse.packed = 1;
			int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);
			if (! amdStatus) 
				return false;
			

			// blow up the permutation to the scalar matrix
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_cholSparseS->ncol);
			if (scalarPermutation.size() < (int)m_cholSparseS->ncol)
				scalarPermutation.resize(2*m_cholSparseS->ncol);
			size_t scalarIdx = 0;

			for ( i = 0; i < m_ncams; ++i)
			{
				const int &pp = blockPermutation(i);
				int base = (pp==0) ? 0 : pp*6;
				int nCols= 6;

				for ( j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;

			}
			assert(scalarIdx == m_cholSparseS->ncol);

			// apply the ordering
			m_cS.nmethods = 1 ;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;
			m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);
		}
	}

	cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS); 
	m_cholSparseR = cholmod_solve (CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS) ;
			
	return true;
}
//高效求解GN中出现的稀疏线性系统
//使用Cholmod库进行矩阵（Cholesky）分解和求解
//Ap是列指针，Aii是行索引，是CSC格式
bool CParallaxBA::pba_solveCholmodGN( int* Ap, int* Aii, bool init, bool ordering)
{
	int i, j;
	int m = m_ncams;
	VectorXi scalarPermutation, blockPermutation;

	ordering = true;
	if (!init)//首次初始化
	{
		if (!ordering)//不需要手动排序
		{
			m_cS.nmethods = 1;
			m_cS.method[0].ordering = CHOLMOD_AMD; //排序方法设置为Approximately Minimum Degree（近似最小度排序）
			m_cholFactorS = cholmod_analyze(m_cholSparseS, &m_cS); // symbolic factorization，生成分解所需的内部数据结构
		}
		else
		{
			// get the ordering for the block matrix
			if (blockPermutation.size() == 0)
				blockPermutation.resize(m_ncams-1);//用于存储每个块的排列顺序

			// prepare AMD call via CHOLMOD
			cholmod_sparse auxCholmodSparse;//构造cholmod_sparse类型的辅助矩阵auxcholmodSparse，将稀疏矩阵元数据（如Ap和Aii）映射到Cholmod格式
			auxCholmodSparse.nzmax = m_nS;
			auxCholmodSparse.nrow = auxCholmodSparse.ncol = m-1;//跳过第1帧
			auxCholmodSparse.p = Ap;
			auxCholmodSparse.i = Aii;
			auxCholmodSparse.nz = 0;
			auxCholmodSparse.x = 0;
			auxCholmodSparse.z = 0;
			auxCholmodSparse.stype = 1;
			auxCholmodSparse.xtype = CHOLMOD_PATTERN;
			auxCholmodSparse.itype = CHOLMOD_INT;
			auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
			auxCholmodSparse.sorted = 1;
			auxCholmodSparse.packed = 1;
			//AMD排序用于减少稀疏矩阵分解中的填充
			int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);//执行AMD排序，并将结果存储到blockPermutation
			if (! amdStatus) {
				return false;
			}

			// blow up the permutation to the scalar matrix
			//将块排列扩展为标量排列
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_cholSparseS->ncol);
			size_t scalarIdx = 0;

			int a = 0;
			for ( i = 0; i < m_ncams-1; ++i)
			{
				const int &pp = blockPermutation(i);
				int base = (pp==0) ? 0 : pp*6-1;
				int nCols= (pp==0) ? 5 : 6;
				for (j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;
			}
			assert(scalarIdx == m_cholSparseS->ncol);//m_cholSparseS->ncol=90*6-7=533

			// apply the ordering
			m_cS.nmethods = 1 ;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;//设置CHOLMOD使用给定的排列顺序
			m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);//调用cholmod_analyze_p进行符号分析
		}
		init = true;//完成初始化
	}
		
	//Cholmod package for solving sparse linear equation              
	cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS); //对稀疏矩阵进行数值分解
	m_cholSparseR = cholmod_solve (CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS) ;//求解线性方程组

	return true;
}
//--------------------------------------------------------

//计算文件中非注释行的数量
int CParallaxBA::findNcameras(FILE *fp)
{
	int lineno, ncams, ch;

	lineno=ncams=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#'){ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);

		SKIP_LINE(fp);
		++lineno;
		if(ferror(fp))
		{
			fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		++ncams;
	}
	return ncams;
}

int CParallaxBA::countNDoubles(FILE *fp)
{
	int lineno, ch, np, i;
	char buf[MAXSTRLEN], *s;
	double dummy;

	lineno=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) return 0;

		ungetc(ch, fp);
		++lineno;
		if(!fgets(buf, MAXSTRLEN-1, fp)){ /* read the line found... */
			fprintf(stderr, "countNDoubles(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		/* ...and count the number of doubles it has */
		for(np=i=0, s=buf; 1 ; ++np, s+=i){
			ch=sscanf(s, "%lf%n", &dummy, &i);
			if(ch==0 || ch==EOF) break;
		}

		rewind(fp);
		return np;
	}
	return 0; // should not reach this point
}

int CParallaxBA::skipNDoubles(FILE *fp, int nvals)
{
	register int i;
	int j;

	for(i=0; i<nvals; ++i)
	{
		j=fscanf(fp, "%*f");
		if(j==EOF) return EOF;

		if(ferror(fp)) return EOF-1;
	}

	return nvals;
}


void CParallaxBA::readNpointsAndNprojections(FILE *fp, int *n3Dpts, int pnp, int *nprojs, int mnp)
{
	int nfirst, lineno, npts, nframes, ch, n;

	/* #parameters for the first line */
	nfirst=countNDoubles(fp);

	*n3Dpts=*nprojs=lineno=npts=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n=readNInts(fp, &nframes, 1);
		if(n!=1)
			exit(1);		

		SKIP_LINE(fp);
		*nprojs+=nframes;
		++npts;
	}

	*n3Dpts=npts;
}


void CParallaxBA::pba_readCablibration(FILE* fp, double *K )
{
	int n = fscanf( fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &K[0], &K[3], &K[6], &K[1], &K[4], &K[7], &K[2], &K[5], &K[8] );
	
	if ( n!= 9 )
	{
		fprintf( stderr, "ParallaxBA error: Format of Calibaration is wrong\n" );
		exit(1);
	}
}


void CParallaxBA::pba_readCameraPose( FILE *fp, double *params )
{
	int n, num, lineno = 0 ;
	double *tofilter;
	double * pPrams = params;

	//the number of element per line is 8, it represents that focal length vary, or it is constant
	num = countNDoubles(fp);
	if ( num==8 )
	{	
		m_bFocal = true;
		m_K = (double*)malloc(m_ncams*2*sizeof(double));
		tofilter=(double *)malloc(8*sizeof(double));
	}
	else
		tofilter=(double *)malloc(6*sizeof(double));

	while(!feof(fp))
	{
		if ( num == 6 )
		{
			n	=	readNDoubles(fp, tofilter, 6);
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2]; 
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5]; 
			/*------------------------------------调试用-------------------------------------------*/
			//for (int i = 0;i < 6;i++)
			//	printf("%f ", tofilter[i]);
			//printf("\n");
			/*-------------------------------------调试用------------------------------------------*/
		}
		if ( num == 8 )
		{
			n	=	readNDoubles(fp, tofilter, 8);
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2]; 
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5]; 

			m_K[lineno*2] = tofilter[6];
			m_K[lineno*2+1] = tofilter[7];
		}		

		if( n==-1 ) 
			break;
		pPrams += 6;
		++lineno;
	}
}


int CParallaxBA::readNInts(FILE *fp, int *vals, int nvals)
{
	register int i;
	int n, j;

	for(i=n=0; i<nvals; ++i){//第一次：nvals就是1
		j=fscanf(fp, "%d", vals+i);//第一次：读取第一个，也就是1割三维点的所有2D投影总数
		if(j==EOF) return EOF;//第一次：一般不会是，会是j=1，表示成功

		if(j!=1 || ferror(fp)) return EOF-1;//第一次：一般都是，j=1，表示成功

		n+=j;//第一次因为n=0，所以就是j，也就是成功标志j=1
	}//之后的话，fscanf每读一次那么就会移动指针，然后%d会会跳过空格、换行符和其他空白字符，直到它遇到%d的输入数据
	//或者遇到了EOF或者j!=0了
	return n;
}


int CParallaxBA::readNDoubles(FILE* fp, double* vals, int nvals)
{
	register int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i)//nvals是2，读取xy坐标
	{
		j = fscanf(fp, "%lf", vals + i);//这里开始读的是float，也就是坐标值
		if (j == EOF) return EOF;

		if (j != 1 || ferror(fp)) return EOF - 1;

		n += j;//读完之后加2因为读了2个
	}

	return n;//返回2
}


void CParallaxBA::pba_readProjectionAndInitilizeFeature(FILE *fp,
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
			

			//Add By Zuo
			//const Eigen::Vector2d pt(ptr1[0],ptr1[1]);//像素坐标
			//double* ptr = m_P + frameno * 12;
			//const Eigen::Matrix<double, 3, 4> Pmat=(Eigen::Matrix<double, 3, 4>() <<
			//	ptr[0], ptr[1], ptr[2], ptr[3],
			//	ptr[4], ptr[5], ptr[6], ptr[7],
			//	ptr[8], ptr[9], ptr[10], ptr[11]).finished();
			//A.row(2 * i) = pt(0) * Pmat.row(2) - Pmat.row(0);
			//A.row(2 * i + 1) = pt(1) * Pmat.row(2) - Pmat.row(1);
			//Add By Zuo

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
				///*调试用*/
				//if (bT == true)
				//	printf("%s\n", "true");
				//else
				//	printf("%s\n", "false");
				//printf("%f %f %f\n", (m_motstruct + m_ncams*cnp + ptno * 3)[0], (m_motstruct + m_ncams*cnp + ptno * 3)[1], (m_motstruct + m_ncams*cnp + ptno * 3)[2]);
				///*调试用*/
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

		//Add By Zuo
		//Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
		//Eigen::Vector4d X = svd.matrixV().col(3);
		//X /= X(3);
		//Add By Zuo
		
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
		/*测试umask*/
		//for (int k1 = 0; k1 < m_ncams; k1++)
		//{
		//	for (int k2 = 0; k2 < m_ncams; k2++)
		//	{
		//		printf("%d", m_smask[k1*m_ncams + k2]);
		//	}
		//	printf("\n");
		//}
		//max_nframes = nframes > max_nframes ? nframes : max_nframes;
		/*测试umask*/
		feastart += nframes;
		ptno++;//3D特征点的序号，point NO. OR photo NO.
	}
	/*测试umask*/
	//for (int k1 = 0; k1 < m_ncams; k1++)
	//{
	//	for (int k2 = 0; k2 < m_ncams; k2++)
	//	{
	//		printf("%d", m_smask[k1 * m_ncams + k2]);// -umask[k1 * m_ncams + k2]);
	//	}
	//	printf("\n");
	//}
	//printf("\n\n\n\n\n");
	/*测试umask*/
	//count number of non-zero element in S matrix
	m_nS = 0;
	for ( i = 0; i < m_ncams; i++ ) 
	{
		for (j = 0; j < m_ncams; j++)
		{
			//printf("%d", m_smask[i*m_ncams + j]);
			if (m_smask[i*m_ncams + j] == 1)
			{
				m_nS++;
			}
		}
		//printf("\n");
	}
	//printf("\n\n");
	//for (i = 0; i < m_ncams; i++)
	//{
	//	for (j = 0; j < m_ncams; j++)
	//	{
	//		printf("%d", umask[i*m_ncams + j]);
	//	}
	//	printf("\n");
	//}
}
void CParallaxBA::pba_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams)
{
	int n;
	int nframes;
	int ptno = 0;
	int frameno;
	double* ptr1 = projs;//记录每个投影点的xy而不是序号的（2个2个这样记录和存储）

	int i, j;
	//read all projection point, initialize three feature angle at the same time
	//读取所有投影点，同时初始化三个特征角***重点***
	while (!feof(fp))//一行一行读Match点
	{
		n = readNInts(fp, &nframes, 1);  //读取Match-FeaturePoint文件每一行的第一列，即同名点的个数为nframes/3D点的2D投影个数
		if (n != 1)
			break;//一般都是1，表示成功
		
		Eigen::MatrixXd A(2 * nframes, 4);
		//if (nframes > 3)
		//	nframes = 3;


		for (i = 0; i < nframes; ++i)//一个投影点一个投影点的读取，按总数nframes
		{
			n = readNInts(fp, &frameno, 1); //第二次：读取第一个投影点所在-image index

			if (frameno >= ncams)//一般不会发生
			{
				fprintf(stderr, "ParallaxBA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles(fp, ptr1, 2); //第二次：读取第一个投影点的像点坐标 ptr1[0] ptr1[1]
			
			//也就是1个id+2个坐标
			if (n != 3)//一般不会发生
			{
				fprintf(stderr, "ParallaxBA:reading image projections wrong!\n");
				return;
			}

			const Eigen::Vector2d pt(ptr1[0], ptr1[1]);//像素坐标
			double* ptr = m_P + frameno * 12;
			const Eigen::Matrix<double, 3, 4> Pmat = (Eigen::Matrix<double, 3, 4>() <<
				ptr[0], ptr[1], ptr[2], ptr[3],
				ptr[4], ptr[5], ptr[6], ptr[7],
				ptr[8], ptr[9], ptr[10], ptr[11]).finished();
			A.row(2 * i) = pt(0) * Pmat.row(2) - Pmat.row(0);
			A.row(2 * i + 1) = pt(1) * Pmat.row(2) - Pmat.row(1);
		}

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
		Eigen::Vector4d X = svd.matrixV().col(3);
		X /= X(3);

		(m_motstruct + m_ncams * 6 + ptno * 3)[0] = X[0];
		(m_motstruct + m_ncams * 6 + ptno * 3)[1] = X[1];
		(m_motstruct + m_ncams * 6 + ptno * 3)[2] = X[2];
		ptno++;//3D特征点的序号，point NO. OR photo NO.
	}
}

void CParallaxBA::pba_readAndInitialize(char *camsfname, char *ptsfname, int *ncams,
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

	pba_readCameraPose(fpc, *motstruct);
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
		pba_constructP(m_P, m_K, *motstruct);
	}
	else
		pba_updateKR(m_KR, m_KdA, m_KdB, m_KdG, m_K, *motstruct);//第5个参数为相机内参，第6个参数为相机外参


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
		pba_readProjectionAndTriangulateFeature(fpp, *imgpts, *ncams);
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


void CParallaxBA::pba_readCameraPoseration(char *fname, double* ical )
{
	FILE *fp;
	int  ch=EOF;

	if((fp=fopen(fname, "r"))==NULL)
	{
		fprintf(stderr, "ParallaxBA: Cannot open calbration file %s, exiting\n", fname);
		return;
	}
	//按列读取1-2-3 || 4-5-6 || 7-8-9
	int num = fscanf( fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &ical[0], &ical[3], &ical[6], &ical[1], &ical[4], &ical[7], &ical[2], &ical[5], &ical[8] );	
	if ( num != 9 )
	{
		fprintf(stderr, "ParallaxBA error: Format of Calibration file is wrong");
		return;
	}

	fclose(fp);
}

void CParallaxBA::pba_updateKR( double *KR, double *KdA, double *KdB, double *KdG, double *K, double *p )
{
	if ( !m_bFocal )
	{
		int i = 0;
		double *ptAngle;
		double *pKR, *pKdA, *pKdB, *pKdG;
		double matR[9];//相对旋转矩阵
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];

		for ( i = 0; i < m_ncams; i++ )
		{
			ptAngle = p + i*6;//指向相机外参矩阵
			/*kappa phi omega系统*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0])-cos(ptAngle[2])*sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0])+cos(ptAngle[2])*cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2])*cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0]) + sin(ptAngle[2])*sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0]) - sin(ptAngle[2])*cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2])*cos(ptAngle[1]);

			//omega旋转矩阵
			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);

			//phi旋转矩阵
			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);

			//kappa旋转矩阵
			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			//matRG矩阵关于omega的一阶导
			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[2]);	matDRG[5] = cos(ptAngle[2]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[2]);	matDRG[8] = -sin(ptAngle[2]);

			//matRB矩阵关于phi的一阶导
			matDRB[0] = -sin(ptAngle[1]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[1]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[1]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[1]);

			//matRA矩阵关于kappa的一阶导
			matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//pKR=KR*matR
			pKR = KR + i*9;

			pKR[0] = K[0]*matR[0]+K[3]*matR[3]+K[6]*matR[6];
			pKR[1] = K[0]*matR[1]+K[3]*matR[4]+K[6]*matR[7];
			pKR[2] = K[0]*matR[2]+K[3]*matR[5]+K[6]*matR[8];
			pKR[3] = K[1]*matR[0]+K[4]*matR[3]+K[7]*matR[6];
			pKR[4] = K[1]*matR[1]+K[4]*matR[4]+K[7]*matR[7];
			pKR[5] = K[1]*matR[2]+K[4]*matR[5]+K[7]*matR[8];
			pKR[6] = K[2]*matR[0]+K[5]*matR[3]+K[8]*matR[6];
			pKR[7] = K[2]*matR[1]+K[5]*matR[4]+K[8]*matR[7];
			pKR[8] = K[2]*matR[2]+K[5]*matR[5]+K[8]*matR[8];		

			//pKR矩阵关于omega的一阶导
			pKdG = KdG + i*9;
			tmp1[0] = K[0]*matDRG[0] + K[3]*matDRG[3] + K[6]*matDRG[6];
			tmp1[1] = K[1]*matDRG[0] + K[4]*matDRG[3] + K[7]*matDRG[6];
			tmp1[2] = K[2]*matDRG[0] + K[5]*matDRG[3] + K[8]*matDRG[6];
			tmp1[3] = K[0]*matDRG[1] + K[3]*matDRG[4] + K[6]*matDRG[7];
			tmp1[4] = K[1]*matDRG[1] + K[4]*matDRG[4] + K[7]*matDRG[7];
			tmp1[5] = K[2]*matDRG[1] + K[5]*matDRG[4] + K[8]*matDRG[7];
			tmp1[6] = K[0]*matDRG[2] + K[3]*matDRG[5] + K[6]*matDRG[8];
			tmp1[7] = K[1]*matDRG[2] + K[4]*matDRG[5] + K[7]*matDRG[8];
			tmp1[8] = K[2]*matDRG[2] + K[5]*matDRG[5] + K[8]*matDRG[8];

			tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			pKdG[0] = tmp2[0]*matRA[0]+tmp2[3]*matRA[3]+tmp2[6]*matRA[6];
			pKdG[3] = tmp2[1]*matRA[0]+tmp2[4]*matRA[3]+tmp2[7]*matRA[6];
			pKdG[6] = tmp2[2]*matRA[0]+tmp2[5]*matRA[3]+tmp2[8]*matRA[6];
			pKdG[1] = tmp2[0]*matRA[1]+tmp2[3]*matRA[4]+tmp2[6]*matRA[7];
			pKdG[4] = tmp2[1]*matRA[1]+tmp2[4]*matRA[4]+tmp2[7]*matRA[7];
			pKdG[7] = tmp2[2]*matRA[1]+tmp2[5]*matRA[4]+tmp2[8]*matRA[7];
			pKdG[2] = tmp2[0]*matRA[2]+tmp2[3]*matRA[5]+tmp2[6]*matRA[8];
			pKdG[5] = tmp2[1]*matRA[2]+tmp2[4]*matRA[5]+tmp2[7]*matRA[8];
			pKdG[8] = tmp2[2]*matRA[2]+tmp2[5]*matRA[5]+tmp2[8]*matRA[8];

			//pKR矩阵关于phi的一阶导
			pKdB = KdB + i*9;
			tmp1[0] = K[0]*matRG[0] + K[3]*matRG[3] + K[6]*matRG[6];
			tmp1[1] = K[1]*matRG[0] + K[4]*matRG[3] + K[7]*matRG[6];
			tmp1[2] = K[2]*matRG[0] + K[5]*matRG[3] + K[8]*matRG[6];
			tmp1[3] = K[0]*matRG[1] + K[3]*matRG[4] + K[6]*matRG[7];
			tmp1[4] = K[1]*matRG[1] + K[4]*matRG[4] + K[7]*matRG[7];
			tmp1[5] = K[2]*matRG[1] + K[5]*matRG[4] + K[8]*matRG[7];
			tmp1[6] = K[0]*matRG[2] + K[3]*matRG[5] + K[6]*matRG[8];
			tmp1[7] = K[1]*matRG[2] + K[4]*matRG[5] + K[7]*matRG[8];
			tmp1[8] = K[2]*matRG[2] + K[5]*matRG[5] + K[8]*matRG[8];

			tmp2[0] = tmp1[0]*matDRB[0]+tmp1[3]*matDRB[3]+tmp1[6]*matDRB[6];
			tmp2[1] = tmp1[1]*matDRB[0]+tmp1[4]*matDRB[3]+tmp1[7]*matDRB[6];
			tmp2[2] = tmp1[2]*matDRB[0]+tmp1[5]*matDRB[3]+tmp1[8]*matDRB[6];
			tmp2[3] = tmp1[0]*matDRB[1]+tmp1[3]*matDRB[4]+tmp1[6]*matDRB[7];
			tmp2[4] = tmp1[1]*matDRB[1]+tmp1[4]*matDRB[4]+tmp1[7]*matDRB[7];
			tmp2[5] = tmp1[2]*matDRB[1]+tmp1[5]*matDRB[4]+tmp1[8]*matDRB[7];
			tmp2[6] = tmp1[0]*matDRB[2]+tmp1[3]*matDRB[5]+tmp1[6]*matDRB[8];
			tmp2[7] = tmp1[1]*matDRB[2]+tmp1[4]*matDRB[5]+tmp1[7]*matDRB[8];
			tmp2[8] = tmp1[2]*matDRB[2]+tmp1[5]*matDRB[5]+tmp1[8]*matDRB[8];

			pKdB[0] = tmp2[0]*matRA[0]+tmp2[3]*matRA[3]+tmp2[6]*matRA[6];
			pKdB[3] = tmp2[1]*matRA[0]+tmp2[4]*matRA[3]+tmp2[7]*matRA[6];
			pKdB[6] = tmp2[2]*matRA[0]+tmp2[5]*matRA[3]+tmp2[8]*matRA[6];
			pKdB[1] = tmp2[0]*matRA[1]+tmp2[3]*matRA[4]+tmp2[6]*matRA[7];
			pKdB[4] = tmp2[1]*matRA[1]+tmp2[4]*matRA[4]+tmp2[7]*matRA[7];
			pKdB[7] = tmp2[2]*matRA[1]+tmp2[5]*matRA[4]+tmp2[8]*matRA[7];
			pKdB[2] = tmp2[0]*matRA[2]+tmp2[3]*matRA[5]+tmp2[6]*matRA[8];
			pKdB[5] = tmp2[1]*matRA[2]+tmp2[4]*matRA[5]+tmp2[7]*matRA[8];
			pKdB[8] = tmp2[2]*matRA[2]+tmp2[5]*matRA[5]+tmp2[8]*matRA[8];

			//pKR矩阵关于kappa的一阶导
			pKdA = KdA + i*9;
			tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			pKdA[0] = tmp2[0]*matDRA[0]+tmp2[3]*matDRA[3]+tmp2[6]*matDRA[6];
			pKdA[3] = tmp2[1]*matDRA[0]+tmp2[4]*matDRA[3]+tmp2[7]*matDRA[6];
			pKdA[6] = tmp2[2]*matDRA[0]+tmp2[5]*matDRA[3]+tmp2[8]*matDRA[6];
			pKdA[1] = tmp2[0]*matDRA[1]+tmp2[3]*matDRA[4]+tmp2[6]*matDRA[7];
			pKdA[4] = tmp2[1]*matDRA[1]+tmp2[4]*matDRA[4]+tmp2[7]*matDRA[7];
			pKdA[7] = tmp2[2]*matDRA[1]+tmp2[5]*matDRA[4]+tmp2[8]*matDRA[7];
			pKdA[2] = tmp2[0]*matDRA[2]+tmp2[3]*matDRA[5]+tmp2[6]*matDRA[8];
			pKdA[5] = tmp2[1]*matDRA[2]+tmp2[4]*matDRA[5]+tmp2[7]*matDRA[8];
			pKdA[8] = tmp2[2]*matDRA[2]+tmp2[5]*matDRA[5]+tmp2[8]*matDRA[8];

		}
	}
	else
	{
		int i = 0;
		double *ptAngle;
		double *pKR, *pKdA, *pKdB, *pKdG;
		double matR[9];
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];
		double K[9];
		memset( K, 0, 9*sizeof(double));
		K[8] = 1;
		for ( i = 0; i < m_ncams; i++ )
		{
			ptAngle = p + i*6;

			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0])-cos(ptAngle[2])*sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0])+cos(ptAngle[2])*cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2])*cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0]) + sin(ptAngle[2])*sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0]) - sin(ptAngle[2])*cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2])*cos(ptAngle[1]);

			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);

			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);

			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[2]);	matDRG[5] = cos(ptAngle[2]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[2]);	matDRG[8] = -sin(ptAngle[2]);

			matDRB[0] = -sin(ptAngle[1]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[1]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[1]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[1]);

			matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//KR

			K[0] = m_K[i*2];
			K[4] = m_K[i*2+1];

			pKR = KR + i*9;
			pKR[0] = K[0]*matR[0]+K[3]*matR[3]+K[6]*matR[6];
			pKR[1] = K[0]*matR[1]+K[3]*matR[4]+K[6]*matR[7];
			pKR[2] = K[0]*matR[2]+K[3]*matR[5]+K[6]*matR[8];
			pKR[3] = K[1]*matR[0]+K[4]*matR[3]+K[7]*matR[6];
			pKR[4] = K[1]*matR[1]+K[4]*matR[4]+K[7]*matR[7];
			pKR[5] = K[1]*matR[2]+K[4]*matR[5]+K[7]*matR[8];
			pKR[6] = K[2]*matR[0]+K[5]*matR[3]+K[8]*matR[6];
			pKR[7] = K[2]*matR[1]+K[5]*matR[4]+K[8]*matR[7];
			pKR[8] = K[2]*matR[2]+K[5]*matR[5]+K[8]*matR[8];		

			//KdG
			pKdG = KdG + i*9;
			tmp1[0] = K[0]*matDRG[0] + K[3]*matDRG[3] + K[6]*matDRG[6];
			tmp1[1] = K[1]*matDRG[0] + K[4]*matDRG[3] + K[7]*matDRG[6];
			tmp1[2] = K[2]*matDRG[0] + K[5]*matDRG[3] + K[8]*matDRG[6];
			tmp1[3] = K[0]*matDRG[1] + K[3]*matDRG[4] + K[6]*matDRG[7];
			tmp1[4] = K[1]*matDRG[1] + K[4]*matDRG[4] + K[7]*matDRG[7];
			tmp1[5] = K[2]*matDRG[1] + K[5]*matDRG[4] + K[8]*matDRG[7];
			tmp1[6] = K[0]*matDRG[2] + K[3]*matDRG[5] + K[6]*matDRG[8];
			tmp1[7] = K[1]*matDRG[2] + K[4]*matDRG[5] + K[7]*matDRG[8];
			tmp1[8] = K[2]*matDRG[2] + K[5]*matDRG[5] + K[8]*matDRG[8];

			tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			pKdG[0] = tmp2[0]*matRA[0]+tmp2[3]*matRA[3]+tmp2[6]*matRA[6];
			pKdG[3] = tmp2[1]*matRA[0]+tmp2[4]*matRA[3]+tmp2[7]*matRA[6];
			pKdG[6] = tmp2[2]*matRA[0]+tmp2[5]*matRA[3]+tmp2[8]*matRA[6];
			pKdG[1] = tmp2[0]*matRA[1]+tmp2[3]*matRA[4]+tmp2[6]*matRA[7];
			pKdG[4] = tmp2[1]*matRA[1]+tmp2[4]*matRA[4]+tmp2[7]*matRA[7];
			pKdG[7] = tmp2[2]*matRA[1]+tmp2[5]*matRA[4]+tmp2[8]*matRA[7];
			pKdG[2] = tmp2[0]*matRA[2]+tmp2[3]*matRA[5]+tmp2[6]*matRA[8];
			pKdG[5] = tmp2[1]*matRA[2]+tmp2[4]*matRA[5]+tmp2[7]*matRA[8];
			pKdG[8] = tmp2[2]*matRA[2]+tmp2[5]*matRA[5]+tmp2[8]*matRA[8];

			//KdB
			pKdB = KdB + i*9;
			tmp1[0] = K[0]*matRG[0] + K[3]*matRG[3] + K[6]*matRG[6];
			tmp1[1] = K[1]*matRG[0] + K[4]*matRG[3] + K[7]*matRG[6];
			tmp1[2] = K[2]*matRG[0] + K[5]*matRG[3] + K[8]*matRG[6];
			tmp1[3] = K[0]*matRG[1] + K[3]*matRG[4] + K[6]*matRG[7];
			tmp1[4] = K[1]*matRG[1] + K[4]*matRG[4] + K[7]*matRG[7];
			tmp1[5] = K[2]*matRG[1] + K[5]*matRG[4] + K[8]*matRG[7];
			tmp1[6] = K[0]*matRG[2] + K[3]*matRG[5] + K[6]*matRG[8];
			tmp1[7] = K[1]*matRG[2] + K[4]*matRG[5] + K[7]*matRG[8];
			tmp1[8] = K[2]*matRG[2] + K[5]*matRG[5] + K[8]*matRG[8];

			tmp2[0] = tmp1[0]*matDRB[0]+tmp1[3]*matDRB[3]+tmp1[6]*matDRB[6];
			tmp2[1] = tmp1[1]*matDRB[0]+tmp1[4]*matDRB[3]+tmp1[7]*matDRB[6];
			tmp2[2] = tmp1[2]*matDRB[0]+tmp1[5]*matDRB[3]+tmp1[8]*matDRB[6];
			tmp2[3] = tmp1[0]*matDRB[1]+tmp1[3]*matDRB[4]+tmp1[6]*matDRB[7];
			tmp2[4] = tmp1[1]*matDRB[1]+tmp1[4]*matDRB[4]+tmp1[7]*matDRB[7];
			tmp2[5] = tmp1[2]*matDRB[1]+tmp1[5]*matDRB[4]+tmp1[8]*matDRB[7];
			tmp2[6] = tmp1[0]*matDRB[2]+tmp1[3]*matDRB[5]+tmp1[6]*matDRB[8];
			tmp2[7] = tmp1[1]*matDRB[2]+tmp1[4]*matDRB[5]+tmp1[7]*matDRB[8];
			tmp2[8] = tmp1[2]*matDRB[2]+tmp1[5]*matDRB[5]+tmp1[8]*matDRB[8];

			pKdB[0] = tmp2[0]*matRA[0]+tmp2[3]*matRA[3]+tmp2[6]*matRA[6];
			pKdB[3] = tmp2[1]*matRA[0]+tmp2[4]*matRA[3]+tmp2[7]*matRA[6];
			pKdB[6] = tmp2[2]*matRA[0]+tmp2[5]*matRA[3]+tmp2[8]*matRA[6];
			pKdB[1] = tmp2[0]*matRA[1]+tmp2[3]*matRA[4]+tmp2[6]*matRA[7];
			pKdB[4] = tmp2[1]*matRA[1]+tmp2[4]*matRA[4]+tmp2[7]*matRA[7];
			pKdB[7] = tmp2[2]*matRA[1]+tmp2[5]*matRA[4]+tmp2[8]*matRA[7];
			pKdB[2] = tmp2[0]*matRA[2]+tmp2[3]*matRA[5]+tmp2[6]*matRA[8];
			pKdB[5] = tmp2[1]*matRA[2]+tmp2[4]*matRA[5]+tmp2[7]*matRA[8];
			pKdB[8] = tmp2[2]*matRA[2]+tmp2[5]*matRA[5]+tmp2[8]*matRA[8];

			//KdA
			pKdA = KdA + i*9;
			tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			pKdA[0] = tmp2[0]*matDRA[0]+tmp2[3]*matDRA[3]+tmp2[6]*matDRA[6];
			pKdA[3] = tmp2[1]*matDRA[0]+tmp2[4]*matDRA[3]+tmp2[7]*matDRA[6];
			pKdA[6] = tmp2[2]*matDRA[0]+tmp2[5]*matDRA[3]+tmp2[8]*matDRA[6];
			pKdA[1] = tmp2[0]*matDRA[1]+tmp2[3]*matDRA[4]+tmp2[6]*matDRA[7];
			pKdA[4] = tmp2[1]*matDRA[1]+tmp2[4]*matDRA[4]+tmp2[7]*matDRA[7];
			pKdA[7] = tmp2[2]*matDRA[1]+tmp2[5]*matDRA[4]+tmp2[8]*matDRA[7];
			pKdA[2] = tmp2[0]*matDRA[2]+tmp2[3]*matDRA[5]+tmp2[6]*matDRA[8];
			pKdA[5] = tmp2[1]*matDRA[2]+tmp2[4]*matDRA[5]+tmp2[7]*matDRA[8];
			pKdA[8] = tmp2[2]*matDRA[2]+tmp2[5]*matDRA[5]+tmp2[8]*matDRA[8];

		}
	}	
}

//By Zuo
void CParallaxBA::pba_constructP(double* P, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9];
		double matT[3];

		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;//指向相机外参矩阵
			/*kappa phi omega系统*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];

			//test
			//double t1 = pP[0] * (2 - ptAngle[3]) + pP[1] * (2 - ptAngle[4]) + pP[2] * (2 - ptAngle[5]);
			//double t2 = pP[4] * (2 - ptAngle[3]) + pP[5] * (2 - ptAngle[4]) + pP[6] * (2 - ptAngle[5]);
			//double t3 = pP[8] * (2 - ptAngle[3]) + pP[9] * (2 - ptAngle[4]) + pP[10] * (2 - ptAngle[5]);

			//double u1 = t1 / t3;
			//double v1 = t2 / t3;

			//t1 = pP[0] * 2 + pP[1] * 2 + pP[2] * 2 + pP[3];
			//t2 = pP[4] * 2 + pP[5] * 2 + pP[6] * 2 + pP[7];
			//t3 = pP[8] * 2 + pP[9] * 2 + pP[10] * 2 + pP[11];
			//double u2 = t1 / t3;
			//double v2 = t2 / t3;
		}
	}
	else
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9],matT[3];
		double K[9];
		memset(K, 0, 9 * sizeof(double));
		K[8] = 1;
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;

			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//KR

			K[0] = m_K[i * 2];
			K[4] = m_K[i * 2 + 1];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];
		}
	}
}

//--------------------------------------------------------
void CParallaxBA::pba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
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

void CParallaxBA::pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
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
/*void CParallaxBA::pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
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

double CParallaxBA::nrmL2xmy(double *const e, const double *const x, const double *const y, const int n)
{
	const int blocksize=8, bpwr=3; /* 8=2^3 */
	register int i;
	int j1, j2, j3, j4, j5, j6, j7;
	int blockn;
	register double sum0=0.0, sum1=0.0, sum2=0.0, sum3=0.0;

	/* n may not be divisible by blocksize, 
	* go as near as we can first, then tidy up.
	*/
	//printf("nobs%d\n", n);//Add by Lu
	blockn = (n>>bpwr)<<bpwr; /* (n / blocksize) * blocksize; */
	//printf("cycle_num%d, blocksize%d", blockn - 1, blocksize);//Add by Lu

	/* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
	for(i=blockn-1; i>0; i-=blocksize){
		e[i ]=x[i ]-y[i ]; sum0+=e[i ]*e[i ];
		j1=i-1; e[j1]=x[j1]-y[j1]; sum1+=e[j1]*e[j1];
		j2=i-2; e[j2]=x[j2]-y[j2]; sum2+=e[j2]*e[j2];
		j3=i-3; e[j3]=x[j3]-y[j3]; sum3+=e[j3]*e[j3];
		j4=i-4; e[j4]=x[j4]-y[j4]; sum0+=e[j4]*e[j4];
		j5=i-5; e[j5]=x[j5]-y[j5]; sum1+=e[j5]*e[j5];
		j6=i-6; e[j6]=x[j6]-y[j6]; sum2+=e[j6]*e[j6];
		j7=i-7; e[j7]=x[j7]-y[j7]; sum3+=e[j7]*e[j7];
	}
	
	i=blockn;
	if(i<n){ 
		switch(n - i){ 
		case 7 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 6 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 5 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 4 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 3 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 2 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 1 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		}
	}

	return sum0+sum1+sum2+sum3;
}

void CParallaxBA::pba_cost(double* p, double* hx, int* archor)
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

void CParallaxBA::readNpointsAndNprojectionsFromProj(FILE *fp, int &n3Dpts, int &nprojs)
{
	int nfirst, lineno, npts, nframes, ch, n;
	nprojs = 0;
	n3Dpts = 0;
	npts = 0;

	/* #parameters for the first line */
	nfirst=countNDoubles(fp);

	//*n3Dpts=*nprojs=lineno=npts=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n=readNInts(fp, &nframes, 1);
		if(n!=1)
		{
			fprintf(stderr, "readNpointsAndNprojections(): error reading input file, line %d: "
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		SKIP_LINE(fp);
		nprojs+=nframes;
		++npts;
	}

	n3Dpts=npts;
}

void CParallaxBA::readPointProjections(FILE *fp,double *imgpts, int *photo,int* imgptsSum, int n3Dpts, int n2Dprojs )
{
	int nframes, ch, lineno, ptno, frameno, n;
	int i;
	int nproj2D = 0;

	lineno=ptno=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			lineno++;

			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);

		n=readNInts(fp, &nframes, 1);  /* read in number of image projections */
		if(n!=1)
		{
			fprintf(stderr, "pba_readProjectionAndInitilizeFeature(): error reading input file, line %d:\n"
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		imgptsSum[ptno] = nframes;

		for(i=0; i<nframes; ++i)
		{
			n=readNInts(fp, &frameno, 1); /* read in frame number... */

			photo[nproj2D] = frameno;

			n+=readNDoubles(fp, imgpts+nproj2D*2, 2); /* ...and image projection */

			nproj2D++;
		}
		fscanf(fp, "\n"); // consume trailing newline

		lineno++;
		ptno++;
	}
}

void CParallaxBA::readImagePts( const char* szProj, double **imgpts, int **photo,int** imgptsSum, int &n3Dpts, int &n2Dprojs )
{
	FILE  *fpp;
	if((fpp=fopen(szProj, "r"))==NULL){
		fprintf(stderr, "cannot open file %s, exiting\n", szProj);
		exit(1);
	}
	readNpointsAndNprojectionsFromProj(fpp, n3Dpts, n2Dprojs);

	*imgpts = (double*)malloc(n2Dprojs*2*sizeof(double));
	if(*imgpts==NULL){
		fprintf(stderr, "memory allocation for 'imgpts' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*photo = (int*)malloc(n2Dprojs*sizeof(int)); 
	if (*photo==NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*imgptsSum = (int*)malloc(n3Dpts*sizeof(int));
	if (*imgptsSum==NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	rewind(fpp);
	readPointProjections(fpp, *imgpts,*photo,*imgptsSum, n3Dpts, n2Dprojs );

	fclose(fpp);
}

/*初始化主锚点*/
bool CParallaxBA::pba_initializeMainArchor( 
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
bool CParallaxBA::pba_initializeAssoArchor( 
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
bool CParallaxBA::pba_initializeOtheArchors( 
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


bool CParallaxBA::pba_initializeOtheArchors_Mindw(
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
int CParallaxBA::pba_angle2xytGN( double *p )
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

void CParallaxBA::pba_saveXYZ( const char* camera, const char* sz3Dpt, double *p,bool gn  )
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

int	CParallaxBA::pba_angle2xytLM( double *p )
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

bool CParallaxBA::pba_parseArgs( int argc, char* argv[] )
{
	int i;
	string param;
	bool bSuccess, bRKF;

	for (i = 1; i < argc; i++) 
	{
		bSuccess = false;
		string name = argv[i];

		if (name[0] != '-') { // each param has to start with at least one dash
			return false;
		}
		
		string::size_type dashPos = name.find_first_not_of('-');
		if (dashPos != string::npos)
			name = name.substr(dashPos);

		if ( strcmp(name.c_str(), "help") == 0 )
		{
			pba_printHelp();
			return false;
		}
		
		if ( strcmp( name.c_str(), "cam") == 0 )
		{
			i++;
			param = argv[i];
			m_szCameraInit = (char*)malloc(param.length());
			strcpy( m_szCameraInit, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "fea") == 0 )
		{
			i++;
			param = argv[i];
			m_szFeatures = (char*)malloc(param.length());
			strcpy( m_szFeatures, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "calib") == 0 )
		{
			i++;
			param = argv[i];
			m_szCalibration = (char*)malloc(param.length());
			strcpy( m_szCalibration, param.c_str() );			
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "pose") == 0 )
		{
			i++;
			param = argv[i];
			m_szCamePose = (char*)malloc(param.length());
			strcpy( m_szCamePose, param.c_str() );			
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "3D") == 0 )
		{
			i++;
			param = argv[i];
			m_sz3Dpts = (char*)malloc(param.length());
			strcpy( m_sz3Dpts, param.c_str() );			
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "report") == 0 )
		{
			i++;
			param = argv[i];
			m_szReport = (char*)malloc(param.length());
			strcpy( m_szReport, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "xyz") == 0 )
		{
			i++;
			param = argv[i];
			m_szXYZ = (char*)malloc(param.length());
			strcpy( m_szXYZ, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "i") == 0 )
		{
			i++;
			param = argv[i];
			m_nMaxIter = atoi(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "robustKernel") == 0 )
		{
			m_bRobustKernel = true;
			bRKF = false;
			i++;
			param = argv[i];
			if ( strcmp(param.c_str(), "Huber") == 0 )
			{
				m_nRobustType = 2;
				bRKF = true;
			}

			if ( strcmp(param.c_str(), "Cauchy") == 0 )
			{
				m_nRobustType = 1;
				bRKF = true;
			}
			
			bSuccess = true;

			if ( !bRKF )
			{
				printf( "ParallaxBA: Must input right robust kernel function!\n" );
				return false;
			}
		}

		if ( strcmp( name.c_str(), "solve") == 0 )
		{
			i++;
			param = argv[i];
			if( strcmp(param.c_str(), "LM") == 0)
				m_bsolverLM = true;

			if( strcmp(param.c_str(), "GN") == 0)
				m_bsolverGN = true;

			bSuccess = true;
		}		

		if ( strcmp( name.c_str(), "t") == 0 )
		{
			i++;
			param = argv[i];
			
			m_Tau = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e1") == 0 )
		{
			i++;
			param = argv[i];

			m_e1 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e2") == 0 )
		{
			i++;
			param = argv[i];

			m_e2 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e3") == 0 )
		{
			i++;
			param = argv[i];

			m_e3 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e4") == 0 )
		{
			i++;
			param = argv[i];

			m_e4 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "robustKernelWidth") == 0 )
		{
			i++;
			param = argv[i];

			m_delt = atof(param.c_str());
			bSuccess = true;
		}

		if (!bSuccess)
		{
			printf( "ParallaxBA error: %s command is wrong!\n", name.c_str() );
			return false;
		}
		
	}

	return true;
			
}

void CParallaxBA::pba_printHelp()
{
	printf( "Parallax Bundle Adjustment General Options\n" );
	printf( "\n" );

	printf( "-cam			Provide initial camera pose.\n" );
	printf( "-fea			Provide features.\n" );
	printf( "-calib			Provide calibration.\n" );
	printf( "-xyz			Provide initial XYZ.\n" );
	printf( "-pose			Output optimal camera pose.\n" );
	printf( "-3D			Output optimal 3D point cloud.\n" );
	printf( "-report			Output report.\n" );
	printf( "-solve			Solve method including LevenbergMarquart(LM) and Gauss-Newton(GN).\n" );
	printf( "-i			Max Iteration.\n" );
	printf( "-t			LevenbergMarquart parameters.\n" );
	printf( "-robustKernel		use Cauchy Robust Kernel Function.\n" );
	printf( "-robustKernelWidth		width for the robust Kernel.\n" );	
}
