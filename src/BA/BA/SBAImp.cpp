
#include "stdafx.h"
#include "SBAImp.h"

//compute reprojection image coordinates of each image points
static void sba_reprojectEachPts( double *KR, double* pa, double* pb, int nP, double n[2] )
{
	double ptXj[3];
	double *pKR;
	
	pKR = KR + nP * 9;//nP�����ͼ��KR����

	ptXj[0] = pb[0] - pa[0];
	ptXj[1] = pb[1] - pa[1];
	ptXj[2] = pb[2] - pa[2];
	n[0] = (pKR[0] * ptXj[0] + pKR[1] * ptXj[1] + pKR[2] * ptXj[2]) /
		(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

	n[1] = (pKR[3] * ptXj[0] + pKR[4] * ptXj[1] + pKR[5] * ptXj[2]) /
		(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);
}
//compute Jacobian for each image point
//pPA��pPB�ֱ���uv����������κ���ά���һ�׵�
static void sba_jacobianEachPts(double* KR, double* KdA, double *KdB, double* KdG, double* pa, double* ppt, int nP, double* pPA, double* pPB)
{
	double matXj[3];
	double matxyt[3];
	double matDuvDxyt[6];
	double matDxytDRA[3], matDxytDRB[3], matDxytDRG[3], matDxytDXc[3], matDxytDYc[3], matDxytDZc[3];
	double matDuvDRA[2], matDuvDRB[2], matDuvDRG[2], matDuvDXc[2], matDuvDYc[2], matDuvDZc[2], matDuvDX[2], matDuvDY[2], matDuvDZ[2];
	double *pKR, *pKdA, *pKdB, *pKdG;//pKdA��pKdB��pKdG�ֱ���pKR�������kappa��phi��omega��һ�׵�

	double* ppa = pa + nP * 6 + 3;//��ǰ��ͼλ��ָ��
	matXj[0] = ppt[0] - ppa[0];
	matXj[1] = ppt[1] - ppa[1];
	matXj[2] = ppt[2] - ppa[2];
	pKR = KR + nP * 9;
	matxyt[0] = pKR[0] * matXj[0] + pKR[1] * matXj[1] + pKR[2] * matXj[2];
	matxyt[1] = pKR[3] * matXj[0] + pKR[4] * matXj[1] + pKR[5] * matXj[2];
	matxyt[2] = pKR[6] * matXj[0] + pKR[7] * matXj[1] + pKR[8] * matXj[2];

	matDuvDxyt[0] = 1 / matxyt[2];
	matDuvDxyt[1] = 0;
	matDuvDxyt[2] = -matxyt[0] / (matxyt[2] * matxyt[2]);
	matDuvDxyt[3] = 0;
	matDuvDxyt[4] = 1 / matxyt[2];
	matDuvDxyt[5] = -matxyt[1] / (matxyt[2] * matxyt[2]);

	matDxytDXc[0] = -pKR[0];
	matDxytDXc[1] = -pKR[3];
	matDxytDXc[2] = -pKR[6];

	matDxytDYc[0] = -pKR[1];
	matDxytDYc[1] = -pKR[4];
	matDxytDYc[2] = -pKR[7];

	matDxytDZc[0] = -pKR[2];
	matDxytDZc[1] = -pKR[5];
	matDxytDZc[2] = -pKR[8];

	matDuvDXc[0] = matDuvDxyt[0] * matDxytDXc[0] + matDuvDxyt[1] * matDxytDXc[1] + matDuvDxyt[2] * matDxytDXc[2];//uv����Xc��һ�׵�
	matDuvDXc[1] = matDuvDxyt[3] * matDxytDXc[0] + matDuvDxyt[4] * matDxytDXc[1] + matDuvDxyt[5] * matDxytDXc[2];

	matDuvDYc[0] = matDuvDxyt[0] * matDxytDYc[0] + matDuvDxyt[1] * matDxytDYc[1] + matDuvDxyt[2] * matDxytDYc[2];//uv����Yc��һ�׵�
	matDuvDYc[1] = matDuvDxyt[3] * matDxytDYc[0] + matDuvDxyt[4] * matDxytDYc[1] + matDuvDxyt[5] * matDxytDYc[2];

	matDuvDZc[0] = matDuvDxyt[0] * matDxytDZc[0] + matDuvDxyt[1] * matDxytDZc[1] + matDuvDxyt[2] * matDxytDZc[2];//uv����Zc��һ�׵�
	matDuvDZc[1] = matDuvDxyt[3] * matDxytDZc[0] + matDuvDxyt[4] * matDxytDZc[1] + matDuvDxyt[5] * matDxytDZc[2];

	matDuvDX[0] = -matDuvDXc[0];
	matDuvDX[1] = -matDuvDXc[1];

	matDuvDY[0] = -matDuvDYc[0];
	matDuvDY[1] = -matDuvDYc[1];

	matDuvDZ[0] = -matDuvDZc[0];
	matDuvDZ[1] = -matDuvDZc[1];

	//camera angles
	pKdG = KdG + nP * 9;
	matDxytDRG[0] = pKdG[0] * matXj[0] + pKdG[1] * matXj[1] + pKdG[2] * matXj[2];//xyt����omega��һ�׵�
	matDxytDRG[1] = pKdG[3] * matXj[0] + pKdG[4] * matXj[1] + pKdG[5] * matXj[2];
	matDxytDRG[2] = pKdG[6] * matXj[0] + pKdG[7] * matXj[1] + pKdG[8] * matXj[2];

	pKdB = KdB + nP * 9;
	matDxytDRB[0] = pKdB[0] * matXj[0] + pKdB[1] * matXj[1] + pKdB[2] * matXj[2];//xyt����phi��һ�׵�
	matDxytDRB[1] = pKdB[3] * matXj[0] + pKdB[4] * matXj[1] + pKdB[5] * matXj[2];
	matDxytDRB[2] = pKdB[6] * matXj[0] + pKdB[7] * matXj[1] + pKdB[8] * matXj[2];

	pKdA = KdA + nP * 9;
	matDxytDRA[0] = pKdA[0] * matXj[0] + pKdA[1] * matXj[1] + pKdA[2] * matXj[2];//xyt����kappa��һ�׵�
	matDxytDRA[1] = pKdA[3] * matXj[0] + pKdA[4] * matXj[1] + pKdA[5] * matXj[2];
	matDxytDRA[2] = pKdA[6] * matXj[0] + pKdA[7] * matXj[1] + pKdA[8] * matXj[2];

	matDuvDRA[0] = matDuvDxyt[0] * matDxytDRA[0] + matDuvDxyt[1] * matDxytDRA[1] + matDuvDxyt[2] * matDxytDRA[2];//uv����kappa��һ�׵�
	matDuvDRA[1] = matDuvDxyt[3] * matDxytDRA[0] + matDuvDxyt[4] * matDxytDRA[1] + matDuvDxyt[5] * matDxytDRA[2];

	matDuvDRB[0] = matDuvDxyt[0] * matDxytDRB[0] + matDuvDxyt[1] * matDxytDRB[1] + matDuvDxyt[2] * matDxytDRB[2];//uv����phi��һ�׵�
	matDuvDRB[1] = matDuvDxyt[3] * matDxytDRB[0] + matDuvDxyt[4] * matDxytDRB[1] + matDuvDxyt[5] * matDxytDRB[2];

	matDuvDRG[0] = matDuvDxyt[0] * matDxytDRG[0] + matDuvDxyt[1] * matDxytDRG[1] + matDuvDxyt[2] * matDxytDRG[2];//uv����omega��һ�׵�
	matDuvDRG[1] = matDuvDxyt[3] * matDxytDRG[0] + matDuvDxyt[4] * matDxytDRG[1] + matDuvDxyt[5] * matDxytDRG[2];

	pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv����ε�һ�׵�
	pPA[3] = matDuvDXc[0];			pPA[4] = matDuvDYc[0];			pPA[5] = matDuvDZc[0];
	pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
	pPA[9] = matDuvDXc[1];			pPA[10] = matDuvDYc[1];			pPA[11] = matDuvDZc[1];

	pPB[0] = matDuvDX[0];		pPB[1] = matDuvDY[0];		pPB[2] = matDuvDZ[0];//uv����ά���һ�׵�
	pPB[3] = matDuvDX[1];		pPB[4] = matDuvDY[1];		pPB[5] = matDuvDZ[1];
}

SBA::SBA(void)
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


SBA::~SBA(void)
{
	if ( m_szCameraInit!= NULL)	free(m_szCameraInit);
	if ( m_szFeatures!= NULL)	free(m_szFeatures);
	if ( m_szCalibration!= NULL)	free(m_szCalibration);
	if ( m_szXYZ!= NULL)		free(m_szXYZ);
	if ( m_szCamePose!= NULL)	free(m_szCamePose);
	if ( m_sz3Dpts!= NULL)		free(m_sz3Dpts);
	if ( m_szReport!= NULL)		free(m_szReport);
}

bool SBA::ba_run( bool bRobust,
	bool bLM,
	int nMaxIter,
	char* szCam,
	char* szFea,
	char* szXYZ,
	char* szCalib,
	char* szReport,
	char* szPose,
	char* sz3D,
	double Tau )
{
	m_szCameraInit = szCam;                                                       //��ʼ�������ļ�·��
	m_szFeatures   = szFea;                                                       //Track�ļ�·��
	m_szCalibration= szCalib;                                                     //����ڲ��ļ�·��
	m_szXYZ        = szXYZ;                                                       //��ʼ��ά���ļ�·��
	m_bRobustKernel= bRobust;                                                     
	m_bsolverGN    = !bLM;                                                        
	m_nMaxIter     = nMaxIter;                                                    
	m_szCamePose   = szPose;                                                      //��������ļ�����·��
	m_sz3Dpts      = sz3D;                                                        //������ά���ļ�����·��
	m_szReport     = szReport;                                                    
	m_Tau          = Tau;                                                         

	BAType ba = sba;
	
	ba_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );     

	if (m_bsolverGN)
		ba_motstr_gn( );
	else
		ba_motstr_levmar();
	return true;
}

bool SBA::ba_run( int argc, char** argv )
{
	bool bTrue = ba_parseArgs( argc, argv );
	if (!bTrue)
	{
		fprintf( stderr, "SparseBA: Input wrong commands, please check them!\n");
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


bool SBA::ba_initialize( char* szCamera, char* szFeature,  char* szCalib, char* szXYZ )
{
	printf("SparseBA: Sparse Bundle Adjustment Version 1.0\n");
	FILE* fp;

	//must input initial initial camera pose file and projection image points file
	fp = fopen( szCamera, "r" );
	if ( fp == NULL )
	{
		fprintf( stderr, "SparseBA: Missing initial camera poses file! \n");
		exit(1);
	}
	else
		fclose(fp);

	fp = fopen( szFeature, "r" );
	if ( fp == NULL )
	{	
		fprintf( stderr, "SparseBA: Missing feature projection points file! \n");
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

	savePara = 0;
	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	sba_readAndInitialize( szCamera, szFeature, &m_ncams, &m_n3Dpts, &m_n2Dprojs,&m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
				 &m_imgpts, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort );
	
	std::string originalPath(szCamera);
	size_t pos = originalPath.find_last_of("/\\");
	std::string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	std::string newPath = parentPath + "/" + "Triangulatedxyz.txt";

	if (savePara == 1)
		ba_saveTriangulatedxyz(newPath.c_str(), m_motstruct);

	printf( "Number of cameras: %d\n", m_ncams );
	printf( "Number of points: %d\n", m_n3Dpts );
	printf( "Number of projections: %d\n", m_n2Dprojs );
	return true;
}

bool SBA::ba_motstr_levmar( )
{
	if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
	{	
		fprintf( stderr, "SparseBA: Missing related source file, please input them first\n");	
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
	esz=mnp; easz=cnp; ebsz=pnp;								//2 �� 6 �� 3
	Sblsz=cnp*cnp;												//6*6
	Sdim=m * cnp;												//�����Ŀ*6
	nvis = m_n2Dprojs;											//�����
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
	sba_cost(p, hx, m_archor ); 

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
			sba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			sba_jacobian(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
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
				
				sba_cost(pdp, hx, m_archor );
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
	sba_saveXYZ( m_szCamePose, m_sz3Dpts, p);

	free(S);	free(W);	free(U);	free(V);	free(IV);
	free(e);	free(eab);	free(E);   	free(dp);	free(hx);	free(pdp);	
	free(m_KR); free(m_KdA);free(m_KdB);free(m_KdG);	

	printf( "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
		nvars, m_n2Dprojs*2, stop, p_eL2/nvis, initialerror, nIter+1, nLinear, tTimeUse / CLOCKS_PER_SEC );
	printf( "SparseBA reasons is listed as following:\n" );
	printf( "reason 1: relative change of state vector is small\n" );
	printf( "reason 2: relative change of projection error is small\n" );
	printf( "reason 3: maximum iteration\n");
	printf( "reason 4: maximum value of b is small\n");
	printf( "reason 5: total reprojection error is small\n");

	if( m_szReport!=NULL )
	{
		fprintf( fpRe, "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, error, initialerror, nIter+1, nLinear, tTimeUse*0.001 );

		fprintf( fpRe, "SparseBA reasons is listed as following:\n" );
		fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
		fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
		fprintf( fpRe, "reason 3: maximum iteration\n");
		fprintf( fpRe, "reason 4: maximum value of b is small\n");
		fprintf( fpRe, "reason 5: total reprojection error is small\n");

		fclose(fpRe);
	}	
	
	return true;
}

bool SBA::ba_motstr_gn( FixType ft )
{	
	if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
	{	
		printf( "SparseBA: Missing related file, please input them first\n");	
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
	struct sba_crsm Sidxij, Uidxij;		// S mask and U mask          �洢ϡ�����
	nuis = ba_ConstructSmask( Sidxij, Uidxij );   //����U�������Ԫ�ظ���

	Usz=cnp*cnp; Vsz=pnp * pnp; Wsz=cnp * pnp; 	
	esz=mnp; easz=cnp; ebsz=pnp; Sblsz=cnp*cnp;	Sdim=m * cnp;	
	nvis = m_n2Dprojs;                                            //�������
	nobs=nvis*mnp;                                                //�۲ⷽ������
	nvars=m*cnp + n*pnp;                                          //δ֪��������

	S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));          //S���󣺷�0Ԫ�ظ���*6*6
	W	=	(double *)malloc(nvis*Wsz*sizeof(double));            //W���������*6*3
	U	=	(double *)malloc(nuis*Usz*sizeof(double));            //U���󣺷�0Ԫ�ظ�����6*6
	V	=	(double *)malloc(n*Vsz*sizeof(double));               //V������ά������3*3
	e	=	(double *)malloc(nobs*sizeof(double));                //����� 
	eab	=	(double *)malloc(nvars*sizeof(double));               //δ֪����
	E	=	(double *)malloc(m*cnp*sizeof(double));	              //�ⷽλԪ��	
	dp	=	(double *)malloc(nvars*sizeof(double));               //δ֪����
	hx	=	(double *)malloc(nobs*sizeof(double));	              //��ͶӰuv����
	pa=p; //paָ��δ֪����
	pb=p+m*cnp; //pbָ����ά�����

	ea=eab;//eaָ�������
	eb=eab+m*cnp;//ebָ����ά�������

	dpa=dp;//dpaָ��δ֪���������� 
	dpb=dp+m*cnp;//dpbָ����ά�������
	
	//Select fix axis for Gauss-Newton
	//����BAֻ��������͵�����λ�ã�������̶�ĳЩ�������Ż������п��ܻᷢ��Drift��������������������ƽ�ơ���ת�����ţ�������ı�ͶӰ�����½ⲻΨһ��
	//�̶���һ֡�͵ڶ�֡������һ��������
	if ( ft == BA_FixDefault )
	{
		tmp = abs(*(m_motstruct+9)-*(m_motstruct+3));//��2֡Xc-��1֡Xc
		ft = BA_FixX;									
		tmp1 = abs(*(m_motstruct+10)-*(m_motstruct+4));//��2֡Yc-��1֡Yc
		if ( tmp1 > tmp )
		{ft= BA_FixY; tmp = tmp1;	}
		tmp1 = abs(*(m_motstruct+11)-*(m_motstruct+5));//��2֡Zc-��1֡Zc
		if ( tmp1 > tmp )
		{ft= BA_FixZ; }
	}
	nft = ft;
	
	// ʹ�� CHOLMOD ����ϡ��������ӷֽ⣨Cholesky �ֽ⣩��׼������
	cholmod_start (&m_cS) ;  	//��ʼ�� CHOLMOD �ڲ������ݽṹ����Ҫ�Ƿ��乤���ռ䣬������Ĭ�ϲ���
	//m_cS.print_function = NULL;	
	//�ڽ��㷽�̵�ʱ���õ�
	int *Ap  = (int*)malloc((m)*sizeof(int));
	int * Aii = (int*)malloc(m_nS*sizeof(int));//m_nS��S�����з���Ԫ������
	ba_constructAuxCSSGN( Ap, Aii );//����ϡ�����������ṹ����1��1�п�ʼ���洢ÿ�е���ʼ����,�洢����Ԫ�ص�������

	//�� CHOLMOD �m_cholSparseE ��һ�� cholmod_dense �ṹ��
	//typedef struct cholmod_dense
	//{
	//	size_t nrow;  // ����
	//	size_t ncol;  // ����
	//	size_t nzmax; // ������Ԫ������ͨ���� nrow * ncol��
	//	void* x;      // ָ�����ݵ�ָ��
	//} cholmod_dense;

	m_cholSparseE = cholmod_zeros( cnp*m-7, 1, CHOLMOD_REAL, &m_cS);//�̶���1֡����κ͵�2֡��Xc�ᣬʵ��δ֪����������6*m-7
	Ex = (double*)m_cholSparseE->x;
	
	//(m_nS-m_ncams)��ʾS�����зǶԽǿ飨���֮���Լ���������������е�ÿ���ǶԽǿ�ռ��6*6������Ԫ�أ����Գƣ�
	//m_ncams��ʾS�����жԽǿ飨������������Լ����������������ÿ���Խǿ�ռ��21������Ԫ�أ��Գƣ�ֻ�洢�����ǲ��֣�
	nMaxS = (m_nS-m_ncams)*36+m_ncams*21;	//maximum non-zero element in S matrix 

	//����ϡ�����S����4������Ϊtrue��ʾ�����������򣬵�5������Ϊtrue������ܴ洢����6������Ϊ1�������Ĵ洢����Ϊʵ����
	//���ܴ洢��Packed Storage��ָ����ֻ�洢�����з���Ԫ�أ�����������������������ڴ��У���Ԥ������ռ䡣�������Լ����ڴ�ռ�ã���߼���Ч�ʡ�
	m_cholSparseS = cholmod_allocate_sparse(m_ncams*6-7,m_ncams*6-7,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_cholSparseS->p;		//column pointer
	Si = (int*)m_cholSparseS->i;		//row pointer
	
	//Compute initial error and initial reprojection error 
	tStart = clock();
	//dpa[0] = dpa[1] = dpa[2] = dpa[3] = dpa[4] = dpa[5] = 0;//δ�õ�
	//dpa[9+nft] = 0;
	sba_cost(p, hx, m_archor );//hxΪ��ͶӰuv����
	//double sum = 0;
	//for (int i = 0; i < nvis; i += 10000)
	//{
	//	printf("%f %f    %f %f\n", hx[2 * i], hx[2 * i + 1], x[2 * i], x[2 * i + 1]);
	//	sum += (pow(hx[2 * i] - x[2 * i], 2) + pow(hx[2 * i + 1] - x[2 * i + 1], 2));
	//}
	p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
	//printf("%f %f\n", sum, p_eL2);
	initialerror = p_eL2/nvis;
	lasterror = initialerror;
	printf("Initial Error %0.1lf [%0.8lf]\n", p_eL2, initialerror);

	if( m_szReport!= NULL )
		fprintf( fpRe, "Initial Error  %0.8lf\n", initialerror );
	init_p_eL2=p_eL2;	
	
	//Iteration 
	for(itno=0; itno<m_nMaxIter && !stop; ++itno)
	{
		//Setup S matrix, include two step
		memset( U, 0, nuis*Usz*sizeof(double) );//nuis��U����ķ���Ԫ�ظ�����Usz=6*6���洢U����ķ���Ԫ��
		memset( ea, 0,m*easz*sizeof(double) );//m�����������easz=6��ea�洢�����������
		memset( V, 0, n*Vsz*sizeof(double));//n����ά��������Vsz=3*3��
		memset( eb, 0, n*ebsz*sizeof(double));//n����ά��������ebsz=3��eb�洢��ά�������
		memset( W, 0, nvis*Wsz*sizeof(double));	//nvis���������Wsz=6*3

		//Step one; compute W V U directly, don't save each projection image Jacobian 
		t0 = clock();
		if ( m_bRobustKernel)
			sba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			sba_jacobian(p, m_archor, &Uidxij, e, U, ea, V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature);


		t1 = clock();

		ba_inverseVGN( U, V, Uidxij ); //compute inverse matrix of V

		//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
		memset( E, 0, m*easz*sizeof(double));
		memset( S, 0, m_nS*Sblsz*sizeof(double) );
		ba_constructSGN( S, E, U, V, W, ea, eb, Sidxij );//SΪ������࣬EΪ�����Ҳ�
		t2 = clock();
	
		//Solve equation
		ba_constructCSSGN( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init, nft ); //set CSS format using S matrix//break
		for ( ii = 0; ii < 3+nft; ii++ )
			Ex[ii] = E[6+ii];//��2֡����̬��

		for ( ii = 0; ii < 6*m-(10+nft); ii++  )
			Ex[3+nft+ii] = E[10+nft+ii];//������2֡��Xc��ֻ�洢Yc��Zc

		//SxΪ��࣬ExΪ�Ҳ�

		ba_solveCholmodGN( Ap, Aii, init, ordering);//break
		init = true;
		rx = (double*)m_cholSparseR->x;//������������

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
				dpa[ii+10+nft] = rx[3+nft+ii];//�洢������������
		}					
		t3 = clock();
		
		//Solve features
		ba_solveFeatures( W, V, ea, eb, dpa, dpb );//break
		
		//update
		for(i=0, dp_L2=0.0; i<nvars; ++i)//nvars��δ֪��������
		{
			p[i]=p[i] + (tmp=dp[i]);//dpaָ��dp�е���θ�������dpbָ��dp�е���ά���������pָ��δ֪����
			dp_L2+=tmp*tmp;//�洢��������ƽ����
		}

		//dp_L2=sqrt(dp_L2);
		if (dp_L2<1E-8*1E-8)//break
		{	stop = 1;	goto iterstop;	}	//�����ı仯�Ѿ���С						

		ba_updateKR( m_KR, m_KdA, m_KdB, m_KdG, m_K, p );
		t4 = clock();
		
		sba_cost(p, hx, m_archor ); 
		pdp_eL2=nrmL2xmy(e, x, hx, nobs); 				
		error = pdp_eL2/nvis;//break ����mse
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
		{	stop = 2;	goto iterstop;		}//���ϴ�mse��ȣ��Ƿ���������
		
		tPerTime = (t5-t0)*0.001;//0.001�����صļ�ʱֵת������
		tTotalTime += tPerTime;//���ε�������

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
	sba_saveXYZ( m_szCamePose, m_sz3Dpts, p );

	free(S);	free(W);	free(U);	free(V);	
	free(e);	free(eab);	free(E);   	free(dp);	free(hx);	
	free(m_KR); free(m_KdA);free(m_KdB);free(m_KdG);

	printf( "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
		nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
	printf( "SparseBA reasons is listed as following:\n " );
	printf( "reason 1: relative change of state vector is small\n" );
	printf( "reason 2: relative change of projection error is small\n " );
	printf( "reason 3: maximum iteration\n");

	if( m_szReport != NULL )
	{
		fprintf( fpRe, "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
		fprintf( fpRe, "SparseBA reasons is listed as following:\n" );
		fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
		fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
		fprintf( fpRe, "reason 3: maximum iteration\n");
		fclose(fpRe);
	}
	
	return true;
}
void SBA::sba_saveXYZ(const char* camera, const char* sz3Dpt, double* p)
{
	static int i = 0;
	double dx, dy, dz;
	FILE* fp = NULL, * fpc = NULL;

	//save camera poss
	if (camera != NULL)
	{
		fpc = fopen(camera, "w");

		for (i = 0; i < m_ncams; i++)
		{
			fprintf(fpc, "%0.5lf     %0.5lf      %0.5lf     %0.5lf       %0.5lf     %0.5lf\n",
				*(p + i * 6), *(p + i * 6 + 1), *(p + i * 6 + 2), *(p + i * 6 + 3), *(p + i * 6 + 4), *(p + i * 6 + 5));
		}
		fclose(fpc);
	}

	//save features xyz
	if (sz3Dpt != NULL)
	{
		fp = fopen(sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			dx = *(p + m_ncams * 6 + i * 3);
			dy = *(p + m_ncams * 6 + i * 3 + 1);
			dz = *(p + m_ncams * 6 + i * 3 + 2);
			fprintf(fp, "%0.5lf     %0.5lf     %0.5lf\n", dx, dy, dz);
		}
		fclose(fp);
	}

}
void SBA::sba_readProjectionAndInitilizeFeature(FILE *fp,
	double *projs, int ncams, int *archor,char* umask,int* nphoto, int* nfeature )
{
	int n;
	int nframes;
	int ptno = 0;

	int nproj2D = 0;

	int frameno;
	int feastart = 0;

	int nP, nP2;
	
	double* ptr1 = projs;//��¼ÿ��ͶӰ���xy��������ŵģ�2��2��������¼�ʹ洢��

	int i, j; 
	int  sum, cnp = 6;
	
	int *ptr2;

	m_smask = (char*)malloc(m_ncams*m_ncams*sizeof(char));
	memset( m_smask, 0, m_ncams*m_ncams*sizeof(char) );

	while(!feof(fp))//һ��һ�ж�Match��
	{
		n = readNInts( fp, &nframes, 1 );  //��ȡMatch-FeaturePoint�ļ�ÿһ�еĵ�һ�У���ͬ����ĸ���Ϊnframes/3D���2DͶӰ����
		if( n!=1 )
			break;//һ�㶼��1����ʾ�ɹ�
		
		archor[ptno*3] = nframes;//��archor[0]��ʼ����¼3D���2DͶӰ����nframes

		for( i=0, sum = 0; i<nframes; ++i )//һ��ͶӰ��һ��ͶӰ��Ķ�ȡ��������nframes
		{
			n = readNInts( fp, &frameno, 1 ); //�ڶ��Σ���ȡ��һ��ͶӰ������-image index

			nphoto[nproj2D] = frameno;//����һ��ͶӰ���id������nphoto��nproj2D��0��ʼ�ģ�
			nfeature[nproj2D] = ptno;//ptnoҲ�Ǵ�0��ʼ�ģ�Ҳ����nfeature[0]=0?
			nproj2D++;

			if(frameno>=ncams)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "SparseBA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles( fp, ptr1, 2 ); //�ڶ��Σ���ȡ��һ��ͶӰ���������� ptr1[0] ptr1[1]
			ptr1+=2;//һ�μ�¼��������+2
		//	//Ҳ����1��id+2������
			if(n!=3)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "SparseBA:reading image projections wrong!\n");
				return;
			}		
		}
		
		//set masks for U and S matrix              umask�洢U����m_smask�洢S����
		for( i = 0; i < nframes; i++ )
		{
			nP = nphoto[feastart+i];                         //��i���۲����ͼ��

			umask[nP*ncams+nP] = 1;//����ǰê�㣬��ǰê�㣩

			for ( j = i; j < nframes; j++  )
			{
				nP2 = nphoto[feastart+j];//��j���۲����ͼ��                     

				if ( nP == nP2 )                              //
					m_smask[nP*m_ncams+nP2] = 1;//����i���۲����ͼ�ţ���j���۲����ͼ�ţ�
				else if ( nP < nP2 )
					m_smask[nP*m_ncams+nP2] = 1;//����i���۲����ͼ�ţ���j���۲����ͼ�ţ�
				else
				{
					m_smask[nP2*m_ncams + nP] = 1;//����j���۲����ͼ�ţ���i���۲����ͼ�ţ�δִ�й�
				}
					
			}
		}					
		feastart += nframes;
		ptno++;//3D���������ţ�point NO. OR photo NO.
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

void SBA::sba_readAndInitialize(char *camsfname, char *ptsfname, int *ncams,
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
	fpc		=	fopen( camsfname, "r" );//���ⷽλԪ���ļ�
	*ncams	=	findNcameras( fpc );//��ȡ��Ƭ������
	m_ncams =	*ncams;

	fpp		=	fopen( ptsfname, "r" );//��ƥ���ͬ�����ļ�
	readNpointsAndNprojections( fpp, n3Dpts, 3, n2Dprojs, 2 );//��ȡ3D������ĸ�����ͶӰ�����

	*motstruct	=	(double *)malloc( (*ncams*6 + *n3Dpts*3)*sizeof(double) );//���ڴ洢���е��ⷽλԪ�� �� �����������3D����
	if(	*motstruct==NULL )
	{
		fprintf(stderr, "SparseBA error: Memory allocation for 'motstruct' failed \n");
		exit(1);
	}

	*imgpts	=	(double *)malloc(*n2Dprojs*2*sizeof(double));//���ڴ洢����ͶӰ���2D����
	if(	*imgpts==NULL )
	{
		fprintf(stderr, "SparseBA error: Memory allocation for 'imgpts' failed\n");
		exit(1);
	}
	//���Ҫ���ļ���ͷ���������ָ���ƶ����ļ�ͷ�����øú���
	rewind(fpc);//ÿ��ȡһ���ַ����ļ��ڲ�λ��ָ������ƶ�һ���ֽڣ���ȡ��ϣ���ָ����ָ���ļ�ĩβ��
	rewind(fpp);

	//allocate indicator of U
	*umask = (char*)malloc(*ncams * *ncams );
	memset(*umask, 0, *ncams * *ncams * sizeof(char));//��1������һ��Ҫ��һ����֪�ģ��Ѿ��������ڴ�ĵ�ַ����3������һ��Ҫʹ��sizeof��������
	//��һ���Ѿ������ַ���ڴ���г�ʼ��������ͨ����ʼ��0�����ַ�'\0'

	//allocate main and associate anchors
	*archor = (int*)malloc(*n3Dpts*3*sizeof(int));//
	memset( *archor, -1, *n3Dpts*3*sizeof(int) ); //��������ê���λ�ó�ʼ��Ϊ-1

	*nphoto		= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*nfeature	= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*archorSort = (int*)malloc(*n3Dpts*3*sizeof(int));

	ba_readCameraPose(fpc, *motstruct);

	fclose(fpc);
	
	//Update KR
	m_KR  = (double*)malloc(m_ncams*9*sizeof(double));//����ڲξ��� �� ��ת����ĳ˻�
	m_KdA = (double*)malloc(m_ncams*9*sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	m_KdB = (double*)malloc(m_ncams*9*sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	m_KdG = (double*)malloc(m_ncams*9*sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	
	if (savePara == 1)
	{
		//����P����
		m_P = (double*)malloc(m_ncams * 12 * sizeof(double));
		ba_constructP(m_P, m_K, *motstruct);
	}
	else
		ba_updateKR(m_KR, m_KdA, m_KdB, m_KdG, m_K, *motstruct);//��5������Ϊ����ڲΣ���6������Ϊ������


	//if XYZ are provided, we can use them as feature initialization.
	if (m_bProvideXYZ)
	{
		fpXYZ = fopen( m_szXYZ, "r");
		m_XYZ = (double*)malloc(m_n3Dpts*3*sizeof(double));
		for (i = 0; i < m_n3Dpts; i++)
		{
			fscanf(fpXYZ, "%lf  %lf  %lf", m_XYZ + i * 3, m_XYZ + i * 3 + 1, m_XYZ + i * 3 + 2);
			(*motstruct)[m_ncams * 6 + i * 3 + 0] = *(m_XYZ + i * 3 + 0);
			(*motstruct)[m_ncams * 6 + i * 3 + 1] = *(m_XYZ + i * 3 + 1);
			(*motstruct)[m_ncams * 6 + i * 3 + 2] = *(m_XYZ + i * 3 + 2);
		}
		fclose(fpXYZ);
	}

	if(savePara ==1)
		ba_readProjectionAndTriangulateFeature(fpp, *imgpts, *ncams);
	else
	{
		//ֻ��������umask��smask
		sba_readProjectionAndInitilizeFeature(fpp,
			*imgpts,
			*ncams,
			*archor,
			*umask,
			*nphoto,
			*nfeature);
	}

	fclose(fpp);//���������ر�

	if(m_bProvideXYZ)
		free(m_XYZ);

	fpXYZ = NULL;
}


void SBA::sba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
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
		numfea = m_archor[i*3];//��ά��Ĺ���ͼ����
		cur = 0;
		nF = i;
		//printf("%d %d\n", m_archorSort[i * 2], m_archorSort[i * 2 + 1]);
		
		for ( j = 0; j < numfea; j++ )
		{
			
			nP = nphoto[pos2];//����ͼ���
			ppt=pb + nF*pnp;//��ά������  
			pe= e + pos2*2;		//error

			//nM = archor[nF*3+1];//��ê��
			//nN = archor[nF*3+2];//��ê��	

			//paָ���������
			//pptָ����ά�����
			//nP��ǰ��ͼ���
			//nM��ê����
			//nA��ê����
			//pPAָ��uv�Ե�ǰ��ͼ������һ�׵�
			//pPBָ��uv����ά���һ�׵�
			//pAMָ��uv������ê����ͼ��uv������ê����ͼλ�Ʋ�����һ�׵�
			//pAAָ��uv��������ê����ͼ��uv���Ը�ê����ͼλ�Ʋ�����һ�׵�

			sba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nP, pPA, pPB );
		
			//����
			//double matXj[3], matxyt[3], n1[2], n2[2];
			//double *ppa = pa + nP * 6 + 3;//��ǰ��ͼλ��ָ��
			//matXj[0] = ppt[0] - ppa[0];
			//matXj[1] = ppt[1] - ppa[1];
			//matXj[2] = ppt[2] - ppa[2];
			//double* pKR = m_KR + nP * 9;
			//matxyt[0] = pKR[0] * matXj[0] + pKR[1] * matXj[1] + pKR[2] * matXj[2];
			//matxyt[1] = pKR[3] * matXj[0] + pKR[4] * matXj[1] + pKR[5] * matXj[2];
			//matxyt[2] = pKR[6] * matXj[0] + pKR[7] * matXj[1] + pKR[8] * matXj[2];
			//n1[0] = matxyt[0] / matxyt[2];
			//n1[1] = matxyt[1] / matxyt[2];
			//double sigma = 1e-6;
			////double* ppa_new = pa + nP * 6;
			//double* ppa_new = ppa;
			//ppa_new[1] += sigma;
			////for (int hh = 0; hh < 90; hh++)
			////{
			////	printf("%f %f %f %f %f %f\n", pa[hh * 6], pa[hh * 6 + 1], pa[hh * 6 + 2], pa[hh * 6 + 3], pa[hh * 6 + 4], pa[hh * 6 + 5]);
			////}
			//sba_updateKR(m_KR, m_KdA, m_KdB, m_KdG, m_K, pa);
			////printf("%f %f %f\n", ppt[0], ppt[1], ppt[2]);
			//matXj[0] = ppt[0] - ppa[0];
			//matXj[1] = ppt[1] - ppa[1];
			//matXj[2] = ppt[2] - ppa[2];
			//matxyt[0] = pKR[0] * matXj[0] + pKR[1] * matXj[1] + pKR[2] * matXj[2];
			//matxyt[1] = pKR[3] * matXj[0] + pKR[4] * matXj[1] + pKR[5] * matXj[2];
			//matxyt[2] = pKR[6] * matXj[0] + pKR[7] * matXj[1] + pKR[8] * matXj[2];
			//n2[0] = matxyt[0] / matxyt[2];
			//n2[1] = matxyt[1] / matxyt[2];
			//double t1 = (n2[0] - n1[0]) / sigma;
			//double t2 = (n2[1] - n1[1]) / sigma;
			//printf("%f %f\n", pPA[0 * 6 + 4], t1);
			//printf("%f %f\n", pPA[1 * 6 + 4], t2);

			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);//����U������nP��nP�е�Ԫ����ϡ��ṹ�е�������
			//printf("%d  ", pos);
			ppUpa = U + pos*(6*6);//������䣬�᲻�����ͬһ�鱻�ظ��������⣿ͬһ���ǲ�ͬ�۲�ƫ������ӹ��̣����ᱻ�ظ����

			//printf("\n\n\n\n");
			for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//�����ǶԳƾ���ֻ��������
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
			pos2++;
			cur++;
		}
	}
}

void SBA::sba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
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

			sba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nP, pPA, pPB );

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
/*void SBA::sba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
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

void SBA::sba_cost(double* p, double* hx, int* archor)
{
	register int i;
	int cnp, pnp, mnp;
	double* pa, * pb, * ppt, * pmeas, * ppa;
	int nF, nP;
	int m=m_ncams;

	cnp = 6, pnp = 3, mnp = 2;
	pa = p; pb = p + m * cnp;

	for (i = 0; i < m_n2Dprojs; i++)
	{
		nF = m_feature[i];
		nP = m_photo[i];

		ppt = pb + nF * pnp;//��ά������
		ppa = pa + nP * cnp + 3;
		pmeas = hx + i * mnp; // set pmeas to point to hx_ij �洢��ͶӰ���ص�
		
		//KR�������ָ�룬��ά��ָ�룬��ǰê�㣬��ͶӰuv
		sba_reprojectEachPts(m_KR, ppa, ppt, nP, pmeas);
		//printf("%f %f\n", pmeas[0], pmeas[1]);
	}
}
