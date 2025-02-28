#include "stdafx.h"
#include "IBA.h"

//cite from SBA code in order to search sub matrix of S according to (i,j)
void IBA::sba_crsm_alloc(struct sba_crsm* sm, int nr, int nc, int nnz)
{
	int msz;
	sm->nr = nr;
	sm->nc = nc;
	sm->nnz = nnz;
	msz = 2 * nnz + nr + 1;
	sm->val = (int*)malloc(msz * sizeof(int));  /* required memory is allocated in a single step */
	if (!sm->val) {
		fprintf(stderr, "memory allocation request failed in sba_crsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
		exit(1);
	}
	sm->colidx = sm->val + nnz;
	sm->rowptr = sm->colidx + nnz;
}

void IBA::sba_crsm_free(struct sba_crsm* sm)
{
	sm->nr = sm->nc = sm->nnz = -1;
	free(sm->val);
	sm->val = sm->colidx = sm->rowptr = NULL;
}

/* returns the index of the (i, j) element. No bounds checking! */
int IBA::sba_crsm_elmidx(struct sba_crsm* sm, int i, int j)//返回i行j列元素的列索引
{
	register int low, high, mid, diff;

	low = sm->rowptr[i];
	high = sm->rowptr[i + 1] - 1;

	/* binary search for finding the element at column j */
	while (low <= high)
	{
		mid = (low + high) >> 1; //(low+high)/2;
		diff = j - sm->colidx[mid];
		if (diff < 0)//在j右侧
			high = mid - 1;
		else if (diff > 0)//在j左侧
			low = mid + 1;
		else
			return mid;
	}

	return -1; /* not found */
}
//计算文件中非注释行的数量
int IBA::findNcameras(FILE* fp)
{
	int lineno, ncams, ch;

	lineno = ncams = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#') { /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);

		SKIP_LINE(fp);
		++lineno;
		if (ferror(fp))
		{
			fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		++ncams;
	}
	return ncams;
}

IBA::IBA(void)
{
}

IBA::~IBA(void)
{
}

int IBA::countNDoubles(FILE* fp)
{
	int lineno, ch, np, i;
	char buf[MAXSTRLEN], * s;
	double dummy;

	lineno = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) return 0;

		ungetc(ch, fp);
		++lineno;
		if (!fgets(buf, MAXSTRLEN - 1, fp)) { /* read the line found... */
			fprintf(stderr, "countNDoubles(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		/* ...and count the number of doubles it has */
		for (np = i = 0, s = buf; 1; ++np, s += i) {
			ch = sscanf(s, "%lf%n", &dummy, &i);
			if (ch == 0 || ch == EOF) break;
		}

		rewind(fp);
		return np;
	}
	return 0; // should not reach this point
}

int IBA::skipNDoubles(FILE* fp, int nvals)
{
	register int i;
	int j;

	for (i = 0; i < nvals; ++i)
	{
		j = fscanf(fp, "%*f");
		if (j == EOF) return EOF;

		if (ferror(fp)) return EOF - 1;
	}

	return nvals;
}

void IBA::readNpointsAndNprojections(FILE* fp, int* n3Dpts, int pnp, int* nprojs, int mnp)
{
	int nfirst, lineno, npts, nframes, ch, n;

	/* #parameters for the first line */
	nfirst = countNDoubles(fp);

	*n3Dpts = *nprojs = lineno = npts = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n = readNInts(fp, &nframes, 1);
		if (n != 1)
			exit(1);

		SKIP_LINE(fp);
		*nprojs += nframes;
		++npts;
	}

	*n3Dpts = npts;
}

void IBA::ba_readCablibration(FILE* fp, double* K)
{
	int n = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &K[0], &K[3], &K[6], &K[1], &K[4], &K[7], &K[2], &K[5], &K[8]);

	if (n != 9)
	{
		fprintf(stderr, "BA error: Format of Calibaration is wrong\n");
		exit(1);
	}
}

void IBA::ba_readCameraPose(FILE* fp, double* params)
{
	int n, num, lineno = 0;
	double* tofilter;
	double* pPrams = params;

	//the number of element per line is 8, it represents that focal length vary, or it is constant
	num = countNDoubles(fp);
	if (num == 8)
	{
		m_bFocal = true;
		m_K = (double*)malloc(m_ncams * 2 * sizeof(double));
		tofilter = (double*)malloc(8 * sizeof(double));
	}
	else
		tofilter = (double*)malloc(6 * sizeof(double));

	while (!feof(fp))
	{
		if (num == 6)
		{
			n = readNDoubles(fp, tofilter, 6);
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];
			/*------------------------------------调试用-------------------------------------------*/
			//for (int i = 0;i < 6;i++)
			//	printf("%f ", tofilter[i]);
			//printf("\n");
			/*-------------------------------------调试用------------------------------------------*/
		}
		if (num == 8)
		{
			n = readNDoubles(fp, tofilter, 8);
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];

			m_K[lineno * 2] = tofilter[6];
			m_K[lineno * 2 + 1] = tofilter[7];
		}

		if (n == -1)
			break;
		pPrams += 6;
		++lineno;
	}
}

int IBA::readNInts(FILE* fp, int* vals, int nvals)
{
	register int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i) {//第一次：nvals就是1
		j = fscanf(fp, "%d", vals + i);//第一次：读取第一个，也就是1割三维点的所有2D投影总数
		if (j == EOF) return EOF;//第一次：一般不会是，会是j=1，表示成功

		if (j != 1 || ferror(fp)) return EOF - 1;//第一次：一般都是，j=1，表示成功

		n += j;//第一次因为n=0，所以就是j，也就是成功标志j=1
	}//之后的话，fscanf每读一次那么就会移动指针，然后%d会会跳过空格、换行符和其他空白字符，直到它遇到%d的输入数据
	//或者遇到了EOF或者j!=0了
	return n;
}

int IBA::readNDoubles(FILE* fp, double* vals, int nvals)
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
void IBA::ba_readCameraPoseration(char* fname, double* ical)
{
	FILE* fp;
	int  ch = EOF;

	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "BA: Cannot open calbration file %s, exiting\n", fname);
		return;
	}
	//按列读取1-2-3 || 4-5-6 || 7-8-9
	int num = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &ical[0], &ical[3], &ical[6], &ical[1], &ical[4], &ical[7], &ical[2], &ical[5], &ical[8]);
	if (num != 9)
	{
		fprintf(stderr, "BA error: Format of Calibration file is wrong");
		return;
	}

	fclose(fp);
}

void IBA::ba_updateKR(double* KR, double* KdA, double* KdB, double* KdG, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pKR, * pKdA, * pKdB, * pKdG;
		double matR[9];//相对旋转矩阵
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];

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
			pKR = KR + i * 9;

			//for (int k1 = 0; k1 < 3; k1++)
			//{
			//	for (int k2 = 0; k2 < 3; k2++)
			//	{
			//		printf("%f ", K[k1 * 3 + k2]);
			//	}
			//	printf("\n");
			//}
			pKR[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pKR[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pKR[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pKR[3] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pKR[4] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pKR[5] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pKR[6] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pKR[7] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pKR[8] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];

			//pKR矩阵关于omega的一阶导
			pKdG = KdG + i * 9;
			tmp1[0] = K[0] * matDRG[0] + K[3] * matDRG[3] + K[6] * matDRG[6];
			tmp1[1] = K[1] * matDRG[0] + K[4] * matDRG[3] + K[7] * matDRG[6];
			tmp1[2] = K[2] * matDRG[0] + K[5] * matDRG[3] + K[8] * matDRG[6];
			tmp1[3] = K[0] * matDRG[1] + K[3] * matDRG[4] + K[6] * matDRG[7];
			tmp1[4] = K[1] * matDRG[1] + K[4] * matDRG[4] + K[7] * matDRG[7];
			tmp1[5] = K[2] * matDRG[1] + K[5] * matDRG[4] + K[8] * matDRG[7];
			tmp1[6] = K[0] * matDRG[2] + K[3] * matDRG[5] + K[6] * matDRG[8];
			tmp1[7] = K[1] * matDRG[2] + K[4] * matDRG[5] + K[7] * matDRG[8];
			tmp1[8] = K[2] * matDRG[2] + K[5] * matDRG[5] + K[8] * matDRG[8];

			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdG[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdG[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdG[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdG[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdG[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdG[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdG[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdG[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdG[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//pKR矩阵关于phi的一阶导
			pKdB = KdB + i * 9;
			tmp1[0] = K[0] * matRG[0] + K[3] * matRG[3] + K[6] * matRG[6];
			tmp1[1] = K[1] * matRG[0] + K[4] * matRG[3] + K[7] * matRG[6];
			tmp1[2] = K[2] * matRG[0] + K[5] * matRG[3] + K[8] * matRG[6];
			tmp1[3] = K[0] * matRG[1] + K[3] * matRG[4] + K[6] * matRG[7];
			tmp1[4] = K[1] * matRG[1] + K[4] * matRG[4] + K[7] * matRG[7];
			tmp1[5] = K[2] * matRG[1] + K[5] * matRG[4] + K[8] * matRG[7];
			tmp1[6] = K[0] * matRG[2] + K[3] * matRG[5] + K[6] * matRG[8];
			tmp1[7] = K[1] * matRG[2] + K[4] * matRG[5] + K[7] * matRG[8];
			tmp1[8] = K[2] * matRG[2] + K[5] * matRG[5] + K[8] * matRG[8];

			tmp2[0] = tmp1[0] * matDRB[0] + tmp1[3] * matDRB[3] + tmp1[6] * matDRB[6];
			tmp2[1] = tmp1[1] * matDRB[0] + tmp1[4] * matDRB[3] + tmp1[7] * matDRB[6];
			tmp2[2] = tmp1[2] * matDRB[0] + tmp1[5] * matDRB[3] + tmp1[8] * matDRB[6];
			tmp2[3] = tmp1[0] * matDRB[1] + tmp1[3] * matDRB[4] + tmp1[6] * matDRB[7];
			tmp2[4] = tmp1[1] * matDRB[1] + tmp1[4] * matDRB[4] + tmp1[7] * matDRB[7];
			tmp2[5] = tmp1[2] * matDRB[1] + tmp1[5] * matDRB[4] + tmp1[8] * matDRB[7];
			tmp2[6] = tmp1[0] * matDRB[2] + tmp1[3] * matDRB[5] + tmp1[6] * matDRB[8];
			tmp2[7] = tmp1[1] * matDRB[2] + tmp1[4] * matDRB[5] + tmp1[7] * matDRB[8];
			tmp2[8] = tmp1[2] * matDRB[2] + tmp1[5] * matDRB[5] + tmp1[8] * matDRB[8];

			pKdB[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdB[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdB[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdB[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdB[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdB[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdB[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdB[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdB[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//pKR矩阵关于kappa的一阶导
			pKdA = KdA + i * 9;
			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdA[0] = tmp2[0] * matDRA[0] + tmp2[3] * matDRA[3] + tmp2[6] * matDRA[6];
			pKdA[3] = tmp2[1] * matDRA[0] + tmp2[4] * matDRA[3] + tmp2[7] * matDRA[6];
			pKdA[6] = tmp2[2] * matDRA[0] + tmp2[5] * matDRA[3] + tmp2[8] * matDRA[6];
			pKdA[1] = tmp2[0] * matDRA[1] + tmp2[3] * matDRA[4] + tmp2[6] * matDRA[7];
			pKdA[4] = tmp2[1] * matDRA[1] + tmp2[4] * matDRA[4] + tmp2[7] * matDRA[7];
			pKdA[7] = tmp2[2] * matDRA[1] + tmp2[5] * matDRA[4] + tmp2[8] * matDRA[7];
			pKdA[2] = tmp2[0] * matDRA[2] + tmp2[3] * matDRA[5] + tmp2[6] * matDRA[8];
			pKdA[5] = tmp2[1] * matDRA[2] + tmp2[4] * matDRA[5] + tmp2[7] * matDRA[8];
			pKdA[8] = tmp2[2] * matDRA[2] + tmp2[5] * matDRA[5] + tmp2[8] * matDRA[8];

		}
	}
	else
	{
		int i = 0;
		double* ptAngle;
		double* pKR, * pKdA, * pKdB, * pKdG;
		double matR[9];
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];
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

			K[0] = m_K[i * 2];
			K[4] = m_K[i * 2 + 1];

			pKR = KR + i * 9;
			pKR[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pKR[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pKR[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pKR[3] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pKR[4] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pKR[5] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pKR[6] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pKR[7] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pKR[8] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];

			//KdG
			pKdG = KdG + i * 9;
			tmp1[0] = K[0] * matDRG[0] + K[3] * matDRG[3] + K[6] * matDRG[6];
			tmp1[1] = K[1] * matDRG[0] + K[4] * matDRG[3] + K[7] * matDRG[6];
			tmp1[2] = K[2] * matDRG[0] + K[5] * matDRG[3] + K[8] * matDRG[6];
			tmp1[3] = K[0] * matDRG[1] + K[3] * matDRG[4] + K[6] * matDRG[7];
			tmp1[4] = K[1] * matDRG[1] + K[4] * matDRG[4] + K[7] * matDRG[7];
			tmp1[5] = K[2] * matDRG[1] + K[5] * matDRG[4] + K[8] * matDRG[7];
			tmp1[6] = K[0] * matDRG[2] + K[3] * matDRG[5] + K[6] * matDRG[8];
			tmp1[7] = K[1] * matDRG[2] + K[4] * matDRG[5] + K[7] * matDRG[8];
			tmp1[8] = K[2] * matDRG[2] + K[5] * matDRG[5] + K[8] * matDRG[8];

			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdG[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdG[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdG[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdG[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdG[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdG[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdG[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdG[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdG[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//KdB
			pKdB = KdB + i * 9;
			tmp1[0] = K[0] * matRG[0] + K[3] * matRG[3] + K[6] * matRG[6];
			tmp1[1] = K[1] * matRG[0] + K[4] * matRG[3] + K[7] * matRG[6];
			tmp1[2] = K[2] * matRG[0] + K[5] * matRG[3] + K[8] * matRG[6];
			tmp1[3] = K[0] * matRG[1] + K[3] * matRG[4] + K[6] * matRG[7];
			tmp1[4] = K[1] * matRG[1] + K[4] * matRG[4] + K[7] * matRG[7];
			tmp1[5] = K[2] * matRG[1] + K[5] * matRG[4] + K[8] * matRG[7];
			tmp1[6] = K[0] * matRG[2] + K[3] * matRG[5] + K[6] * matRG[8];
			tmp1[7] = K[1] * matRG[2] + K[4] * matRG[5] + K[7] * matRG[8];
			tmp1[8] = K[2] * matRG[2] + K[5] * matRG[5] + K[8] * matRG[8];

			tmp2[0] = tmp1[0] * matDRB[0] + tmp1[3] * matDRB[3] + tmp1[6] * matDRB[6];
			tmp2[1] = tmp1[1] * matDRB[0] + tmp1[4] * matDRB[3] + tmp1[7] * matDRB[6];
			tmp2[2] = tmp1[2] * matDRB[0] + tmp1[5] * matDRB[3] + tmp1[8] * matDRB[6];
			tmp2[3] = tmp1[0] * matDRB[1] + tmp1[3] * matDRB[4] + tmp1[6] * matDRB[7];
			tmp2[4] = tmp1[1] * matDRB[1] + tmp1[4] * matDRB[4] + tmp1[7] * matDRB[7];
			tmp2[5] = tmp1[2] * matDRB[1] + tmp1[5] * matDRB[4] + tmp1[8] * matDRB[7];
			tmp2[6] = tmp1[0] * matDRB[2] + tmp1[3] * matDRB[5] + tmp1[6] * matDRB[8];
			tmp2[7] = tmp1[1] * matDRB[2] + tmp1[4] * matDRB[5] + tmp1[7] * matDRB[8];
			tmp2[8] = tmp1[2] * matDRB[2] + tmp1[5] * matDRB[5] + tmp1[8] * matDRB[8];

			pKdB[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdB[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdB[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdB[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdB[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdB[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdB[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdB[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdB[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//KdA
			pKdA = KdA + i * 9;
			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdA[0] = tmp2[0] * matDRA[0] + tmp2[3] * matDRA[3] + tmp2[6] * matDRA[6];
			pKdA[3] = tmp2[1] * matDRA[0] + tmp2[4] * matDRA[3] + tmp2[7] * matDRA[6];
			pKdA[6] = tmp2[2] * matDRA[0] + tmp2[5] * matDRA[3] + tmp2[8] * matDRA[6];
			pKdA[1] = tmp2[0] * matDRA[1] + tmp2[3] * matDRA[4] + tmp2[6] * matDRA[7];
			pKdA[4] = tmp2[1] * matDRA[1] + tmp2[4] * matDRA[4] + tmp2[7] * matDRA[7];
			pKdA[7] = tmp2[2] * matDRA[1] + tmp2[5] * matDRA[4] + tmp2[8] * matDRA[7];
			pKdA[2] = tmp2[0] * matDRA[2] + tmp2[3] * matDRA[5] + tmp2[6] * matDRA[8];
			pKdA[5] = tmp2[1] * matDRA[2] + tmp2[4] * matDRA[5] + tmp2[7] * matDRA[8];
			pKdA[8] = tmp2[2] * matDRA[2] + tmp2[5] * matDRA[5] + tmp2[8] * matDRA[8];

		}
	}
}
void IBA::ba_constructP(double* P, double* K, double* p)
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
		double matR[9], matT[3];
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

double IBA::nrmL2xmy(double* const e, const double* const x, const double* const y, const int n)
{
	const int blocksize = 8, bpwr = 3; /* 8=2^3 */
	register int i;
	int j1, j2, j3, j4, j5, j6, j7;
	int blockn;
	register double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

	/* n may not be divisible by blocksize,
	* go as near as we can first, then tidy up.
	*/
	//printf("nobs%d\n", n);//Add by Lu
	blockn = (n >> bpwr) << bpwr; /* (n / blocksize) * blocksize; */
	//printf("cycle_num%d, blocksize%d", blockn - 1, blocksize);//Add by Lu

	/* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
	for (i = blockn - 1; i > 0; i -= blocksize) {
		e[i] = x[i] - y[i]; sum0 += e[i] * e[i];
		j1 = i - 1; e[j1] = x[j1] - y[j1]; sum1 += e[j1] * e[j1];
		j2 = i - 2; e[j2] = x[j2] - y[j2]; sum2 += e[j2] * e[j2];
		j3 = i - 3; e[j3] = x[j3] - y[j3]; sum3 += e[j3] * e[j3];
		j4 = i - 4; e[j4] = x[j4] - y[j4]; sum0 += e[j4] * e[j4];
		j5 = i - 5; e[j5] = x[j5] - y[j5]; sum1 += e[j5] * e[j5];
		j6 = i - 6; e[j6] = x[j6] - y[j6]; sum2 += e[j6] * e[j6];
		j7 = i - 7; e[j7] = x[j7] - y[j7]; sum3 += e[j7] * e[j7];
	}

	i = blockn;
	if (i < n) {
		switch (n - i) {
		case 7: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
		case 6: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
		case 5: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
		case 4: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
		case 3: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
		case 2: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
		case 1: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
		}
	}

	return sum0 + sum1 + sum2 + sum3;
}

void IBA::readNpointsAndNprojectionsFromProj(FILE* fp, int& n3Dpts, int& nprojs)
{
	int nfirst, lineno, npts, nframes, ch, n;
	nprojs = 0;
	n3Dpts = 0;
	npts = 0;

	/* #parameters for the first line */
	nfirst = countNDoubles(fp);

	//*n3Dpts=*nprojs=lineno=npts=0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n = readNInts(fp, &nframes, 1);
		if (n != 1)
		{
			fprintf(stderr, "readNpointsAndNprojections(): error reading input file, line %d: "
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		SKIP_LINE(fp);
		nprojs += nframes;
		++npts;
	}

	n3Dpts = npts;
}

void IBA::readPointProjections(FILE* fp, double* imgpts, int* photo, int* imgptsSum, int n3Dpts, int n2Dprojs)
{
	int nframes, ch, lineno, ptno, frameno, n;
	int i;
	int nproj2D = 0;

	lineno = ptno = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			lineno++;

			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);

		n = readNInts(fp, &nframes, 1);  /* read in number of image projections */
		if (n != 1)
		{
			fprintf(stderr, "sba_readProjectionAndInitilizeFeature(): error reading input file, line %d:\n"
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		imgptsSum[ptno] = nframes;

		for (i = 0; i < nframes; ++i)
		{
			n = readNInts(fp, &frameno, 1); /* read in frame number... */

			photo[nproj2D] = frameno;

			n += readNDoubles(fp, imgpts + nproj2D * 2, 2); /* ...and image projection */

			nproj2D++;
		}
		fscanf(fp, "\n"); // consume trailing newline

		lineno++;
		ptno++;
	}
}
void IBA::readImagePts(const char* szProj, double** imgpts, int** photo, int** imgptsSum, int& n3Dpts, int& n2Dprojs)
{
	FILE* fpp;
	if ((fpp = fopen(szProj, "r")) == NULL) {
		fprintf(stderr, "cannot open file %s, exiting\n", szProj);
		exit(1);
	}
	readNpointsAndNprojectionsFromProj(fpp, n3Dpts, n2Dprojs);

	*imgpts = (double*)malloc(n2Dprojs * 2 * sizeof(double));
	if (*imgpts == NULL) {
		fprintf(stderr, "memory allocation for 'imgpts' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*photo = (int*)malloc(n2Dprojs * sizeof(int));
	if (*photo == NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*imgptsSum = (int*)malloc(n3Dpts * sizeof(int));
	if (*imgptsSum == NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	rewind(fpp);
	readPointProjections(fpp, *imgpts, *photo, *imgptsSum, n3Dpts, n2Dprojs);

	fclose(fpp);
}
bool IBA::ba_parseArgs(int argc, char* argv[])
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

		if (strcmp(name.c_str(), "help") == 0)
		{
			ba_printHelp(pba);
			return false;
		}

		if (strcmp(name.c_str(), "cam") == 0)
		{
			i++;
			param = argv[i];
			m_szCameraInit = (char*)malloc(param.length());
			strcpy(m_szCameraInit, param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "fea") == 0)
		{
			i++;
			param = argv[i];
			m_szFeatures = (char*)malloc(param.length());
			strcpy(m_szFeatures, param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "calib") == 0)
		{
			i++;
			param = argv[i];
			m_szCalibration = (char*)malloc(param.length());
			strcpy(m_szCalibration, param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "pose") == 0)
		{
			i++;
			param = argv[i];
			m_szCamePose = (char*)malloc(param.length());
			strcpy(m_szCamePose, param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "3D") == 0)
		{
			i++;
			param = argv[i];
			m_sz3Dpts = (char*)malloc(param.length());
			strcpy(m_sz3Dpts, param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "report") == 0)
		{
			i++;
			param = argv[i];
			m_szReport = (char*)malloc(param.length());
			strcpy(m_szReport, param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "xyz") == 0)
		{
			i++;
			param = argv[i];
			m_szXYZ = (char*)malloc(param.length());
			strcpy(m_szXYZ, param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "i") == 0)
		{
			i++;
			param = argv[i];
			m_nMaxIter = atoi(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "robustKernel") == 0)
		{
			m_bRobustKernel = true;
			bRKF = false;
			i++;
			param = argv[i];
			if (strcmp(param.c_str(), "Huber") == 0)
			{
				m_nRobustType = 2;
				bRKF = true;
			}

			if (strcmp(param.c_str(), "Cauchy") == 0)
			{
				m_nRobustType = 1;
				bRKF = true;
			}

			bSuccess = true;

			if (!bRKF)
			{
				printf("BA: Must input right robust kernel function!\n");
				return false;
			}
		}

		if (strcmp(name.c_str(), "solve") == 0)
		{
			i++;
			param = argv[i];
			if (strcmp(param.c_str(), "LM") == 0)
				m_bsolverLM = true;

			if (strcmp(param.c_str(), "GN") == 0)
				m_bsolverGN = true;

			bSuccess = true;
		}

		if (strcmp(name.c_str(), "t") == 0)
		{
			i++;
			param = argv[i];

			m_Tau = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e1") == 0)
		{
			i++;
			param = argv[i];

			m_e1 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e2") == 0)
		{
			i++;
			param = argv[i];

			m_e2 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e3") == 0)
		{
			i++;
			param = argv[i];

			m_e3 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e4") == 0)
		{
			i++;
			param = argv[i];

			m_e4 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "robustKernelWidth") == 0)
		{
			i++;
			param = argv[i];

			m_delt = atof(param.c_str());
			bSuccess = true;
		}

		if (!bSuccess)
		{
			printf("BA error: %s command is wrong!\n", name.c_str());
			return false;
		}

	}

	return true;

}

void IBA::ba_printHelp(BAType ba)
{
	string name;
	if (ba == pba)
		name = "Parallax";
	else if (ba == sba)
		name = "Sparse";
	printf(name.append("Sparse Bundle Adjustment General Options\n").c_str());
	printf("\n");

	printf("-cam			Provide initial camera pose.\n");
	printf("-fea			Provide features.\n");
	printf("-calib			Provide calibration.\n");
	printf("-xyz			Provide initial XYZ.\n");
	printf("-pose			Output optimal camera pose.\n");
	printf("-3D			Output optimal 3D point cloud.\n");
	printf("-report			Output report.\n");
	printf("-solve			Solve method including LevenbergMarquart(LM) and Gauss-Newton(GN).\n");
	printf("-i			Max Iteration.\n");
	printf("-t			LevenbergMarquart parameters.\n");
	printf("-robustKernel		use Cauchy Robust Kernel Function.\n");
	printf("-robustKernelWidth		width for the robust Kernel.\n");
}

void IBA::ba_saveTriangulatedxyz(const char* sz3Dpt, double* p)
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

int IBA::ba_ConstructSmask(sba_crsm& Sidxij, sba_crsm& Uidxij)//分别是S矩阵和U矩阵的稀疏索引
{
	int i, j, k, ii, jj;
	int nuis, m = m_ncams;//nuis是U矩阵的非零元素个数
	//compute total smask
	for (i = 0; i < m; i++) for (j = 0; j < m; j++)
	{
		//这段代码未执行
		if (m_umask[i * m + j] == 1 && m_smask[i * m + j] == 0)//m_umask[i*m+j] == 1 表示 U 矩阵在该位置是非零的，如果 S 还没标记，就把它标记上。
		{
			m_smask[i * m + j] = 1;//S 矩阵的稀疏掩码 (Sidxij 的结构),m_smask[i*m+j] 是一个 m×m 维的二维数组（用 1D 存储），表示 S 矩阵在 (i,j) 位置是否有非零元素。
			m_nS += 1;
		}
	}
	//分配 S 矩阵的存储空间,Sidxij是S 矩阵的结构体,m, m为矩阵大小
	sba_crsm_alloc(&Sidxij, m, m, m_nS);//m_nS是S矩阵的非零元素个数
	for (i = k = 0; i < m; ++i)
	{
		Sidxij.rowptr[i] = k;// 记录 S 矩阵第 i 行的起始索引
		ii = i * m;
		for (j = 0; j < m; ++j)
			if (m_smask[ii + j])// 如果该位置是非零元素
			{
				Sidxij.val[k] = k;// 这里存的值是索引 k
				Sidxij.colidx[k++] = j; // 记录该非零元素的列索引
			}
	}
	Sidxij.rowptr[m] = m_nS;// 最后一行指向末尾

	for (i = nuis = 0, jj = m * m; i < jj; ++i)
		nuis += (m_umask[i] != 0);//U 矩阵的稀疏掩码 (Uidxij 的结构)

	sba_crsm_alloc(&Uidxij, m, m, nuis);
	for (i = k = 0; i < m; ++i)
	{
		Uidxij.rowptr[i] = k;
		ii = i * m;
		for (j = 0; j < m; ++j)
			if (m_umask[ii + j])
			{
				Uidxij.val[k] = k;
				Uidxij.colidx[k++] = j;
			}
	}
	Uidxij.rowptr[m] = nuis;

	//测试
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

void IBA::ba_inverseVLM(double* V, double* IV, sba_crsm& Uidxij, double mu)
{
	int i, j;
	int m = m_ncams, n = m_n3Dpts;
	int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
	double* ptr1, * ptr2;
	Matrix3d MatInv;

	//IV save inverse V matrix, V must unchange for the next step
	memcpy(IV, V, n * Vsz * sizeof(double));
	for (i = 0; i < n; ++i)
	{
		ptr1 = V + i * Vsz;
		ptr2 = IV + i * Vsz;

		for (j = 0; j < pnp; ++j)
			ptr2[j * pnp + j] += mu;

		Eigen::Matrix3d matV(ptr2);
		MatInv = matV.inverse();
		ptr2[0] = MatInv(0, 0);
		ptr2[4] = MatInv(1, 1);
		ptr2[8] = MatInv(2, 2);
		ptr2[1] = ptr2[3] = MatInv(0, 1);
		ptr2[2] = ptr2[6] = MatInv(0, 2);
		ptr2[5] = ptr2[7] = MatInv(1, 2);
	}
}
void IBA::ba_inverseVGN(double* U, double* V, sba_crsm& Uidxij)
{
	int i;
	int m = m_ncams, n = m_n3Dpts;
	int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
	double* ptr1;
	Matrix3d matV, MatInv;

	//compute V inverse matrix using Eigen that has better performance than Lapack
	for (i = 0; i < n; ++i)
	{
		ptr1 = V + i * Vsz; // set ptr1 to point to V_i
		matV << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];
		MatInv = matV.inverse();
		ptr1[0] = MatInv(0, 0);
		ptr1[4] = MatInv(1, 1);
		ptr1[8] = MatInv(2, 2);
		ptr1[1] = ptr1[3] = MatInv(0, 1);
		ptr1[2] = ptr1[6] = MatInv(0, 2);
		ptr1[5] = ptr1[7] = MatInv(1, 2);
	}
}

double IBA::ba_computeInitialmu(double* U, double* V, sba_crsm& Uidxij, double tau, int nvars)
{
	int i, j;
	int pos, m = m_ncams, n = m_n3Dpts, cnp = 6, pnp = 3, Usz = 36, Vsz = 9;
	double tmp = 0;
	double* ptr1, * ptr2;
	double mu;

	double* diagUV = (double*)malloc(nvars * sizeof(double));

	double* diagU = diagUV;
	double* diagV = diagUV + m * cnp;

	for (j = 0; j < m; ++j)
	{
		pos = sba_crsm_elmidx(&Uidxij, j, j);
		ptr1 = U + pos * Usz;
		ptr2 = diagU + j * cnp;
		for (i = 0; i < cnp; ++i)
			ptr2[i] = ptr1[i * cnp + i];
	}
	for (i = 0; i < n; ++i)
	{
		ptr1 = V + i * Vsz; // set ptr1 to point to V_i
		ptr2 = diagV + i * pnp; // set ptr2 to point to diagV_i
		for (j = 0; j < pnp; ++j)
			ptr2[j] = ptr1[j * pnp + j];
	}

	/* find max diagonal element */
	for (i = 0, tmp = DBL_MIN; i < m * cnp; ++i)
		if (diagUV[i] > tmp)
			tmp = diagUV[i];
	for (i = m * cnp; i < nvars; ++i) /* tmp is not re-initialized! */
		if (diagUV[i] > tmp)
			tmp = diagUV[i];

	mu = m_Tau * tmp;

	free(diagUV);
	diagU = diagV = NULL;

	return mu;
}
void IBA::ba_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams)
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
				fprintf(stderr, "BA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles(fp, ptr1, 2); //第二次：读取第一个投影点的像点坐标 ptr1[0] ptr1[1]

			//也就是1个id+2个坐标
			if (n != 3)//一般不会发生
			{
				fprintf(stderr, "BA:reading image projections wrong!\n");
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

void IBA::ba_solveFeatures(double* W, double* IV, double* ea, double* eb, double* dpa, double* dpb)
{
	int i, j, ii, jj, pos, numfea;
	int nP1, cnp = 6, pnp = 3;
	double* ptr1, * ptr2, * ptr3, * ptr4, * ptr5;
	double sum, eb2[6];
	pos = 0;
	for (i = 0; i < m_n3Dpts; i++)
	{
		ptr1 = eb + i * 3;//指向三维点误差向量
		ptr2 = IV + i * 3 * 3;//指向V矩阵
		ptr5 = dpb + i * 3;//指向三维点未知参数
		memset(eb2, 0, sizeof(double) * cnp);
		numfea = m_archor[i * 3];

		for (j = 0; j < numfea; j++)
		{
			nP1 = m_photo[pos];//视图编号
			ptr3 = W + pos * cnp * 3;//指向W矩阵
			ptr4 = dpa + nP1 * cnp;//指向相机未知参数
			//Wta
			for (ii = 0; ii < pnp; ++ii)
			{
				for (jj = 0, sum = 0; jj < cnp; ++jj)
					sum += ptr3[jj * 3 + ii] * ptr4[jj];
				eb2[ii] += sum;
			}
			pos++;
		}

		//V*(eb-Wta)
		for (ii = 0; ii < pnp; ++ii)
		{
			for (jj = 0, sum = 0; jj < pnp; jj++)
				sum += ptr2[ii * 3 + jj] * (ptr1[jj] - eb2[jj]);
			ptr5[ii] = sum;
		}
	}
}
//高效求解GN中出现的稀疏线性系统
//使用Cholmod库进行矩阵（Cholesky）分解和求解
//Ap是列指针，Aii是行索引，是CSC格式
bool IBA::ba_solveCholmodGN(int* Ap, int* Aii, bool init, bool ordering)
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
				blockPermutation.resize(m_ncams - 1);//用于存储每个块的排列顺序

			// prepare AMD call via CHOLMOD
			cholmod_sparse auxCholmodSparse;//构造cholmod_sparse类型的辅助矩阵auxcholmodSparse，将稀疏矩阵元数据（如Ap和Aii）映射到Cholmod格式
			auxCholmodSparse.nzmax = m_nS;
			auxCholmodSparse.nrow = auxCholmodSparse.ncol = m - 1;//跳过第1帧
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
			if (!amdStatus) {
				return false;
			}

			// blow up the permutation to the scalar matrix
			//将块排列扩展为标量排列
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_cholSparseS->ncol);
			size_t scalarIdx = 0;

			int a = 0;
			for (i = 0; i < m_ncams - 1; ++i)
			{
				const int& pp = blockPermutation(i);
				int base = (pp == 0) ? 0 : pp * 6 - 1;
				int nCols = (pp == 0) ? 5 : 6;
				for (j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;
			}
			assert(scalarIdx == m_cholSparseS->ncol);//m_cholSparseS->ncol=90*6-7=533

			// apply the ordering
			m_cS.nmethods = 1;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;//设置CHOLMOD使用给定的排列顺序
			m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);//调用cholmod_analyze_p进行符号分析
		}
		init = true;//完成初始化
	}

	//Cholmod package for solving sparse linear equation              
	cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS); //对稀疏矩阵进行数值分解
	m_cholSparseR = cholmod_solve(CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS);//求解线性方程组

	return true;
}
//learning this skill from G2O
bool IBA::ba_solveCholmodLM(int* Ap, int* Aii, bool init, bool ordering)
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
				blockPermutation.resize(2 * m_ncams);

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
			if (!amdStatus)
				return false;


			// blow up the permutation to the scalar matrix
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_cholSparseS->ncol);
			if (scalarPermutation.size() < (int)m_cholSparseS->ncol)
				scalarPermutation.resize(2 * m_cholSparseS->ncol);
			size_t scalarIdx = 0;

			for (i = 0; i < m_ncams; ++i)
			{
				const int& pp = blockPermutation(i);
				int base = (pp == 0) ? 0 : pp * 6;
				int nCols = 6;

				for (j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;

			}
			assert(scalarIdx == m_cholSparseS->ncol);

			// apply the ordering
			m_cS.nmethods = 1;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;
			m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);
		}
	}

	cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS);
	m_cholSparseR = cholmod_solve(CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS);

	return true;
}
void IBA::ba_constructCSSLM(int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init)
{
	int ii, jj, jjj, k;
	int pos1, m = m_ncams;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double* ptr5;
	int nZ = 0;

	Sx = (double*)m_cholSparseS->x;
	if (!init)
	{
		for (ii = 0; ii < m; ii++)  //colum
		{
			for (k = 0; k < 6; k++)
			{
				*Sp = nZ;
				for (jj = 0; jj <= ii; jj++)	//row
				{
					if (m_smask[jj * m + ii] == 1)
					{
						pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
						ptr5 = S + pos1 * 36;

						if (ii == jj)
						{
							for (jjj = 0; jjj <= k; jjj++)
							{
								*Si++ = jj * 6 + jjj;
								*Sx++ = ptr5[jjj * 6 + k];
								nZ++;
							}
						}
						else
						{
							for (jjj = 0; jjj < 6; jjj++)
							{
								*Si++ = jj * 6 + jjj;
								*Sx++ = ptr5[jjj * 6 + k];
								nZ++;
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp = nZ;
	}
	else
	{
		for (ii = 0; ii < m; ii++)  //colum
		{
			for (k = 0; k < 6; k++)
			{
				for (jj = 0; jj <= ii; jj++)	//row
				{
					if (m_smask[jj * m + ii] == 1)
					{
						pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
						ptr5 = S + pos1 * 36;

						if (ii == jj)
						{
							for (jjj = 0; jjj <= k; jjj++)
								*Sx++ = ptr5[jjj * 6 + k];
						}
						else
						{
							for (jjj = 0; jjj < 6; jjj++)
								*Sx++ = ptr5[jjj * 6 + k];
						}
					}
				}
			}
		}
	}
}
void IBA::ba_constructCSSGN(int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init, int nft)
{
	int ii, jj, jjj, k;
	int pos1, m = m_ncams;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double* ptr5;
	int nZ = 0;
	Sx = (double*)m_cholSparseS->x;//存储稀疏矩阵中所有非零元素的实际值
	//printf("\n\n\n\n");
	if (!init)
	{
		for (ii = 1; ii < m; ii++)  //column，第0个视图是参考帧，通常被固定，不作为优化变量
		{
			for (k = 0; k < 6; k++)
			{
				*Sp = nZ;//表示稀疏矩阵中每一列的起始位置（非零元素的索引）,0,1,3,6,
				//printf("%d ", *Sp);
				if ((ii * 6 + k) == (9 + nft))//k=3代表Xc，4代表Yc，5代表Zc。nft=0代表固定Xc，1代表固定Yc，2代表固定Zc
					continue;//固定第2帧的Xc作为尺度约束

				for (jj = 1; jj <= ii; jj++)	//row
				{
					if ((m_smask[jj * m + ii] == 1))
					{
						pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
						//printf("%d ", pos1);
						ptr5 = S + pos1 * 36;//指向S矩阵对应（jj,ii）的子矩阵数据

						if (ii == jj)//对角线块
						{
							for (jjj = 0; jjj <= k; jjj++)
							{
								if ((jj * 6 + jjj) != (9 + nft))//除了jj=1，jjj=3（代表第2个相机的Xc固定）外，
								{
									if (jj * 6 + jjj < 9 + nft)//第2个相机的姿态角
										*Si++ = jj * 6 + jjj - 6;
									else//第2个相机的Yc，Zc
										*Si++ = jj * 6 + jjj - 7;

									*Sx++ = ptr5[jjj * 6 + k];
									//printf("%d ", ptr5[jjj * 6 + k]);
									nZ++;
								}
							}
						}
						else
						{
							for (jjj = 0; jjj < 6; jjj++)
							{
								if ((jj * 6 + jjj) != (9 + nft))
								{
									if (jj * 6 + jjj < 9 + nft)
										*Si++ = jj * 6 + jjj - 6;
									else
										*Si++ = jj * 6 + jjj - 7;

									*Sx++ = ptr5[jjj * 6 + k];
									nZ++;
								}
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp = nZ;
	}
	else
	{
		for (ii = 1; ii < m; ii++)  //column
		{
			for (k = 0; k < 6; k++)
			{
				if ((ii * 6 + k) == (9 + nft))
					continue;

				for (jj = 1; jj <= ii; jj++)	//row
				{
					if ((m_smask[jj * m + ii] == 1))
					{
						pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
						ptr5 = S + pos1 * 36;

						if (ii == jj)
						{
							for (jjj = 0; jjj <= k; jjj++)
							{
								if ((jj * 6 + jjj) != (9 + nft))
									*Sx++ = ptr5[jjj * 6 + k];
							}
						}
						else
						{
							for (jjj = 0; jjj < 6; jjj++)
							{
								if ((jj * 6 + jjj) != (9 + nft))
									*Sx++ = ptr5[jjj * 6 + k];
							}
						}
					}
				}
			}
		}
	}
}
void IBA::ba_constructAuxCSSLM(int* Ap, int* Aii)
{
	int* Cp = Ap;
	int* Ci = Aii;
	int ii, jj;
	int m = m_ncams, nZ = 0;
	for (ii = 0; ii < m; ii++)
	{
		*Cp = nZ;
		for (jj = 0; jj <= ii; jj++)
		{
			if (m_smask[jj * m + ii] == 1)
			{
				*Ci++ = jj;
				nZ++;
			}
		}
		Cp++;
	}
	*Cp = nZ;
}
//构造稀疏矩阵的索引结构，从1行1列开始，存储每列的起始索引,存储非零元素的行索引
void IBA::ba_constructAuxCSSGN(int* Ap, int* Aii)
{
	//从0行0列开始是因为固定第一帧
	int* Cp = Ap;//存储每列的起始索引
	int* Ci = Aii;//存储非零元素的行索引
	int ii, jj;
	int m = m_ncams, nZ = 0;
	//printf("\n\n\n");
	for (ii = 1; ii < m; ii++) //列
	{
		*Cp = nZ;
		for (jj = 1; jj <= ii; jj++)//行
		{
			if (m_smask[jj * m + ii] == 1)//按列遍历S矩阵（对称矩阵）的上三角
			{
				*Ci++ = jj - 1;//不能从0行开始
				//printf("%d ", jj-1);
				nZ++;
			}
		}
		Cp++;
	}
	*Cp = nZ;

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
void IBA::ba_constructSLM(double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij, double mu)
{
	int i, j, ii, jj, k, l;
	int m = m_ncams, cnp = 6, pnp = 3, Usz = 36, ebsz = 3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double* ptr1, * ptr2, * ptr3, * ptr4, * ptr5, * ptrS, * ptrE;
	double WV[6 * 3], sum;

	//Copy U matrix to S matrix 
	pos = 0;
	for (i = 0; i < m; i++) for (j = 0; j < m; j++)
	{
		if (m_umask[i * m + j] == 1)// save upper triangle for diagonal element S
		{
			pos1 = sba_crsm_elmidx(&Sidxij, i, j);
			ptr2 = S + pos1 * 36;
			if (i == j)
			{
				ptr1 = U + pos * Usz;
				for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
				{
					ptr2[ii] = ptr1[ii * cnp + ii] + mu;
					for (jj = ii + 1; jj < cnp; ++jj)
						ptr2[jj] = ptr1[ii * cnp + jj];
				}
				pos++;
			}
			else
			{
				ptr1 = U + pos * Usz;
				for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
					for (jj = 0; jj < cnp; ++jj)
						ptr2[jj] = ptr1[ii * cnp + jj];
				pos++;
			}
		}
	}

	for (i = 0; i < m * cnp; i++)
		E[i] = ea[i];

	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for (i = 0; i < m_n3Dpts; i++)
	{
		numfea = m_archor[i * 3];
		for (j = 0; j < numfea; j++)
		{
			nF1 = m_feature[pos];
			nP1 = m_photo[pos];
			memset(WV, 0, sizeof(double) * cnp * 3);

			ptr1 = W + pos * cnp * 3;
			ptr2 = V + nF1 * 3 * 3;
			ptrE = E + nP1 * cnp;


			//WV
			for (ii = 0; ii < cnp; ++ii)
			{
				ptr3 = ptr1 + ii * pnp;
				for (jj = 0; jj < pnp; ++jj)
				{
					for (k = 0, sum = 0.0; k <= jj; ++k)
						sum += ptr3[k] * ptr2[jj * pnp + k];
					for (; k < pnp; ++k)
						sum += ptr3[k] * ptr2[k * pnp + jj];
					for (k = 0, sum = 0.0; k < pnp; k++)
						sum += ptr3[k] * ptr2[jj * pnp + k];
					WV[ii * pnp + jj] = sum;
				}
			}

			for (k = j; k < numfea; k++)
			{
				nP2 = m_photo[pos + (k - j)];

				//W(V^-1)W^T
				ptr3 = W + (pos + (k - j)) * cnp * 3;
				//ptrS = S + (nP1*m*36) + nP2*cnp;

				if (nP1 == nP2)
				{
					pos1 = sba_crsm_elmidx(&Sidxij, nP1, nP2);
					ptrS = S + pos1 * 36;
					for (ii = 0; ii < cnp; ++ii, ptrS += 6)
					{
						ptr5 = WV + ii * pnp;
						for (jj = ii; jj < cnp; ++jj)
						{
							ptr4 = ptr3 + jj * pnp;

							for (l = 0, sum = 0.0; l < pnp; ++l)
								sum += ptr5[l] * ptr4[l];

							ptrS[jj] -= sum;
						}
					}
				}
				else
				{
					pos1 = sba_crsm_elmidx(&Sidxij, nP1, nP2);
					ptrS = S + pos1 * 36;
					for (ii = 0; ii < cnp; ++ii, ptrS += 6)
					{
						ptr5 = WV + ii * pnp;
						for (jj = 0; jj < cnp; ++jj)
						{
							ptr4 = ptr3 + jj * pnp;

							for (l = 0, sum = 0.0; l < pnp; ++l)
								sum += ptr5[l] * ptr4[l];

							ptrS[jj] -= sum;
						}
					}
				}
			}
			//-W^tb
			ptr5 = eb + nF1 * ebsz;
			for (ii = 0; ii < cnp; ++ii)
			{
				ptr4 = WV + ii * pnp;
				for (jj = 0, sum = 0.0; jj < pnp; ++jj)
					sum += ptr4[jj] * ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
				ptrE[ii] -= sum;
			}
			pos++;
		}
	}
}

void IBA::ba_constructSGN(double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij)
{
	int i, j, ii, jj, k, l;
	int m = m_ncams, cnp = 6, pnp = 3, Usz = 36, ebsz = 3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double* ptr1, * ptr2, * ptr3, * ptr4, * ptr5, * ptrS, * ptrE;
	double WV[6 * 3], sum;

	//Copy U matrix to S matrix 
	pos = 0;
	for (i = 0; i < m; i++) for (j = 0; j < m; j++)
	{
		if (m_umask[i * m + j] == 1)// save upper triangle for diagonal element S
		{
			pos1 = sba_crsm_elmidx(&Sidxij, i, j);
			//printf("%d ", pos1);
			ptr2 = S + pos1 * 36;
			if (i == j)
			{
				ptr1 = U + pos * Usz;
				for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
				{
					ptr2[ii] = ptr1[ii * cnp + ii];
					for (jj = ii + 1; jj < cnp; ++jj)
						ptr2[jj] = ptr1[ii * cnp + jj];
				}
				pos++;
			}
			else
			{
				ptr1 = U + pos * Usz;
				for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
					for (jj = 0; jj < cnp; ++jj)
						ptr2[jj] = ptr1[ii * cnp + jj];
				pos++;
			}
		}
	}

	for (i = 0; i < m * cnp; i++)
		E[i] = ea[i];//外方位元素

	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for (i = 0; i < m_n3Dpts; i++)
	{
		numfea = m_archor[i * 3];
		for (j = 0; j < numfea; j++)
		{
			nF1 = m_feature[pos];
			nP1 = m_photo[pos];
			memset(WV, 0, sizeof(double) * cnp * 3);

			ptr1 = W + pos * cnp * 3;
			ptr2 = V + nF1 * 3 * 3;
			ptrE = E + nP1 * cnp;

			//W(V^-1)，V已经变成了V^-1
			for (ii = 0; ii < cnp; ++ii)
			{
				ptr3 = ptr1 + ii * pnp;//W矩阵第ii行
				for (jj = 0; jj < pnp; ++jj)
				{
					for (k = 0, sum = 0.0; k <= jj; ++k)
						sum += ptr3[k] * ptr2[jj * pnp + k];
					for (; k < pnp; ++k)
						sum += ptr3[k] * ptr2[k * pnp + jj];
					for (k = 0, sum = 0.0; k < pnp; k++)//(V^-1)是对称矩阵
						sum += ptr3[k] * ptr2[jj * pnp + k];//此处的sum与上面的sum相等
					WV[ii * pnp + jj] = sum;
				}
			}

			for (k = j; k < numfea; k++)
			{
				nP2 = m_photo[pos + (k - j)];

				//W(V^-1)W^T
				ptr3 = W + (pos + (k - j)) * cnp * 3;
				pos1 = sba_crsm_elmidx(&Sidxij, nP1, nP2);
				ptrS = S + pos1 * 36;

				if (nP1 == nP2)
				{
					for (ii = 0; ii < cnp; ++ii, ptrS += 6)
					{
						ptr5 = WV + ii * pnp;	//WV矩阵第ii行								
						for (jj = ii; jj < cnp; ++jj)
						{
							ptr4 = ptr3 + jj * pnp;//W矩阵第jj行

							for (l = 0, sum = 0.0; l < pnp; ++l)
								sum += ptr5[l] * ptr4[l]; //-W(V^-1)W^T

							ptrS[jj] -= sum; //S-W(V^-1)W^T
						}
					}
				}
				else
				{
					for (ii = 0; ii < cnp; ++ii, ptrS += 6)
					{
						ptr5 = WV + ii * pnp;
						for (jj = 0; jj < cnp; ++jj)
						{
							ptr4 = ptr3 + jj * pnp;

							for (l = 0, sum = 0.0; l < pnp; ++l)
								sum += ptr5[l] * ptr4[l];

							ptrS[jj] -= sum; //方程左侧
						}
					}
				}
			}
			//-W^tb 方程右侧
			ptr5 = eb + nF1 * ebsz;//三维点
			for (ii = 0; ii < cnp; ++ii)
			{
				ptr4 = WV + ii * pnp;//W(V^-1)第ii行
				for (jj = 0, sum = 0.0; jj < pnp; ++jj)
					sum += ptr4[jj] * ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];   W(V^-1)eb
				ptrE[ii] -= sum; //ea-W(V^-1)eb
			}
			pos++;
		}
	}
}