#pragma warning(disable:4996)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "juvf.h"
#include "dcf.h"
#define tmin 0.6995	//周期起始值
#define tmax 0.7005	//周期结束值
#define tstep 0.000001	//周期步长

int main()
{
	char olcname[50];
	printf("请输入光变曲线文件路径：");
	scanf("%s", olcname);
	FILE* olc;
	olc = fopen(olcname, "r");

	int np = 0;
	int nlc;	//光变曲线数量
	(void)fscanf(olc, "%d", &nlc);
	int** nlp;	//数据点数量
	nlp = (int**)malloc(nlc * sizeof(int));
	for (int i = 0; i < nlc; i++)
	{
		nlp[i] = (int*)malloc(2 * sizeof(int));
	}
	double*** lp;	//数据点
	lp = (double***)malloc(nlc * sizeof(double));
	double*** lpc;
	lpc = (double***)malloc(nlc * sizeof(double));

	for (int i = 0; i < nlc; i++)
	{
		(void)fscanf(olc, "%d %d", &nlp[i][0], &nlp[i][1]);
		lp[i] = (double**)malloc(nlp[i][0] * sizeof(double));
		lpc[i] = (double**)malloc(nlp[i][0] * sizeof(double));
		for (int j = 0; j < nlp[i][0]; j++)
		{
			lp[i][j] = (double*)malloc(8 * sizeof(double));
			lpc[i][j] = (double*)malloc(8 * sizeof(double));
			for (int k = 0; k < 8; k++)
			{
				(void)fscanf(olc, "%lf", &lp[i][j][k]);
				lpc[i][j][k] = lp[i][j][k];
			}
			np++;
		}
	}

	int M = 0;
	double dcf_r = 0;
	double d = 0;
	double temp = -1e26;
	double T = 0;
	for (double t = tmin; t < tmax; t = t + tstep)
	{
		dcf_r = 0;
		M = 0;
		d = 0;
		for (int i = 0; i < nlc; i++)
		{
			for (int j = 0; j < nlp[i][0]; j++)
			{
				lpc[i][j][0] = lp[i][j][0] - (double)((int)(lp[i][j][0] / t) * t);
			}
		}
		for (int i = 0; i < nlc - 1; i++)
		{
			for (int j = 1; j < nlc; j++)
			{
				if (min(lpc[j][nlp[j][0] - 1][0], lpc[i][nlp[i][0] - 1][0]) - max(lpc[j][0][0], lpc[i][0][0] > 0.01))
				{
					d = DCF(lpc[i], lpc[j], nlp[i][0], nlp[j][0], 0, t * 0.01);
					if (d != -13)
					{
						dcf_r += d;
						M++;
					}
				}
			}
		}
		dcf_r = dcf_r / M;
		if (temp < dcf_r)
		{
			temp = dcf_r;
			T = t;
		}
	}
	printf("%f %f\n", T, temp);

	return 0;
}