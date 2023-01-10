#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<string.h>
#include<iostream>  
#include <iomanip>  

using namespace std;
#define  N  4                //定义方阵的最大阶数为4

//函数的声明部分  
float MatDet(float* p, int n);                    //求矩阵的行列式  
float Creat_M(float* p, int m, int n, int k);    //求矩阵元素A(m, n)的代数余之式  
void print(float* p, int n);                    //输出矩阵n*n  
bool Gauss(float A[][N], float B[][N], int n);    //采用部分主元的高斯消去法求方阵A的逆矩阵B  



// 0 - 1
float t0 = 0, x0 = 30, v0 = 75;
float t1 = 5, x1 = 60, v1 = 10;
//[x01, v01, a01] = plan(t0, x0, v0, t1, x1, v1);
// 1 - 2
float t2 = 13, x2 = 80, v2 = 10;
//[x12, v12, a12] = plan(0, x1, v1, t2 - t1, x2, v2);
// 2 - 3
float t3 = 20, x3 = 10, v3 = 0;
//[x23, v23, a23] = plan(0, x2, v2, t3 - t2, x3, v3);


float  ts, start_x, start_v, te, end_x, end_v;

struct Node
{
	float x;
	float v;
	float a;
	struct Node* next;
};
typedef struct Node List;




//-----------------------------------------------  
//功能: 求矩阵(n*n)的行列式  
//入口参数: 矩阵的首地址，矩阵的行数  
//返回值: 矩阵的行列式值  
//----------------------------------------------  
float MatDet(float* p, int n)
{
	int r, c, m;
	int lop = 0;
	float result = 0;
	float mid = 1;

	if (n != 1)
	{
		lop = (n == 2) ? 1 : n;            //控制求和循环次数,若为2阶，则循环1次，否则为n次  
		for (m = 0; m < lop; m++)
		{
			mid = 1;            //顺序求和, 主对角线元素相乘之和  
			for (r = 0, c = m; r < n; r++, c++)
			{
				mid = mid * (*(p + r * n + c % n));
			}
			result += mid;
		}
		for (m = 0; m < lop; m++)
		{
			mid = 1;            //逆序相减, 减去次对角线元素乘积  
			for (r = 0, c = n - 1 - m + n; r < n; r++, c--)
			{
				mid = mid * (*(p + r * n + c % n));
			}
			result -= mid;
		}
	}
	else
		result = *p;
	return result;
}

//----------------------------------------------------------------------------  
//功能: 求k*k矩阵中元素A(m, n)的代数余之式  
//入口参数: k*k矩阵的首地址，矩阵元素A的下标m,n,矩阵行数k  
//返回值: k*k矩阵中元素A(m, n)的代数余之式  
//----------------------------------------------------------------------------  
float Creat_M(float* p, int m, int n, int k)
{
	int len;
	int i, j;
	float mid_result = 0;
	int sign = 1;
	float* p_creat, * p_mid;

	len = (k - 1) * (k - 1);            //k阶矩阵的代数余之式为k-1阶矩阵  
	p_creat = (float*)calloc(len, sizeof(float)); //分配内存单元  
	p_mid = p_creat;
	for (i = 0; i < k; i++)
	{
		for (j = 0; j < k; j++)
		{
			if (i != m && j != n) //将除第i行和第j列外的所有元素存储到以p_mid为首地址的内存单元  
			{
				*p_mid++ = *(p + i * k + j);
			}
		}
	}
	sign = (m + n) % 2 == 0 ? 1 : -1;    //代数余之式前面的正、负号  
	mid_result = (float)sign * MatDet(p_creat, k - 1);
	free(p_creat);
	return mid_result;
}

//-----------------------------------------------------  
//功能: 打印n*n矩阵  
//入口参数: n*n矩阵的首地址,矩阵的行数n  
//返回值: 无返回值  
//-----------------------------------------------------  
void print(float* p, int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		cout << setw(4);
		for (j = 0; j < n; j++)
		{
			cout << setiosflags(ios::right) << *p++ << setw(10);
		}
		cout << endl;
	}
}

//------------------------------------------------------------------  
//功能: 采用部分主元的高斯消去法求方阵A的逆矩阵B  
//入口参数: 输入方阵，输出方阵，方阵阶数  
//返回值: true or false  
//-------------------------------------------------------------------  
bool Gauss(float A[][N], float B[][N], int n)
{
	int i, j, k;
	float max, temp;
	float t[N][N];                //临时矩阵  
	//将A矩阵存放在临时矩阵t[n][n]中  
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			t[i][j] = A[i][j];
		}
	}
	//初始化B矩阵为单位阵  
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			B[i][j] = (i == j) ? (float)1 : 0;
		}
	}
	for (i = 0; i < n; i++)
	{
		//寻找主元  
		max = t[i][i];
		k = i;
		for (j = i + 1; j < n; j++)
		{
			if (fabs(t[j][i]) > fabs(max))
			{
				max = t[j][i];
				k = j;
			}
		}
		//如果主元所在行不是第i行，进行行交换  
		if (k != i)
		{
			for (j = 0; j < n; j++)
			{
				temp = t[i][j];
				t[i][j] = t[k][j];
				t[k][j] = temp;
				//B伴随交换  
				temp = B[i][j];
				B[i][j] = B[k][j];
				B[k][j] = temp;
			}
		}
		//判断主元是否为0, 若是,则矩阵A不是满秩矩阵,不存在逆矩阵  
		if (t[i][i] == 0)
		{
			cout << "There is no inverse matrix!";
			return false;
		}
		//消去A的第i列除去i行以外的各行元素  
		temp = t[i][i];
		for (j = 0; j < n; j++)
		{
			t[i][j] = t[i][j] / temp;        //主对角线上的元素变为1  
			B[i][j] = B[i][j] / temp;        //伴随计算  
		}
		for (j = 0; j < n; j++)        //第0行->第n行  
		{
			if (j != i)                //不是第i行  
			{
				temp = t[j][i];
				for (k = 0; k < n; k++)        //第j行元素 - i行元素*j列i行元素  
				{
					t[j][k] = t[j][k] - t[i][k] * temp;
					B[j][k] = B[j][k] - B[i][k] * temp;
				}
			}
		}
	}
	return true;
}




float x=0, v=0, a=0;
float  c0, c1, c2, c3;
float z[10][3];

List* CreatList(int n, float *zq )
{
	List* pre, * cur, * p;
	p = NULL;
	for (int j = 0; j < n; j++)
	{
		cur = (List*)malloc(sizeof(List));
		cur->next = NULL;
		if (p == NULL)
			p = cur;
		else
			pre->next = cur;
		cur->x = *zq, cur->v = *zq+1, cur->a = *zq+2;
		pre = cur;
	}
	return p;
}

float*  plan(int k,float ts, float start_x, float start_v, float te, float end_x, float end_v)
{
	float* p;
	float i, t;
	for ( i = ts;i<= te; i = i + 0.01)
	{
		t = i;

		x = x + (c0 + c1 * t + c2 * pow(t , 2) + c3 * pow(t , 3));
		v = v + (c1 + 2 * c2 * t + 3 * c3 * pow(t , 2));
		a = a + (2 * c2 + 6 * c3 * t);
	}
	
	//printf("%.2f %.2f %.2f\n", x, v, a);
	z[k][0] = x, z[k][1] = v, z[k][2] = a;
	p = *z;
	return p ;
}

void Ci(float a[4][4], float b[4][1], float c[4][1])
{
	int i, j, k;
	float sum = 0;
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
			sum += a[i][j] * b[j][0];
		c[i][0] = sum;
		sum = 0;
	}

}
void Pr(List* h)
{
	List* p = h;
	while (p != NULL)
	{
		printf("%.2f %.2f %.2f \n", p->x, p->v, p->a);
		p = p->next;
	}
}
int main()
{

	//int X[], V[], A[];
	List* p;
	float a[N][N], b[N][N];
	int i, j, n;


	printf("输入阶段数：\n");
	scanf("%d", &n);
	for(int k=0;k<n;k++)
	{
		printf("请输入 ts, start_x, start_v, te, end_x, end_v ：\n");
		scanf("%f%f%f%f%f%f", &ts, &start_x, &start_v, &te, &end_x, &end_v);

		float para[4][1] = { start_x, end_x, start_v, end_v };
		float c[4][1];
		float Tran[4][4] =
		{ 1, ts,  pow(ts ,2), pow(ts , 3),
		 1, te,  pow(te , 2), pow(te , 3),
		 0, 1,   2 * ts,  3 * pow(ts , 2),
		 0, 1,   2 * te, 3 * pow(te , 2) };



		//运用高斯消去法求该矩阵的逆矩阵并输出  
		if (Gauss(Tran, b, 4))
		{
			cout << "Tran方阵的逆矩阵为: \n";
			for (i = 0; i < 4; i++)
			{
				//cout << setw(4);
				for (j = 0; j < 4; j++)
				{
					printf("%.2f ", b[i][j]);
				}
				cout << endl;
			}
		}

		Ci(b, para, c);

		c0 = c[0][0];
		c1 = c[1][0];
		c2 = c[2][0];
		c3 = c[3][0];
	//printf("jieguo:\n");
	//for (i = 0; i < 4; i++)
	//	printf("%.2f \n", c[i][0]);


	 plan(k, ts, start_x, start_v, te, end_x, end_v);

     }
	printf("输出第i阶段的x  v  a  \n");
	for (i = 0; i < n; i++)
		printf("第%d阶段 :%.2f %.2f %.2f\n",i, z[i][0], z[i][1], z[i][2]);
	
	return 0;
}