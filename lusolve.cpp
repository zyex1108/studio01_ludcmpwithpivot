/*\file  lusolve.cpp
*
*   Solve linear equations by LU decomposition
*/
//#include <cstdio>
#include "cstdio"
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include <string>
#include <errno.h>
using namespace std;

char text[128];
void print_vector(double *v, int n);


void lusolve(Matrix &mx, FILE *in);
int ludcmp(Matrix &a, int order[]);
void solvlu(const Matrix &a, const Vec &b, Vec &x, const int order[]);

int main(int argc, char *argv[])
{
    Matrix mx;
    FILE *in;
    char *filename;

    if (argc>1)
    {
        filename = argv[1];
        printf("opening: %s\n",filename);
    }
    else
    {
        printf("This program requires a data file.\n");
        printf("Enter filename: ");
// 	gets_s(text,72);
//		gets(text,72);//by depengchen@126.com
        gets(text);
/*gets是从标准输入设备读字符串函数。
函数原型：char * gets ( char * str );
功能为：从stdin流中读取字符串，直至接受到换行符或EOF时停止，并
将读取的结果存放在buffer指针所指向的字符数组中。换行符不作为读取串的内容，
读取的换行符被转换为‘\0’空字符，并由此来结束字符串。
*/



        if (*text==0 || *text=='\n')
            return -1;

        filename = text;
    }

    //errno_t nerr = fopen_s(&in, filename,"rt");
    //////////////////////////////////////////////////
    //error_t nerr = fopen(filename,"r");r  readonly,t-text
/* 函数原型：FILE * fopen(const char * path,const char * mode);
*  返回值：文件顺利打开后，指向该流的文件指针就会被返回。如果文件打开失败则返回NULL，并把错误代码存在errno中。
*  一般而言，打开文件后会做一些文件读取或写入的动作，若打开文件失败，接下来的读写动作也无法顺利进行，
*  所以一般在fopen()后作错误判断及处理。
*/
    in=fopen(filename,"rt");
    if (in==NULL)
    {
        perror(filename);
        return -1;
    }

    string title;
    //当文本文件之矩阵数据读取成功的是和，返回整数0.


    while (read_matrix(mx,title,in)==0)
    {
        //printf("\n%s\n",title); 打印读出eqn5.dat文件中的文本注释。。。。
        printf("\n%s\n",title.c_str());
        //打印矩阵的内容，通过调用对象mx 的成员函数print（）
        printf("\nMatrix A\n");
        mx.print();

////////////////////////////////////////////
        lusolve(mx, in);



    }
    return 0;
}



void lusolve(Matrix &mx, FILE *in)
{
    int i, j, k, n, flag, index, nrhs;
    int *ipvt;
    double sum;
    Vec b, x, save;
    double det;
    /*
    *   Create matrix a as a copy of input matrix.
    */
    assert(mx.rows()==mx.cols()); // input must be square matrix
    n = mx.rows();

    Matrix a = mx;

//
//int *ivector(int n)
//{
//    int *iv;
//    iv = new int[n];
//    assert(iv);
//    return iv;
//}



    ipvt = ivector(n);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
    flag = ludcmp(a,ipvt);
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
    /* Calculate determinant */

    det = flag;

    for (i=0; i<n; i++)
        det *= a[i][i];

    printf("\ndeterminant: %g\n",det);
/////////////////////////////////////////////////////////////////////////////////


    printf("pivot vector: ");



    for (i=0; i<n; i++)
        printf("%3d",ipvt[i]);

    printf("\n");

    /* print LU Matrix */

    Matrix lower = a;

    Matrix upper = a;


    Matrix product(n, n);

    for (i=0; i<n; i++)
    {
        for (j=i+1; j<n; j++)
            lower[i][j] = 0.0;
        for (j=0; j<i; j++)
            upper[i][j] = 0.0;

        upper[i][i] = 1.0;//U对角线置 1.
    }
    printf("\nLower triangular Matrix\n");
    lower.print();

    printf("\nUpper triangular Matrix\n");
    upper.print();

    /* Multiply the lower and upper matrices */

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            sum = 0.0;
            for (k=0; k<n; k++)
                sum += lower[i][k] * upper[k][j];

            product[i][j] = sum;
        }
    }
    printf("\n\nProduct of L U\n\n");

    product.print();
    //应该得到矩阵A。。。。。。LU的乘积。。。
    /*
    *   Read in the right-hand side vectors.
    */
    //fscanf_s(in," %d",&nrhs);
    fscanf(in," %d",&nrhs);

    fgets(text,72,in);
    //创建2个矢量数组。。。。
    x.create(n);

    b.create(n);

    index = 0;

    while (nrhs--)
    {
        for (i=0; i<n; i++)
            // fscanf_s(in," %lf", &b[i]);
            fscanf(in," %lf", &b[i]);

        if (feof(in))
            break;

        printf("\nRight-hand-side number %d\n\n",++index);

        b.print();
        //输出矢量数组。。。b

        save = b;
//如果行列式值为0，则判断次矩阵是奇异的。。。。。。。。。。。
        if (fabs(det)<1e-10)
        {
            printf("\nCoefficient Matrix is singular.\n");
            continue;
        }



        //调用solvlu函数，求解x矢量数组。。。。。。。。。。。。。。。。。。。。。
        solvlu(a,b,x,ipvt);



        printf("\nSolution vector\n\n");
        x.print();


        //验证所求结果矢量与矩阵的乘积矢量是否相同，或者相差非常小。。。。。
        k = 0;
        for (i=0; i<n; i++)
        {
            sum = 0.0;
            for (j=0; j<n; j++)
                sum += mx[i][j]*x[j];

            if (fabs(sum-save[i])>1e-8)//对得到的结果进行判断是否相等。。。。如果等；如果大于容差，则意味不等，让标记k增1.
                k++;
        }
        printf("\nsolution %s.\n",(k? "does not check": "checks"));//如果k大于1，为真，则表示没有通过checks。



    }

//in matrix.cpp
//int *ivector(int n)
//{
//    int *iv;
//    iv = new int[n];
//    assert(iv);
//    return iv;
//}


    free(ipvt);
}

void print_vector(double *v, int n)
{
    while (n--)
        printf("%10.4f", *v++);

    printf("\n");
}
