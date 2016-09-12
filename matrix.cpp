/* matrix.cpp
*
*    Utility routines for allocating space for matrices and vectors,
*    printing matrices and vectors, and copying one matrix to another.
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "matrix.h"

void Matrix::create(int rows, int cols)
{
    //如果准被创建的矩阵，规模与默认构造函数规模相同，就返回，利用默认构造函数创建矩阵对象。
    if (rows==nrows && cols==ncols)
        return;
    //如果要重新创建对象。先删除指向每行的指针，在删除指向矩阵的指针。
    while (nrows)
        delete [] body[--nrows];
    if (body)
        delete body;

    nrows = rows;
    ncols = cols;
    body = new double * [nrows];
    assert(body);//判断分配内存是否成功，body如果指向NULL，则断言

    for (int i=0; i<nrows; i++)
    {
        body[i]= new double[ncols];
        //body是指针的指针。。。。。
        assert(body[i]);
    }
}

/*!
*   Copy Matrix src to the current Matrix
*/
Matrix &Matrix::copy(const Matrix &src)
{
    int i, j;
    double *up, *vp;
    create(src.rows(),src.cols());
    for (i=0; i<nrows; i++)
    {
        up = src[i];//源对象之笫i行之首元素之指针
        vp = body[i];//目标对象之第i行之首元素之指针

        for (j=0; j<ncols; j++)
            *vp++ = *up++;//深复制。。。。
    }
    return *this;
}
//析构函数
Matrix::~Matrix()
{
    while (nrows) delete [] body[--nrows];
    delete [] body;
}
//对矩阵赋值，全部为双精度浮点数v
Matrix &Matrix::set(double v)
{
    int i, j;
    double *p;
    for (i=0; i<nrows; i++)
    {
        p = body[i];

        for (j=0; j<ncols; j++)
            *p++ = v;
    }
    return *this;
}

//交换i和j两行。。。。。。

void Matrix::swaprows(int i,int j)
{
    double *p = body[i];
    body[i]=body[j];
    body[j]=p;
}


//转置矩阵
Matrix Matrix::transpose()
{
    int i, j;
    Matrix out(ncols,nrows);//定义一个临时矩阵，为存储转置服务。
    for (i=0; i<nrows; i++)
    {
        for (j=0; j<ncols; j++)
        {
            out[j][i] = body[i][j];
        }
    }
    return out;
}

//! Matrix multiplication两个矩阵相乘。。。。。

Matrix Matrix::operator *(const Matrix &mx) const
{
    int i, j, k;
    double sum;
    double *vp;
    assert(ncols == mx.rows()); // check for non-conforming Matrix multiplication

    Matrix out(nrows,mx.cols());

    for (i=0; i<nrows; i++)
    {
        vp = out[i];
        for (j=0; j<mx.cols(); j++)
        {
            sum = 0.0;
            for (k=0; k<ncols; k++)
                sum += body[i][k] * mx[k][j];
            *vp++ = sum;
        }
    }
    return out;
}

//矩阵与矢量的乘积。返回一个矩阵
Vec operator *(const Matrix &mx, const Vec &v)
{
    int i, j;
    int m = mx.rows();
    int n = mx.cols();
    double sum;
    assert(m>0 && n>0);
    assert(n==v.length());
    Vec vout(m);

    for (i=0; i<m; i++)
    {
        sum = 0.0;
        double *ptr = mx[i];

        for (j=0; j<n; j++)
            sum += (*ptr++) * v[j];
        vout[i] = sum;
    }
    return vout;
}

//矩阵与矩阵的乘积，返回一个矢量

Vec operator *(const Vec &v, const Matrix &mx)
{
    int i, j;
    int m = mx.rows();
    int n = mx.cols();
    double sum;
    assert(m>0 && n>0);
    assert(m==v.length());
    Vec vout(n);
    for (j=0; j<n; j++)
    {
        sum = 0.0;
        for (i=0; i<m; i++) sum += v[i] * mx[i][j];
        vout[j] = sum;
    }
    return vout;
}

//!    create an Identity Matrix，创建一个单位矩阵。。。。
Matrix identity(int size)
{
    int i, j;
    double *vp;
    Matrix out(size,size);
    for (i=0; i<size; i++)
    {
        vp = out[i];
        for (j=0; j<size; j++)
            vp[j]=0.0;

        vp[i] = 1.0;//仅对角线元素为1
    }
    return out;
}


//! create a diagonal Matrix,创建一个对角阵
/*!
 *  @param v is an array of diagonal elements
 *  @param size is the length of the array
 *
 */
Matrix diagonal(double *v, int size)
{
    int i, j;
    double *vp;
    Matrix out(size,size);
    for (i=0; i<size; i++)
    {
        vp = out[i];
        for (j=0; j<size; j++)
            vp[j]=0.0;

        vp[i] = v[i];
    }
    return out;
}

//! read Matrix from text file,从文本文件中读取矩阵
/*! @param mx destination Matrix，mx为码表矩阵
 *  @param title c-string comment (remainder of header line)，c字串注释，保存在头行
 *  @param in input file输入文件
 *
*  The first two entries are number of rows and columns.\n
*  The rest of the line is a title.\n
*  The body of the Matrix follows.\n
*/

int read_matrix(Matrix &mx, string &title, FILE *in)
{
    int i, j;
    int nrows, ncols;
    double *v;
    double vin;
    char *cp;
    char buf[256];

    //fscanf_s(in," %d %d",&nrows,&ncols);
    fscanf(in," %d %d",&nrows,&ncols);
    if (feof(in))
        return -1;

    fgets(buf,255,in);

    if (feof(in))
    return -1;

    cp = buf;
    while (*cp)
    {
        if (*cp=='\n') *cp = 0;
        else cp++;
    }
    //读完tittle，如果读莉255个字符，还没有到\n则，继续，知道字符串尾部

    title = buf;
    mx.create(nrows,ncols);
    for (j=0; j<nrows; j++)
    {
        v  = mx[j];
        for (i=0; i<ncols; i++)
        {
            fscanf(in," %lf",&vin);
            //fscanf_s(in," %lf",&vin);
            *v++ = vin;
        }
        if (feof(in))
        {
            printf("\nerror reading %s\n",title.c_str());
            printf("Matrix has %d rows and %d columns\n",
                   mx.rows(),mx.cols());
            printf("Unexpected EOF reading row %d",j+1);
            return -1;
        }
    }
    return 0;
}




//class Vec类成员函数的实现。。。。。。。
//! create memory for elements of vector
void Vec::create(int size)
{
    assert(size);
    if (size == n)
        return;

    if (body != 0)
        delete [] body;

    body = new double[size];
    assert(body!=0);
    n = size;//用size初始化更新Vec的数据成员n。。。。。

    //将已开辟的空间body,用（size_t）(n*sizeof(double))个字符 0 来填充。
    //memset((void *)body, 0, (size_t)(n*sizeof(double));
}
//! copy constructor
Vec &Vec::copy(const Vec &src)
{
    int size = src.length();

    create(size);

    for (int i=0; i<n; i++)
        body[i] = src[i];

    return *this;
}

double Vec::max()
{
    double vmax = body[0];

    for (int i=1; i<n; i++)
        if (body[i]>vmax)
            vmax = body[i];

    return vmax;
}

double Vec::min()
{
    double vmin = body[0];

    for (int i=1; i<n; i++)
        if (body[i]<vmin)
        vmin = body[i];

    return vmin;
}

//得到二范数，即矢量的几何长度--模
double Vec::norm()
{
    double sum = 0.0;

    for (int i=0; i<n; i++)
        sum += body[i]*body[i];

    return sqrt(sum);
}

//矢量归一化

void Vec::normalize()
{
    double vnorm = norm();

    for (int i=0; i<n; i++)
        body[i] /= vnorm;
}

//对矢量中的元素，施加吗中被fct指针指向的函数，即用该函数算法对元素进行处理，再返回处理后的函数。
//相当于标准模板库里面的算法，施加于容器Vec，很有意思的编程思想。
Vec& Vec::apply(V_FCT_PTR fct)
{
    for (int i=0; i<n; i++)
        body[i] = (*fct)(body[i]);
        //利用fct指向的函数，对矢量数组为参数进行操作，再返回矢量数组。

    return *this;
}

//! returns a + b
Vec operator +(const Vec &a, const Vec &b)
{
    int n = a.length();
    assert(b.length()==n);
    Vec sum(n);
    for (int i=0; i<n; i++)
        sum[i] = a[i] + b[i];
    return sum;
}

//! returns a - b
Vec operator -(const Vec &a, const Vec &b)
{
    int n = a.length();
    assert(b.length()==n);
    Vec diff(n);
    for (int i=0; i<n; i++)
        diff[i] = a[i] - b[i];
    return diff;
}

//! returns scalar multiplied by a vector
Vec operator *(double c, const Vec &a)
{
    int n = a.length();
    Vec vout(n);
    for (int i=0; i<n; i++)
        vout[i] = c*a[i];
    return vout;
}

//! returns vector mutliplied by a scalar
Vec operator *(const Vec &a, double c)
{
    int n = a.length();
    Vec vout(n);
    for (int i=0; i<n; i++)
        vout[i] = c*a[i];
    return vout;
}

//! returns vector divided by a scalar
Vec operator /(const Vec &a, double c)
{
    int n = a.length();
    Vec vout(n);
    for (int i=0; i<n; i++)
        vout[i] = a[i]/c;
    return vout;
}

double *vector(int n)
{
    double *v;
    v = new double[n];
    assert(v);
    return v;
}


int *ivector(int n)
{
    int *iv;
    iv = new int[n];
    assert(iv);
    return iv;
}


void Matrix::print() const
{
    int i, j;
    const double tiny = 1e-13;
    double *vp;
    double v;
    for (j=0; j<nrows; j++)
    {
        vp = body[j];
        for (i=0; i<ncols; i++)
        {
            v = *vp++;
            if (fabs(v)<tiny)
                v = 0.0;

            printf("% 10.4g", v);
        }
        printf("\n");
    }
}

void Vec::print() const
{
    double *v = body;
    int size = n;
    while (size--)
    {
        printf(" %10.4g", *v++);
    }
    printf("\n");
}

//! copy one Matrix to another
/*!
 * @param dst destination Matrix
 * @param src source Matrix
 * if destination is smaller than the source, the source is truncated.\n
 * if destination is larger than the source, the remainder is
 * zero-filled\n
 */

void copy_matrix(Matrix &dst, Matrix &src)
{
    int i, j;
    int nrows, ncols;
    double *s, *t;;
    nrows = (dst.rows()<src.rows()? dst.rows(): src.rows());
    ncols = (dst.cols()<src.cols()? dst.cols(): src.cols());
    for (j=0; j<nrows; j++)
    {
        s = src[j];
        t = dst[j];
        for (i=0; i<ncols; i++)
            *t++ = *s++;
    }
    // Append zero-filled rows to destination
    for (j=src.rows(); j<nrows; j++)
    {
        t = dst[j];
        for (i=0; i<ncols; i++)
            t[i] = 0.0;
    }
    // Append zero-filled columns to destination
    if (src.cols() < ncols)
    {
        for (j=0; j<dst.rows(); j++)
        {
            t = dst[j];
            for (i=src.cols(); i<ncols; i++)
                t[i] = 0.0;
        }
    }
}

// overloaded insertion operator <<
std::ostream& operator << (std::ostream& s, const Vec& v)
{
   char buf[16];
   // char* buf= new char[16];
    int n = v.length();
    if (n > 0)
    {
        //int old_precision = s.precision() ; // get current precision
        s << "[" ;
        for (int i=0; i<n; i++)
        {
            //sprintf_s(buf,15,"%.4g",v[i]);
            ///////////////////////////////////////////////
            sprintf(buf,"%.4g",15,v[i]);//by depengchen@126.com
            //s << setprecision (2)
            //<< setiosflags (ios_base::showpoint|ios_base::fixed)
            s << buf ;
            i!=n-1 ? s << ", " : s << "]" ;
        }
        //s << endl ;
        // reset precision and ios flags
        //s.precision (old_precision) ;
        //std::resetiosflags (ios::showpoint|std::ios::fixed) ;
    }
    return s ;
}
