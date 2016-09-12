/*\file ludcmp.cpp
*
*!  Finds solution to set of linear equations A x = b by LU decomposition.
*
*  Chapter 2, Programs 3-5, Fig. 2.8-2.10
*  Gerald/Wheatley, APPLIED NUMERICAL ANALYSIS (fourth edition)
*  Addison-Wesley, 1989
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"

int ludcmp(Matrix &a, int order[]);

void solvlu(const Matrix &a, const Vec &b, Vec &x, const int order[]);

static int pivot(Matrix &a, int order[], int jcol);

#define TINY 1e-20
//! finds LU decomposition of Matrix
/*!
*   The function ludcmp computes the lower L and upper U triangular
*   matrices equivalent to the A Matrix, such that L U = A.  These
*   matrices are returned in the space of A, in compact form.
*   The U Matrix has ones on its diagonal.  Partial pivoting is used
*   to give maximum valued elements on the diagonal of L.  The order of
*   the rows after pivoting is returned in the integer vector "order".
*   This should be used to reorder the right-hand-side vectors before
*   solving the system A x = b.
*
*
* \param       a     - n by n Matrix of coefficients
* \param       order - integer vector holding row order after pivoting.
*
*/
/*ludcmp计算下L和上U三角阵，LU=A。L和U存储在A中，用紧凑格式。U在其对角线上为1.
 *部分选主元用于给出位于L对角郑上的最大值。选主元后的行的顺序用整型矢量order返回。
 *在求解方程组Ax=b前，应该重排右手边的矢量。




*/
int ludcmp(Matrix &a, int order[])
{
    int i, j, k, n, nm1;
    int flag = 1;    /* changes sign with each row interchange */
    double sum, diag;

    n = a.rows();
    assert(a.cols()==n);

    /* establish initial ordering in order vector */
    //建立初始顺序。。。。。。。
    for (i=0; i<n; i++)
        order[i] = i;

    /* do pivoting for first column and check for singularity */
    //调用了选主元函数。。。。。。。。。如果进行了行交换，选主元，则pivot（）返回值为1.
    if (pivot(a,order,0))//先处理 第0列，即最后一个参数。。。。
        flag = -flag;
 //
    diag = 1.0/a[0][0];

///////////////////////////////////////////////////////////////////////////////
    for (i=1; i<n; i++)
        a[0][i] *= diag;
        // 对第0行，全部乘以系数
///////////////////////////////////////////////////////////////////////////////
    /*
    *  Now complete the computing of L and U elements.
    *  The general plan is to compute a column of L's, then
    *  call pivot to interchange rows, and then compute
    *  a row of U's.
    */
    //上面的思想很重要，先计算L，再选主元，决定是否交换，然后计算u行。

    nm1 = n - 1;
    for (j=1; j<nm1; j++)
    {
        /* column of L's */
        for (i=j; i<n; i++)
        {
            sum = 0.0;

            for (k=0; k<j; k++)
                sum += a[i][k]*a[k][j];

            a[i][j] -= sum;
        }
        /* pivot, and check for singularity */
        //调用莉选主元函数。。。。。。。。。对行交换标记进行标识。。。。有行交换，则返回1，否则返回值为0；
        if (pivot(a,order,j))
            flag = -flag;

        /* row of U's */

        diag = 1.0/a[j][j];

        for (k=j+1; k<n; k++)
        {
            sum = 0.0;

            for (i=0; i<j; i++)
                sum += a[j][i]*a[i][k];

            a[j][k] = (a[j][k]-sum)*diag;
        }
    }

    /* still need to get last element in L Matrix */
    //仍然需要得到L矩阵的最后一个元素。。。。。

    sum = 0.0;

    for (k=0; k<nm1; k++)
        sum += a[nm1][k]*a[k][nm1];

    a[nm1][nm1] -= sum;

    return flag;
}

//!  Find pivot element  部分选主元。。。
/*!
*   The function pivot finds the largest element for a pivot in "jcol"
*   of Matrix "a", performs interchanges of the appropriate
*   rows in "a", and also interchanges the corresponding elements in
*   the order vector.
*
*
*  \param     a      -  n by n Matrix of coefficients
*  \param     order  - integer vector to hold row ordering
*  \param     jcol   - column of "a" being searched for pivot element
*
*/
//寻找矩阵a中jcol列的最大元素，且对数组order的元素进行交换。
int pivot(Matrix &a, int order[], int jcol)
{
    int i, ipvt,n;
    double big, anext;
    n = a.rows();

    /*
    *  Find biggest element on or below diagonal.
    *寻找再对角线或者对角线以下的最大单元。。。。
    *  This will be the pivot row.这将成为主行。。。。
    */

    ipvt = jcol;
    //先令矩阵的第ipvt行和ipvt列的绝对值为最大值。。。。
    big = fabs(a[ipvt][ipvt]);
     //对第ipvt行以下的各行进行循环，得到jcol列中的最大值
    for (i = ipvt+1; i<n; i++)
    {
        anext = fabs(a[i][jcol]);
        if (anext>big)
        {
            //令最大值等于big，且令i行的行号等于主元所在行号。。。。。。
            big = anext;
            ipvt = i;
        }
    }
    //断言最大值之绝对值大于最小容差，
    assert(fabs(big)>TINY); // otherwise Matrix is singular

    /*
    *   Interchange pivot row (ipvt) with current row (jcol).
    */
    //如果主元所在的行ipvt与jcol相同，则返回该函数被调用处。
    if (ipvt==jcol)
        return 0;
    //交换jcol行和ipvt行的各个元素。。。。
    a.swaprows(jcol,ipvt);

//////////////////////////////////
//对order数组的内容进行更新。。。。。。。。。。。。。。
    i = order[jcol];
    order[jcol] = order[ipvt];
    order[ipvt] = i;
/////////////////////////////////
    return 1;
}

//!  This function is used to find the solution to a system of equations,
/*!   A x = b, after LU decomposition of A has been found.
*    Within this routine, the elements of b are rearranged in the same way
*    that the rows of a were interchanged, using the order vector.
*    The solution is returned in x.
*   在本子程序中，b的元素重排后与矩阵a行交换后的顺序一致，用order矢量数组。。。
*
*  \param  a     - the LU decomposition of the original coefficient Matrix.
*  \param  b     - the vector of right-hand sides
*  \param  x     - the solution vector
*  \param  order - integer array of row order as arranged during pivoting
*
*/
void solvlu(const Matrix &a, const Vec &b, Vec &x, const int order[])
{
    int i,j,n;
    double sum;
    n = a.rows();

    /* rearrange the elements of the b vector. x is used to hold them. */

    for (i=0; i<n; i++)
    {
        j = order[i];
        x[i] = b[j];
    }

    /* do forward substitution, replacing x vector. */

    x[0] /= a[0][0];
    for (i=1; i<n; i++)
    {
        sum = 0.0;

        for (j=0; j<i; j++)
            sum += a[i][j]*x[j];

        x[i] = (x[i]-sum)/a[i][i];
    }

    /* now get the solution vector, x[n-1] is already done */

    for (i=n-2; i>=0; i--)
    {
        sum = 0.0;

        for (j=i+1; j<n; j++)
            sum += a[i][j] * x[j];

        x[i] -= sum;
    }
}
