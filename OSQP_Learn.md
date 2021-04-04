# OSQP

## 主要链接

[官网]( https://osqp.org/ )

[github]( https://github.com/oxfordcontrol/osqp )

[C example]( https://github.com/oxfordcontrol/osqp/tree/master/examples )

[各个语言求解示例地址]( https://osqp.org/docs/examples/setup-and-solve.html )

## 求解问题概述

OSQP求解以下形式的凸二次规划问题：
$$
\begin{split}\begin{array}{ll}  \mbox{minimize} & \frac{1}{2} x^T P x + q^T x \\  \mbox{subject to} & l \leq A x \leq u\end{array}\end{split}
$$
式中，$$ x\in\mathbf{R}^{n}$$是优化变量，目标函数由半正定矩阵$$ P\in\mathbf{S}^{n}_{+}$$和$$ q\in\mathbf{R}^{n}$$定义，线性约束由矩阵$$ A \in\mathbf{R}^{m \times n} $$和$$ l $$, $$ u $$定义，$$ l_{i} \in\ mathbf{R} \cup \{ -\infty\} $$和$$ u_{i} \in\ mathbf{R} \cup \{ +\infty\} $$ 对于所有的$$ i \in \{1,\ldots,m\} $$。



半正定矩阵一般来说都是对称矩阵，QSQP在使用C语言求解时，是将 $$ P $$ 和 $$ A $$ 都存储为csc_matrix( compressed sparse colmun matrix) 基于列的稀疏矩阵，需要注意的是，因为矩阵都为对称矩阵，在构造csc_matix时，仅存储上三角矩阵元素。

## C语言求解疑问

$$
\begin{split}\begin{array}{ll}
  \mbox{minimize} & \frac{1}{2} x^T \begin{bmatrix}4 & 1\\ 1 & 2 \end{bmatrix} x + \begin{bmatrix}1 \\ 1\end{bmatrix}^T x \\
  \mbox{subject to} & \begin{bmatrix}1 \\ 0 \\ 0\end{bmatrix} \leq \begin{bmatrix} 1 & 1\\ 1 & 0\\ 0 & 1\end{bmatrix} x \leq  \begin{bmatrix}1 \\ 0.7 \\ 0.7\end{bmatrix}
\end{array}\end{split}
$$



**C语言**

>```c
>#include "osqp.h"
>
>int main(int argc, char **argv) {
>    // Load problem data
>    c_float P_x[3] = {4.0, 1.0, 2.0, };
>    c_int P_nnz = 3;
>    c_int P_i[3] = {0, 0, 1, };
>    c_int P_p[3] = {0, 1, 3, };
>    c_float q[2] = {1.0, 1.0, };
>    c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };
>    c_int A_nnz = 4;
>    c_int A_i[4] = {0, 1, 0, 2, };
>    c_int A_p[3] = {0, 2, 4, };
>    c_float l[3] = {1.0, 0.0, 0.0, };
>    c_float u[3] = {1.0, 0.7, 0.7, };
>    c_int n = 2;
>    c_int m = 3;
>
>    // Exitflag
>    c_int exitflag = 0;
>
>    // Workspace structures
>    OSQPWorkspace *work;
>    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
>    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));
>
>    // Populate data
>    if (data) {
>        data->n = n;
>        data->m = m;
>        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
>        data->q = q;
>        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
>        data->l = l;
>        data->u = u;
>    }
>
>    // Define solver settings as default
>    if (settings) {
>        osqp_set_default_settings(settings);
>        settings->alpha = 1.0; // Change alpha parameter
>    }
>
>    // Setup workspace
>    exitflag = osqp_setup(&work, data, settings);
>
>    // Solve Problem
>    osqp_solve(work);
>
>    // Cleanup
>    if (data) {
>        if (data->A) c_free(data->A);
>        if (data->P) c_free(data->P);
>        c_free(data);
>    }
>    if (settings) c_free(settings);
>
>    return exitflag;
>};
>```

**csc_matrix构造函数**

>```c
>csc* csc_matrix(c_int m, c_int n, c_int nzmax, c_float *x, c_int *i, c_int *p)
>{
>  csc *M = (csc *)c_malloc(sizeof(csc));
>
>  if (!M) return OSQP_NULL;
>
>  M->m     = m;
>  M->n     = n;
>  M->nz    = -1;
>  M->nzmax = nzmax;
>  M->x     = x;
>  M->i     = i;
>  M->p     = p;
>  return M;
>}
>```

**csc结构体**

Matrix in compressed-column form.

The structure is used internally to store matrices in the triplet form as well,

but the API requires that the matrices are in the CSC format.

>```C
>typedef struct {
>  c_int    nzmax; ///< maximum number of entries
>  c_int    m;     ///< number of rows
>  c_int    n;     ///< number of columns
>  c_int   *p;     ///< column pointers (size n+1); col indices (size nzmax) start from 0 when using triplet format (direct KKT matrix formation)
>  c_int   *i;     ///< row indices, size nzmax starting from 0
>  c_float *x;     ///< numerical values, size nzmax
>  c_int    nz;    ///< number of entries in triplet matrix, -1 for csc
>} csc;
>```



**CSC matrix解释参考资料**

[SciPy.org官网说明]( https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html )

[SciPy课堂]( http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html )

[苹果开发者指南说明]( https://developer.apple.com/documentation/accelerate/creating_sparse_matrices )

[英文博客]( https://www.cnblogs.com/rollenholt/p/5960523.html )
