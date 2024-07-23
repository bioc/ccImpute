#define EIGEN_USE_BLAS
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::plugins(openmp)]]

#include <algorithm>
#include <random>
#include <vector>
#include <unordered_map>

using Eigen::Map;
using Eigen::Ref;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;

// Index sort of vector x with s elements
//
// @param x Input vector.
// @param s Number of elements in x.
// @return Indexes of elements corresponding to the sorted order of elements.
VectorXi index_sort(const Ref<VectorXd> x, unsigned int s) {
    VectorXi idx(s);
    std::iota(idx.data(), idx.data()+s, 0);
    std::sort(idx.data(), idx.data()+s, [&](int i,int j){return x[i] < x[j];});
    return idx;
}

// Compute ranking of a single column vector.
// Helper function needed to compute Spearman correlation.
//
// @param x Input vector to be ranked.
// @param ranked Reference to a vector where the rankings will be stored.
void rank_fast(Ref<VectorXd> x, Ref<VectorXd> ranked) {
    unsigned int sz = x.size();
    VectorXi is = index_sort(x, sz);
    
    for (unsigned int n, i = 0; i < sz; i += n) {
        n = 1;
        while (i + n < sz && x[is[i]] == x[is[i + n]]) ++n;
        double tied_val = i + (n + 1) / 2.;
        for (unsigned int k = 0; k < n; k++) {
            ranked[is[i + k]] = tied_val;
        }
    }
}

//' Computes rankings for each column of a matrix in parallel.
//'
//' @param x The input matrix to be ranked.
//' @param n_cores The number of CPU cores to use for parallel processing.
//' @return  A matrix where each column contains the rankings for the
//' corresponding column in the input matrix.
// [[Rcpp::export]]
Eigen::MatrixXd colRanks_fast(Eigen::Map<Eigen::MatrixXd> x, 
                  const unsigned int n_cores) {
    Eigen::MatrixXd rank_m(x.rows(), x.cols());
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores) shared(rank_m, x)
    #endif
    for(unsigned int i = 0; i < x.cols(); ++i){
        rank_fast(x.col(i), rank_m.col(i));
    }
    return rank_m;
}

//' Computes rankings for each column of a matrix in parallel.
//'
//' @param x The input matrix to be ranked.
//' @param n_cores The number of CPU cores to use for parallel processing.
//' @return  A matrix where each column contains the rankings for the
//' corresponding column in the input matrix.
// [[Rcpp::export]]
Eigen::MatrixXd sparseColRanks_fast(const Eigen::MappedSparseMatrix<double> x,
                                    const unsigned int n_cores) {
    Eigen::MatrixXd rank_m(x.rows(), x.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) shared(rank_m, x)
#endif
    for(unsigned int i = 0; i < x.cols(); ++i){
        Eigen::VectorXd temp = x.col(i);
        rank_fast(temp, rank_m.col(i));
    }
    return rank_m;
}

//' Computes a Weighted Pearson
//'
//' This function calculates weighted Pearson correlation matrix
//'
//' @param x The input matrix (dense), where each column represents a set of
//' observations.
//' @param w A vector of weights, one for each observation (must have the same
//' number of elements as rows in `x`).
//' @param n_cores The number of CPU cores to utilize for parallel computation.
//' @return A weighted Pearson correlation matrix 
// [[Rcpp::export]]
Eigen::MatrixXd wCor_fast(Eigen::Map<Eigen::MatrixXd> x,
                         const Eigen::Map<Eigen::VectorXd> w,
                         const unsigned int n_cores){
    Eigen::setNbThreads(n_cores);
    
    Eigen::VectorXd w_s = w.array()/w.array().sum();
    VectorXd center = (w_s.adjoint() * x);
    x.rowwise() -= center.adjoint();
    x.array().colwise() *= w_s.array().sqrt();
    x.array().rowwise() /= x.colwise().norm().array();
    return x.adjoint() * x;
}

//' Computes a Pearson
//'
//' This function calculates a Pearson correlation matrix
//' 
//' @param x The input matrix, where each column represents a set of
//' observations.
//' @param n_cores The number of CPU cores to utilize for parallel computation
//' (optional, defaults to 1).
//' @return A Pearson correlation matrix if `useRanks` is `false`. If
//' `useRanks` is `true`, returns a Spearman correlation matrix.
// [[Rcpp::export]]
Eigen::MatrixXd cor_fast(Eigen::Map<Eigen::MatrixXd> x,
                         const unsigned int n_cores){
    Eigen::setNbThreads(n_cores);
    x.rowwise() -=  x.colwise().mean();
    x.array().rowwise() /= x.colwise().norm().array();
    return x.adjoint() * x;
}

//' Computes Row Variances Efficiently
//' @param x A numeric dense matrix for which to compute row variances.
//' @param n_cores The number of cores to utilize for parallel processing.
//' @return A numeric vector containing the variance for each row of the input
//' matrix.
//' 
//' @examples
//'
//' library(Matrix)
//' rand_vals <- sample(0:10,1e4,replace=TRUE, p=c(0.99,rep(0.001,10)))
//' x <- as.matrix(Matrix(rand_vals,ncol=5))
//' cores <- 2
//' vars_vector <- rowVars_fast(x, cores)
//' 
//' @export
// [[Rcpp::export]]
Eigen::VectorXd rowVars_fast(const Eigen::Map<Eigen::MatrixXd> x, int n_cores){
    unsigned int rc = x.rows(), cc = x.cols();
    Eigen::VectorXd vars(rc);
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores)
    #endif
    for(unsigned int row=0; row < rc; ++row){
        const Eigen::VectorXd & r = x.row(row);
        double sum1 = r.squaredNorm();
        double sum2 = r.sum();
        vars[row] = (sum1-sum2*sum2/cc)/(cc-1);
    }
    return vars;
}

//' This function calculates an average consensus matrix from a set of
//' clustering solutions. It filters out values below a specified minimum
//' threshold (`consMin`) and normalizes the remaining non-zero columns to sum
//' to 1.
//'
//' @param dat An integer matrix where each column represents a different
//' clustering solution (cluster assignments for each data point).
//' @param consMin The minimum consensus value to retain. Values below this
//' threshold are set to zero.
//' @param n_cores The number of cores to use for parallel processing.
//' This can speed up the normalization step.
//' @return A processed consensus matrix where each element (i, j) represents
//' the proportion of times data points i and j were assigned to the same
//' cluster with filtering and normalization applied.
// [[Rcpp::export]]
Eigen::MatrixXd getConsMtx(const Eigen::Map<Eigen::MatrixXi> dat, 
                           double consMin, int n_cores) {
    unsigned int n = dat.rows();
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n, n);
    for (unsigned int j = 0; j < dat.cols(); ++j) {
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int k = i + 1; k < n; ++k) {
                if (dat(i, j) == dat(k, j)) {
                    ++res(i, k);
                    ++res(k, i);
                }
            }
        }
    }
    res /= dat.cols();
    
    if(consMin < 1)
        res = (res.array() < consMin).select(0, res);
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores)
    #endif
    for (int i = 0; i < res.cols(); i++){
        double sum = res.col(i).sum();
        if(sum > 0)
            res.col(i)/=res.col(i).sum();
    }
    return res;
}
 
//' Computes Means and Standard Deviations for Scaling
//' 
//' @param x A numeric matrix representing the gene expression data, where rows
//' are genes and columns are samples.
//' @param n_cores The number of cores to use for parallel processing.
//' @return A list containing:
//'   * `means`: A numeric vector of column means.
//'   * `sds`: A numeric vector of column standard deviations.
// [[Rcpp::export]]
Rcpp::List getScale(const Eigen::Map<Eigen::MatrixXd> x, int n_cores){
    unsigned int cols = x.cols();
    unsigned int rows = x.rows();
    VectorXd means(cols);
    VectorXd sds(cols);
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores)
    #endif
    for(unsigned int i=0; i < cols; ++i){
        means[i] = x.col(i).mean();
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores)
    #endif
    for(unsigned int i=0; i < cols; ++i){
        VectorXd diff = x.col(i).array() - means[i];
        sds[i] = sqrt(diff.dot(diff)/rows);
    }
    return Rcpp::List::create(Rcpp::Named("means") = means, 
                              Rcpp::Named("sds") = sds); 
}

//' Computes imputed expression matrix using linear eq solver
//' @param cm processed consensus matrix
//' @param em expression matrix
//' @param ids location of values determined to be dropout events
//' @param n_cores number of cores to use for parallel computation.
//' @return imputed expression matrix
// [[Rcpp::export]]
Eigen::MatrixXd solver(const Eigen::Map<Eigen::MatrixXd> cm,
                           Eigen::Map<Eigen::MatrixXd> em,
                           const Eigen::Map<Eigen::MatrixXi> ids,
                           const int n_cores) {
    //map row to columns
    std::unordered_map<unsigned int, std::vector<unsigned int>> row2cols;
    std::vector<unsigned int> keys;
    
    // Map out the rows for each column
    for(unsigned int i=0; i<ids.rows(); ++i){
        // account for the fact that r indices are 1 based
        unsigned int row = ids(i,0) - 1;
        if (row2cols.find(row) == row2cols.end()) {
            keys.push_back(row);
            row2cols[row] = std::vector<unsigned int>();
        }
        // account for the fact that r indices are 1 based
        row2cols[row].push_back(ids(i,1)-1);
    }
    Eigen::MatrixXd em_t = em.transpose();

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores)
    #endif
    for(unsigned int i=0; i < keys.size(); ++i){
        unsigned int row_i = keys[i];
        std::vector<unsigned int> cidx = row2cols[row_i];
        unsigned int col_count = cidx.size();
        if(col_count == 1){
            em_t(cidx[0],row_i) = cm.col(cidx[0]).dot(em_t.col(row_i));
        }
        else{
            Eigen::MatrixXd A(col_count, col_count);

            for(unsigned int j=0; j < col_count; ++j){
                for(unsigned int k=0; k < col_count; ++k){
                    A(j,k) = j==k ? 1:-cm(cidx[j], cidx[k]);
                }
            }
            Eigen::VectorXd b(col_count);
            for(unsigned int j=0; j < col_count; ++j){
                b[j] = cm.col(cidx[j]).dot(em_t.col(row_i));
            }
            Eigen::VectorXd solution = A.llt().solve(b);
            for(unsigned int j = 0; j < col_count; ++j){
                em_t(cidx[j],row_i) = solution(j);
            }
        }
    }
    return em_t.transpose();
}

//' Fast Calculation of "Dropout" values
//' @param cm A numeric matrix representing the consensus matrix.
//' @param em A dense numeric matrix representing the gene expression data,
//' where rows are genes and columns are samples.
//' @param ids An integer matrix specifying the row and column indices of
//' entries for which to calculate importance scores. Each row of `ids` should
//' contain two integers: the row index (gene) and column index (sample) in the
//' `em` matrix.
//' @param n_cores The number of cores to use for parallel processing.
//' @return A numeric vector of imputed dropout values, corresponding to the
//' entries specified in the `ids` matrix.
// [[Rcpp::export]]
Eigen::VectorXd solver2(const Eigen::Map<Eigen::MatrixXd> cm,
                            const Eigen::Map<Eigen::MatrixXd> em,
                            const Eigen::Map<Eigen::MatrixXi> ids,
                            const int n_cores) {
    Eigen::MatrixXd em_t = em.transpose();
    Eigen::VectorXd imp(ids.rows());
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores)
    #endif
    for(unsigned int i=0; i < ids.rows(); ++i){
        int row = ids(i,0)-1, col = ids(i,1)-1;
        Eigen::VectorXd em_ =  (em_t.col(row).array() > 0).cast<double>();
        double div2 = em_.dot(cm.col(col));
        imp(i) = em_t.col(row).dot(cm.col(col))/div2;
    }
    return imp;
}

//' Fast Calculation of "Dropout" values
//' @param cm A numeric matrix representing the consensus matrix.
//' @param em A sparse numeric matrix representing the gene expression data,
//' where rows are genes and columns are samples.
//' @param ids An integer matrix specifying the row and column indices of
//' entries for which to calculate importance scores. Each row of `ids` should
//' contain two integers: the row index (gene) and column index (sample) in the
//' `em` matrix.
//' @param n_cores The number of cores to use for parallel processing.
//' @return A numeric vector of imputed dropout values, corresponding to the
//' entries specified in the `ids` matrix.
// [[Rcpp::export]]
Eigen::VectorXd sparseSolver2(const Eigen::Map<Eigen::MatrixXd> cm,
                                    const Eigen::MappedSparseMatrix<double> em,
                                    const Eigen::Map<Eigen::MatrixXi> ids, 
                                    const int n_cores){
    Eigen:: SparseMatrix<double> em_t = em.transpose();

    Eigen::VectorXd imp(ids.rows());

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores)
    #endif
    for(unsigned int i=0; i < ids.rows(); ++i){
        int row = ids(i,0)-1, col = ids(i,1)-1;
        double n = 0, d = 0;
        Eigen::SparseVector<double> temp = em_t.col(row);
        
        for (Eigen::SparseVector<double>::InnerIterator it(temp); it; ++it) {
            int i = it.index();
            double cm_val = cm.col(col)(i);
            n += cm_val*it.value();
            d += cm_val;
        }
        imp(i) = n/d;
    }
    return imp;
}
