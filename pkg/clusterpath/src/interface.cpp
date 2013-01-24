#include <R.h>
#include <Rinternals.h>
#include <Rcpp.h>
#include "l1.h"
#include "l2.h"

using namespace Rcpp ;

extern "C" {
 
SEXP tree2list(Cluster *c){
  int nrow=c->total;
  //printf("%d total\n",nrow);
  SEXP alpha,lambda,i,result,names;
  PROTECT(alpha = allocVector(REALSXP, nrow));
  PROTECT(lambda = allocVector(REALSXP, nrow));
  PROTECT(i = allocVector(INTSXP, nrow));
  int row=0;
  add_results(c,REAL(alpha),REAL(lambda),INTEGER(i),&row);
  PROTECT(names = allocVector(STRSXP, 3));
  SET_STRING_ELT(names,0,mkChar("alpha"));
  SET_STRING_ELT(names,1,mkChar("lambda"));
  SET_STRING_ELT(names,2,mkChar("i"));
  PROTECT(result = allocVector(VECSXP, 3));
  namesgets(result, names);
  SET_VECTOR_ELT(result, 0, alpha);
  SET_VECTOR_ELT(result, 1, lambda);
  SET_VECTOR_ELT(result, 2, i);
  
  UNPROTECT(5);
  return(result);
}

// we process each dimension individually using this function
SEXP join_clusters_convert(SEXP xR){
  //NumericVector x(xR);
  Cluster *tree = make_clusters_l1(REAL(xR),length(xR));
  SEXP L=tree2list(tree);
  delete_tree(tree);
  return L;
}

//construct a list that will be used in R to return a data.frame
SEXP res2list(Results *r){
  unsigned int nrow=r->iterations.size() * r->n,i,k,row=0;
  SEXP alpha, ivec, lambda, grad, result, names;
  PROTECT(alpha = allocMatrix(REALSXP,nrow,r->p));
  PROTECT(ivec = allocVector(INTSXP, nrow));
  PROTECT(lambda = allocVector(REALSXP,nrow));
  PROTECT(grad = allocVector(REALSXP, nrow));
  std::list<Iteration*>::iterator l;
  //r->print();
  for(l=r->iterations.begin();l!=r->iterations.end();l++){
    for(i=0;i<r->n;i++){
      for(k=0;k<r->p;k++){
	//printf("%d %d %d %d %d\n",nrow,k,row,i,r->n);
	REAL(alpha)[nrow*k+row] = (*l)->alpha[i+k*r->n];
      }
      REAL(lambda)[row] = (*l)->lambda;
      REAL(grad)[row] = (*l)->grad;
      INTEGER(ivec)[row] = i;
      row++;
    }
  }
  PROTECT(names = allocVector(STRSXP, 4));
  SET_STRING_ELT(names,0,mkChar("alpha"));
  SET_STRING_ELT(names,1,mkChar("i"));
  SET_STRING_ELT(names,2,mkChar("lambda"));
  SET_STRING_ELT(names,3,mkChar("grad"));
  PROTECT(result = allocVector(VECSXP, 4));
  namesgets(result, names);
  SET_VECTOR_ELT(result,0,alpha);
  SET_VECTOR_ELT(result,1,ivec);
  SET_VECTOR_ELT(result,2,lambda);
  SET_VECTOR_ELT(result,3,grad);
  
  UNPROTECT(6);
  return result;
}

RcppExport SEXP calcW_convert(SEXP x_R,SEXP gamma){
  NumericMatrix x(x_R);
  unsigned int i, n = x.nrow();
  NumericVector w_R(n*(n-1)/2);
  SymNoDiag W(n);
  W.calcW(&x[0],x.ncol(),NumericVector(gamma)[0]);
  for(i=0;i<w_R.length();i++){
    w_R[i]=W.data[i];
  }
  return w_R;
}

RcppExport SEXP join_clusters2_restart_convert
(SEXP x_R,SEXP w_R,
 SEXP lambda,SEXP join_thresh,SEXP opt_thresh,SEXP lambda_factor,SEXP smooth,
 SEXP maxit,SEXP linesearch_freq,SEXP linesearch_points,SEXP check_splits,
 SEXP target_cluster,
 SEXP verbose
 ){
  NumericMatrix x(x_R);
  unsigned int i, N=x.nrow();
  SymNoDiag W(N);
  NumericVector wvec(w_R);
  for(i=0;i<W.length;i++){
    W.data[i] = wvec[i];
  }
  Results *r = join_clusters2_restart
    (&x[0],
     &W,
     x.ncol(),
     NumericVector(lambda)[0],
     NumericVector(join_thresh)[0],
     NumericVector(opt_thresh)[0],
     NumericVector(lambda_factor)[0],
     NumericVector(smooth)[0],
     IntegerVector(maxit)[0],
     IntegerVector(linesearch_freq)[0],
     IntegerVector(linesearch_points)[0],
     IntegerVector(check_splits)[0],
     IntegerVector(target_cluster)[0],
     IntegerVector(verbose)[0]);
  SEXP L = res2list(r);
  delete r;
  return L;
}

}
