#include <R.h>
#include <Rinternals.h>
#include "l1.h"
#include "l2.h"

// #include <Rcpp.h>
// using namespace Rcpp ;

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
  Rprintf("%f %f %f %d\n",REAL(alpha)[0],REAL(lambda)[0],REAL(grad)[0],
	  INTEGER(ivec)[0]);
  UNPROTECT(6);
  return result;
}
 
SEXP join_clusters2_restart_convert
(SEXP x_R,SEXP w_R,
 SEXP lambda,SEXP join_thresh,SEXP opt_thresh,SEXP lambda_factor,SEXP smooth,
 SEXP maxit,SEXP linesearch_freq,SEXP linesearch_points,SEXP check_splits,
 SEXP target_cluster,
 SEXP verbose
 ){
  SEXP x;
  PROTECT(x = allocMatrix(REALSXP, nrows(x_R), ncols(x_R)));
  //NumericMatrix x(x_R);
  unsigned int i, N=nrows(x_R);
  SymNoDiag W(N);
  //NumericVector wvec(w_R);
  for(i=0; i<W.length; i++){
    W.data[i] = REAL(w_R)[i];
  }
  Results *r = join_clusters2_restart
    (REAL(x),
     &W,
     ncols(x),
     REAL(lambda)[0],
     REAL(join_thresh)[0],
     REAL(opt_thresh)[0],
     REAL(lambda_factor)[0],
     REAL(smooth)[0],
     INTEGER(maxit)[0],
     INTEGER(linesearch_freq)[0],
     INTEGER(linesearch_points)[0],
     INTEGER(check_splits)[0],
     INTEGER(target_cluster)[0],
     INTEGER(verbose)[0]);
  SEXP L = res2list(r);
  delete r;

  UNPROTECT(1);
  return L;
}

}

