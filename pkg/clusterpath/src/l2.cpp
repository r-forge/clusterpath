/* -*- compile-command: "make l2.out"  -*- */ 
#include "l2.h"
#include <R_ext/Utils.h>
#include <R.h>

Vector operator-( Vector u, Vector v ) {
  int s = u.size();
  Vector ans( s );
  for( int i = 0 ; i < s ; i++ ) {
    ans[i] = u[i] - v[i];
  }
  return ans;
}

Vector operator+( Vector u, Vector v ) {
  int s = u.size();
  Vector ans( s );
  for( int i = 0 ; i < s ; i++ ) {
    ans[i] = u[i] + v[i];
  }
  return ans;
}

Vector operator*(Vector v,double x){
  unsigned int k;
  for(k=0;k<v.size();k++){
    v[k] *= x;
  }
  return v;
}

Vector operator/(Vector v,double x){
  unsigned int k;
  for(k=0;k<v.size();k++){
    v[k] /= x;
  }
  return v;
}

bool dojoin(JoinPair p){
  if(p.first==p.second)return false;
  return true;
}

JoinPair check_clusters_thresh
(std::list< std::vector<int> > *clusters,
 SymNoDiag const &diffs,
 double thresh
 ){
  std::vector<int> rows,rowsj;
  std::list< std::vector<int> >::iterator it,cj,penultimate;
  int i,j;
  penultimate=clusters->end();
  penultimate--;
  for(it=clusters->begin();it!=penultimate;it++){
    rows = *it;
    i = rows[0];
    for(cj=it,cj++;cj!=clusters->end();cj++){
      rowsj = *cj;
      j = rowsj[0];
      if(diffs.getval(i,j) < thresh)return JoinPair(it,cj);
    }
  }
  return JoinPair(it,it);
}

void take_step
(std::list< std::vector<int> > const &clusters,
 double *alpha,
 double *dir,
 unsigned int N,
 unsigned int P,
 double step
 ){
  std::list< std::vector<int> >::const_iterator it;
  std::vector<int> rows;
  unsigned int i,k;
  for(it=clusters.begin();it!=clusters.end();it++){
    rows = *it;
    i=rows[0];
    for(k=0;k<P;k++){
      alpha[i+k*N] += step * dir[i+k*N];
    }
  }
}

double calc_cost
(std::list< std::vector<int> > const &clusters,
 Matrix const &alpha,
 Matrix const &x,
 SymNoDiag *W,
 SymNoDiag const &diffs,
 double lambda
 ){
  std::list< std::vector<int> >::const_iterator it,cj;
  std::vector<int> rows,rowsj;
  unsigned int i,j,k,l;
  double loss=0,penalty=0;
  for(it=clusters.begin();it!=clusters.end();it++){
    //first calc loss
    rows = *it;
    i=rows[0]; //index to use for alpha/diffs
    Array a=alpha.row(i);
    for(k=0;k<rows.size();k++){//index to use for x/W
      loss += nrm2sq(a-x.row(rows[k]));
      //then calc penalty
      for(cj=it,cj++;cj!=clusters.end();cj++){
	rowsj = *cj;
	j=rowsj[0];//alpha/diffs
	for(l=0;l<rowsj.size();l++){//x/W
	  penalty += W->getval(rows[k],rowsj[l]) * diffs.getval(i,j);
	}
      }
    }
  }
  //Rprintf("loss %f lambda %f penalty %f\n",loss,lambda,penalty);
  //return 0.5 * loss + lambda / ((double)(W->N-1)) * penalty;
  //TDH alternate parameterization
  return 0.5 * loss + lambda * penalty;
}

Results* join_clusters2_restart
(double *x,//array/matrix of data
 SymNoDiag *W,//lower triangle of weight matrix
 unsigned int Px,//problem size
 double lambda,//starting point in regularization path
 double join_thresh, //tolerance for equality of points
 double opt_thresh, //tolerance for optimality
 double lambda_factor,//increase of lambda after optimality
 double smooth,//smoothing parameter
 int maxit,
 int linesearch_freq,//how often to do a linesearch? if 0, never. if
		     //n>0, do n-1 linesearch steps for every
		     //decreasing step size step. set this to 2 if
		     //unsure.
 int linesearch_points,//how many points to check along the gradient
		       //direction. set to 10 if unsure.
 int check_splits,
 int target_cluster,
 int verbose
 ){
  unsigned int N = W->N;
  //W->print();
  double old_lambda=0;
  std::vector<int> rows,rowsj;
  std::vector<int>::iterator rowit,ri,rj;
  std::list< std::vector<int> > clusters,tocheck;
  std::list< std::vector<int> >::iterator it,cj;
  unsigned int i,k,j;
  int tried_restart;
  for(i=0;i<N;i++){
    rows.assign(1,i);
    clusters.push_back(rows);
  }
  double *old_alpha = new double[N*Px];
  double *alpha = new double[N*Px];
  double *xbar = new double[N*Px];
  double *dir = new double[N*Px];
  for(i=0;i<N*Px;i++){
    alpha[i]=xbar[i]=x[i];
  }
  Matrix amat(alpha,N,Px),xmat(x,N,Px);
  SymNoDiag diffs(N);
  diffs.calc_diffs(clusters,amat,nrm2);
  //store initial trivial solution
  Results *results = new Results(N,Px,opt_thresh);
  if(target_cluster==0)results->add(alpha,0,0);
  double weight,diff,step;
  while(clusters.size()>1){
    double grad=opt_thresh;
    int iteration=1;
    tried_restart=0;
    //if we use the general (slower) algorithm for any weights, then
    //split the clusters to individual points
    if(check_splits){
      clusters.clear();
      //reassign original clusters
      for(i=0;i<N;i++){
	rows.assign(1,i);
	clusters.push_back(rows);
      }
      //recopy original xbar
      for(i=0;i<N*Px;i++){
	xbar[i]=x[i];
      }
    }
    while(grad>=opt_thresh){
      R_CheckUserInterrupt();
      //first calc gradients
      grad = 0;
      for(it=clusters.begin();it!=clusters.end();it++){
	rows = *it;
	i = rows[0];
	for(k=0;k<Px;k++){
	  dir[i+k*N] = xbar[i+k*N] - alpha[i+k*N];
	}
	for(cj=clusters.begin();cj!=clusters.end();cj++){
	  if(it!=cj){
	    rowsj = *cj;
	    j=rowsj[0];
	    weight=0;
	    diff = *diffs(i,j);
	    if(diff!=0){
	      if(smooth!=0){
		diff *= diff; //now squared l2 norm
		diff += smooth; //add smoothing parameter under sqrt
		diff = sqrt(diff);//put sqrt back
	      }
	      for(ri=rows.begin();ri!=rows.end();ri++){
		for(rj=rowsj.begin();rj!=rowsj.end();rj++){
		  weight += W->getval(*ri,*rj);
		}
	      }
	      //weight *= lambda / diff / ((double)(N-1)) / ((double)rows.size());
	      weight *= lambda / diff / ((double)rows.size());
	      for(k=0;k<Px;k++){
		dir[i+k*N] += weight * (alpha[j+k*N]-alpha[i+k*N]);
	      }
	    }
	  }
	}
	grad += nrm2(Array(dir+i,N,Px));
      }
      //store this iteration
      //results->add(alpha,lambda,grad);
      //then take a step
      if(linesearch_freq==0 || (iteration % linesearch_freq)==0 ){
	//Decreasing step size
	//TDH and pierre 18 jan 2011 try sqrt dec step size
	step=1/((double)iteration);
	//step=1/sqrt((double)iteration);
	if(verbose>=2)Rprintf("grad %f step %f it %d\n",grad,step,iteration);
	take_step(clusters,alpha,dir,N,Px,step);
      }else{
	double cost_here,cost_step;
	std::map<double,double> cost_steps;
	std::map<double,double>::iterator step1,step2;
	for(i=0;i<N*Px;i++)old_alpha[i]=alpha[i];//copy alpha
	//compare current cost to cost after stepping in gradient direction
	cost_here=cost_step=calc_cost(clusters,amat,xmat,W,diffs,lambda);
	step = 0;
	cost_steps.insert(std::pair<double,double>(cost_here,0));
	while(cost_step<=cost_here){
	  take_step(clusters,alpha,dir,N,Px,1);
	  step += 1;
	  diffs.calc_diffs(clusters,amat,nrm2);
	  cost_step=calc_cost(clusters,amat,xmat,W,diffs,lambda);
	  if(verbose>=2)
	Rprintf("cost %.10f step %f cost_here %f\n",cost_step,step,cost_here);
	R_CheckUserInterrupt();
	  cost_steps.insert(std::pair<double,double>(cost_step,step));
	}
	for(int cuts=0;cuts<linesearch_points;cuts++){
	  step1=step2=cost_steps.begin();
	  step2++;
	  step = (step1->second + step2->second)/2;
	  for(i=0;i<N*Px;i++){
	    alpha[i]=old_alpha[i];
	  }
	  take_step(clusters,alpha,dir,N,Px,step);
	  diffs.calc_diffs(clusters,amat,nrm2);
	  cost_step=calc_cost(clusters,amat,xmat,W,diffs,lambda);
	  if(verbose>=2)Rprintf("cost %.10f step %f %d\n",cost_step,step,cuts);
	  cost_steps.insert(std::pair<double,double>(cost_step,step));
	}
	cost_steps.clear();
      }
      if(iteration++ > maxit){
	if(tried_restart){
	  Rprintf("max iteration %d exit\n",maxit);
	  delete old_alpha;
	  delete alpha;
	  delete xbar;
	  delete dir;
	  return results;
	}else{
	  if(verbose>=1)Rprintf("max iterations, trying restart from x\n");
	  tried_restart=1;
	  iteration=1;
	  for(i=0;i<N*Px;i++)alpha[i]=x[i];
	}
      }
      //calculate differences
      diffs.calc_diffs(clusters,amat,nrm2);
      //check for joins
      JoinPair tojoin;
      while(dojoin(tojoin=check_clusters_thresh(&clusters,diffs,join_thresh))){
	//if(verbose>=1)
	//  Rprintf("join: %d %d\n",tojoin.first->front(),tojoin.second->front());
	int ni=tojoin.first->size();
	int nj=tojoin.second->size();
	i=tojoin.first->front();
	j=tojoin.second->front();
	tojoin.first->insert(tojoin.first->end(),
			    tojoin.second->begin(),
			    tojoin.second->end());
	for(k=0;k<Px;k++){
	  alpha[i+k*N] = (alpha[i+k*N]*ni + alpha[j+k*N]*nj)/(ni+nj);
	  xbar[i+k*N] = (xbar[i+k*N]*ni + xbar[j+k*N]*nj)/(ni+nj);
	}
	clusters.erase(tojoin.second);
	iteration=1;
	if(clusters.size()>1){
	  diffs.calc_diffs(clusters,amat,nrm2);//inefficient
	}else{
	  grad=0;//so we can escape from the last optimization loop
	}
      }
    }//while(grad>=opt_thresh)
    if(verbose>=1)
    Rprintf("solution iteration %d lambda %f nclusters %d\n",
	   iteration,lambda,(int)clusters.size());
    
    if(target_cluster == 0){
      //for each cluster, there may be several points. we store the
      //alpha value just in the row of the first point. thus here we
      //copy this value to the other rows before copying the optimal
      //alpha to results.
      for(it=clusters.begin();it!=clusters.end();it++){
	rows = *it;
	if(rows.size()>1){
	  for(i=1;i<rows.size();i++){
	    for(k=0;k<Px;k++){
	      alpha[rows[i]+k*N] = alpha[rows[0]+k*N];
	    }
	  }
	}
      }
      results->add(alpha,lambda,grad);
    }
    //haven't yet reached the target number of clusters, multiply
    //lambda by lambda_factor and continue along the path
    if((int)clusters.size()>target_cluster){
      old_lambda=lambda;
      lambda *= lambda_factor;
    }
    //if we have passed the target cluster number then decrease
    //lambda and go look for it!
    if((int)clusters.size()<target_cluster){
      if(verbose>=1){
	Rprintf("missed target %d, going back for it\n",target_cluster);
      }
      lambda = (lambda+old_lambda)/2;
      clusters.clear();
      //reassign original clusters
      for(i=0;i<N;i++){
	rows.assign(1,i);
	clusters.push_back(rows);
      }
      //recopy original xbar
      for(i=0;i<N*Px;i++){
	xbar[i]=x[i];
      }
    }
    //this is the number of clusters that we were looking for,
    //save and quit!
    if((int)clusters.size()==target_cluster){
      for(it=clusters.begin();it!=clusters.end();it++){
	rows = *it;
	if(rows.size()>1){
	  for(i=1;i<rows.size();i++){
	    for(k=0;k<Px;k++){
	      alpha[rows[i]+k*N] = alpha[rows[0]+k*N];
	    }
	  }
	}
      }
      results->add(alpha,lambda,grad);
      if(verbose>=1)Rprintf("got target cluster %d exit\n",target_cluster);
      delete old_alpha;
      delete alpha;
      delete xbar;
      delete dir;
      return results;
    }
  }	
  //TODO: consolidate cleanup... just use data structures that
  //automatically clean themselves up when the function exits.
  delete old_alpha;
  delete alpha;
  delete xbar;
  delete dir;
  return results;
}

