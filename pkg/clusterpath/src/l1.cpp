#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "l1.h"

#define NEW_EVENT_THRESH 1e-10
#define SIGN(x)( ((x)==0)?(0):( (((x)>0)?(1):(-1)) ) )

//calculate lambda where cluster c1 would hit cluster c2
double lnew(Cluster *c1,Cluster *c2){
  return((c2->alpha-c2->v*c2->lambda+c1->lambda*c1->v-c1->alpha)/(c1->v-c2->v));
}

// used for adding/deleting previous and next joining events
void insert_new_event(Events *events_ptr,Cluster *newc,
		      Cluster *neighbor_cluster,Clusters::iterator top_it,
		      Cluster *old_cluster){
  if(old_cluster->merge_with_next != (*events_ptr).end()){
    (*events_ptr).erase(old_cluster->merge_with_next);//constant
  }
  Cluster *top_cluster = *top_it;
  Events::iterator new_event;
  Event e=Event(lnew(neighbor_cluster,newc),top_it);
  //printf("event:%f cluster:%f\n",e.first,newc->lambda);
  if(e.first + NEW_EVENT_THRESH >= newc->lambda){
    //printf("inserting!\n");
    new_event=(*events_ptr).insert(e);//logarithmic
  }else{
    new_event=(*events_ptr).end();
  }
  top_cluster->merge_with_next = new_event;
}

// assume clusters are sorted already, with v already calculated, and
// that no splits will occur.
Cluster* join_clusters(Clusters clusters){
  int ni,nj;
  Clusters::iterator it,prev_it,next_it;
  Events events;
  // set to 0 here to avoid compiler warnings.
  Cluster *c=0,*ci,*cj,*prev_cluster=0;
  Event e;
  Events::iterator e_it,new_event;
  double lambda;
  for(prev_it=clusters.begin(),it=clusters.begin(),it++;
      it!=clusters.end();
      it++,prev_it++){
    c=*it;
    prev_cluster=*prev_it;
    // calculate initial events list from intersection of all adjacent
    // lines
    lambda=(prev_cluster->alpha - c->alpha)/(c->v - prev_cluster->v);
    if(lambda>0){//only insert if greater than 0!
      e=Event(lambda,prev_it);
      e_it = events.insert(e);//logarithmic
      prev_cluster->merge_with_next = e_it;
    }
    prev_cluster=c;
  }

  while(events.size()>0){
    do{
      e_it=events.begin();//constant
      lambda=e_it->first;
      // follow links and define current clusters
      it=next_it=e_it->second;
      ci= *it;
      ++next_it;
      cj= *next_it;
      ni=ci->i.size();
      nj=cj->i.size();

      //print_events(events);
      //print_clusters(clusters);

      //Make the new cluster object
      c = new Cluster;
      c->i.reserve(ni+nj);
      c->i.insert(c->i.end(),ci->i.begin(),ci->i.end());
      c->i.insert(c->i.end(),cj->i.begin(),cj->i.end());
      c->total=ci->total+cj->total+c->i.size();
      //calculate new velocity of joined lines
      c->v=(ni*ci->v + nj*cj->v)/(ni+nj);
      //store other data...
      c->lambda=lambda;
      c->alpha= ci->alpha + (c->lambda - ci->lambda)*ci->v;
      c->child1=ci;
      c->child2=cj;

      //set up pointers to neighboring clusters
      prev_it=e_it->second;
      prev_it--;
      ++next_it;

      //We have prev-ci-cj-next, we remove ci==it from the list (erase)
      //and set cj==it after to a pointer to new cluster
      it=clusters.erase(it);//linear on the number of elements erased==constant
      *it=c;

      if(it!=clusters.begin()){//constant
	prev_cluster= *prev_it;
	insert_new_event(&events,c,prev_cluster,prev_it,prev_cluster);
      }
      if(next_it != clusters.end()){//constant
	insert_new_event(&events,c,*next_it,it,cj);
      } 
      events.erase(e_it);//constant
    }while( events.begin()->first == lambda );
  }
  return c;
}

bool compare_alpha(Cluster *lhs, Cluster *rhs){
  if(lhs->alpha<rhs->alpha)return true;
  else return false;
}

Cluster* make_clusters_l1(double *x, int N){
  double v,sign;
  Cluster *c,*next;
  Clusters::iterator it,next_it,other;
  Clusters clusters;
  //create new list of clusters for this dimension
  for(int i=0;i<N;i++){
    c = new Cluster;
    c->alpha=x[i];
    (c->i).push_back(i);
    c->lambda=0.0;
    c->child1=NullCluster;
    c->child2=NullCluster;
    clusters.push_back(c);//constant
  }

  //print_clusters(clusters);
  clusters.sort(compare_alpha);
  //print_clusters(clusters);

  next_it=it=clusters.begin();
  next_it++;
  while(next_it!=clusters.end()){
    c=*it;
    next=*next_it;
    // check if it and next_it have the same alpha value
    if(next->alpha == c->alpha){
      //then merge the 2
      c->i.insert(c->i.end(),next->i.begin(),next->i.end());
      next_it=clusters.erase(next_it);
    }else{
      next_it++;
      it++;
    }
  }
  
  //print_clusters(clusters);

  // calculate cluster velocities for identity weights. first element
  // is smallest alpha value, with N-1 other clusters above it. thus
  // it has a velocity of N-1. Next cluster up has N-2 clusters above
  // and 1 cluster below = velocity of N-3, etc. UNLESS there are
  // multiple points with the same exact value, in which case we need
  // to scale the velocity by cluster size.
  for(it=clusters.begin();it!=clusters.end();it++){
    (*it)->total = (*it)->i.size();
    v=0.0;
    sign=-1.0;
    for(other=clusters.begin();other!=clusters.end();other++){
      if(other==it){
	sign=1.0;
      }else{
	v += sign * (*other)->i.size();
      }
    }
    //TDH alternate parameterization
    (*it)->v = v;
    //(*it)->v = v/((double)N-1);
  }
  //print_clusters(clusters);
  return join_clusters(clusters);
}
// This will allocate memory for new clusters, then it will be the job
// of the caller to delete them later!
Clusters clustermat_l1(double *x, int N, int P){
  int k;
  Clusters dims;
  // dims is an array of cluster pointers, which will be the final
  // cluster in each of the dimensions in the clustering.
  Clusters clusters;
  for(k=0;k<P;k++){
    dims.push_back(make_clusters_l1(x+N*k,N));
  }
  return dims;
}

void delete_tree(Cluster *c){
  Cluster *c1=c->child1,*c2=c->child2;
  delete c;
  if(c1 != NullCluster)delete_tree(c1);
  if(c2 != NullCluster)delete_tree(c2);
}

void add_results(Cluster *c,double *alpha,
		 double *lambda,int *i,int *row){
  unsigned int k;
  for(k=0;k<c->i.size();k++){
    i[*row] = c->i[k];
    alpha[*row]=c->alpha;
    lambda[*row]=c->lambda;
    (*row)++;
  }
  if(c->child1 != NullCluster){
    add_results(c->child1,alpha,lambda,i,row);
    add_results(c->child2,alpha,lambda,i,row);
  }
}

