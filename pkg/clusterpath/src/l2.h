/* -*- compile-command: "make l2.out"; compilation-read-command: nil; gdb-many-windows: 1 -*- */ 
#include <map>
#include <vector>
#include <cstdio>
#include <cmath>
#include <list>

#define MIN(I,J) (((I)<(J))?(I):(J))
#define MAX(I,J) (((I)>(J))?(I):(J))

typedef std::vector<double> Vector;
typedef std::pair< std::list< std::vector<int> >::iterator,
  std::list< std::vector<int> >::iterator > JoinPair;

Vector operator*(Vector,Vector);
Vector operator*(Vector,double);
Vector operator+(Vector,Vector);

template <class T> double nrm2sq(T x){
  unsigned int k;
  double norm=0.0;
  for(k=0;k<x.size();k++){
    norm += x[k] * x[k];
  }
  return norm;
}

template <class T> double nrm2(T x){
  return sqrt(nrm2sq(x));
}
  
class Array {
 public:
  double *data;
  unsigned int length;
  int offset;
  Array(double *x,int elem_offset,int len){
    this->data = x;
    this->length = len;
    this->offset = elem_offset;
  }
  unsigned int size(){
    return this->length;
  }
  double operator[](int i)const{
    return this->data[i*this->offset];
  }
};

//thin wrapper around a pointer
class Matrix {
 public:
  double *data;
  unsigned int n,p;
  Matrix(double *x,unsigned int N,unsigned int P){
    this->data = x;
    this->n = N;
    this->p = P;
  }
  Array row(unsigned int i)const{
    return Array(this->data+i,this->n,this->p);
  }
  Array col(unsigned int j)const{
    return Array(this->data+j*this->n,1,this->n);
  }
};

inline
Vector operator-(Array const &u,Array const &v ){
  int s=u.length;
  Vector ans(s);
  for(int k=0;k<s;k++){
    ans[k] = u[k] - v[k];
  }
  return ans;
}

class SymNoDiag {
 public:
  double *data;
  unsigned int N;
  unsigned int length;
  int delete_ptr;
  SymNoDiag(int size, double *ptr){
      this->N = size;
      this->length = N*(N-1)/2;
      this->data = ptr;
      this->delete_ptr = 0;
  }
  SymNoDiag(int size){
    this->N = size;
    this->length = N*(N-1)/2;
    this->data = new double[this->length];
    this->delete_ptr = 1;
  }
  ~SymNoDiag(){
      if(this->delete_ptr){
	  delete this->data;
      }
  }
  int index(int first,int second)const{
    int I,J;
    if(first<second){
      I=first;
      J=second;
    }else{
      J=first;
      I=second;
    }
    return I*this->N - (I+1)*I/2 +J-I-1;
  }
  double getval(int first,int second)const{
    return this->data[this->index(first,second)];
  }
  double* getptr(int first,int second){
    return this->data + this->index(first,second);
  }
  double* operator()(int first,int second){
    return this->getptr(first,second);
  }
  void calcW(double *x,unsigned int P,double gamma){
    std::list< std::vector<int> > clusters;
    std::vector<int> rows;
    unsigned int i;
    for(i=0;i<this->N;i++){
      rows.assign(1,i);
      clusters.push_back(rows);
    }
    this->calc_diffs(clusters,Matrix(x,this->N,P),nrm2sq);
    for(i=0;i<this->length;i++){
      this->data[i]=exp(-gamma * this->data[i]);
    }
  }
  void calc_diffs
    (std::list< std::vector<int> > const &clusters,
     Matrix const &x,
     double(*FUN)(Vector)
     ){
    std::list< std::vector<int> >::const_iterator it,cj,penultimate;
    std::vector<int> rows,rowsj;
    int i,j;
    penultimate = clusters.end();
    penultimate--;
    for(it=clusters.begin();it!=penultimate;it++){
      rows = *it;
      i=rows[0];
      for(cj=it,cj++;cj!=clusters.end();cj++){
	rowsj = *cj;
	j=rowsj[0];
	*this->getptr(i,j) = FUN(x.row(i)-x.row(j));
      }
    }
  }
  void print()const{
    unsigned int i,j;
    for(i=0;i<this->N-1;i++){
      for(j=i+1;j<this->N;j++){
	printf("%5.2f ",this->getval(i,j));
      }
      printf("\n");
    }
  }
};


class Iteration {
 public:
  double lambda;
  double grad; /* l2 norm of the gradient at this iteration. compare
		  to thresh to see if this is solved. */
  std::vector<double> alpha;//matrix at this iteration
  Iteration(double l,double g,unsigned int s){
    this->lambda = l;
    this->grad = g;
    this->alpha.reserve(s);
  }
};

struct Results {
  unsigned int n,p;
  double thresh; // tolence for solved
  std::list<Iteration*> iterations;
  Results(unsigned int N,unsigned int P,double opt_thresh){
    this->n = N;
    this->p = P;
    this->thresh = opt_thresh;
  }
  ~Results(){
    std::list<Iteration*>::iterator it;
    for(it=this->iterations.begin();it!=this->iterations.end();it++){
      delete *it;
    }
  }
  Iteration* newit(double l,double g){
    return new Iteration(l,g,this->n*this->p);
  }
  void add(double *alpha,double l,double grad){
    Iteration *r = this->newit(l,grad);
    for(unsigned int i=0;i< this->n * this->p;i++){
      r->alpha[i]=alpha[i];
    }
    this->iterations.push_back(r);
  }
  void print(){
    std::list<Iteration*>::iterator l;
    for(l=this->iterations.begin();l!=this->iterations.end();l++){
      printf("lambda=%f grad=%f %f ...\n",
	     (*l)->lambda,(*l)->grad,(*l)->alpha[0]);
    }
  }
};

Results* join_clusters2_restart
(double*,
 SymNoDiag*,
 unsigned int,
 double,double,double,double,double,
 int,int,int,int,int,int);

