#include <list>
#include <vector>
#include <map>

#define NullCluster ((Cluster*)0)

struct Cluster;

// Doubly-linked list of clusters
typedef std::list<Cluster*> Clusters;

// events are represented as a map, thus e->first is lambda and
// e->second is the pointer to the top cluster of the merge event.
typedef std::multimap<double,Clusters::iterator> Events;
typedef std::pair<double,Clusters::iterator> Event;

struct Cluster {
  int total;
  double v;
  std::vector<int> i;//row indices and cluster size
  double alpha;
  double lambda;
  Cluster *child1;
  Cluster *child2;
  Events::iterator merge_with_next;
};

Cluster* join_clusters(double*,double*,int);
void delete_tree(Cluster*);
Clusters clustermat_l1(double*,int,int);
void print_clusters(Clusters);
void print_path(Cluster*);
Cluster* make_clusters_l1(double*,int);
void add_results(Cluster*,double*,double*,int*,int*);
