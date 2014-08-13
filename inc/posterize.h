//posterize.h
#ifndef __POSTERIZE_H__
#define __POSTERIZE_H__ 1

/* posterize.h
 * Kmeans implementation by:
 * Ethan Brodsky
 * October 2011
 *
 * please refer posterize.cpp for more detail
 */	
#include <string>
#include "util.h"

void fail(std::string);
double calc_distance(int, double*, double*);
void calc_all_distances(int, int, int, double*, double*, double*);
double calc_total_distance(int, int, int, double*, double*, int*);
void choose_all_clusters_from_distances(int, int, int, double*, int*);
void calc_cluster_centroids(int, int, int, double*, int*, double*);
void get_cluster_member_count(int, int, int*, int*);
void update_delta_score_table(int, int, int, double*, int*, double*, int*, double*, int);
void  perform_move(int, int, int, double*, int*, double*, int*, int, int);
void cluster_diag(int, int, int, double*, int*, double*);
void copy_assignment_array(int, int*, int*);
int assignment_change_count(int, int[], int[]);
void kmeans(int, double*, int, int, double*, int*, uint);
void posterize(ImageMatrix*, uint, uint);

#endif