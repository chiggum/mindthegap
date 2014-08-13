//posterize.cpp
#include "posterize.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#define sqr(x) ((x)*(x))

#define MAX_CLUSTERS 100

#define BIG_double (INFINITY)

void fail(std::string str)
{
  std::cout << str << std::endl;
  exit(EXIT_FAILURE);
}

//Calculates square distance between two points
double calc_distance(int dim, double *p1, double *p2)
{	
	double distance_sq_sum = 0;

	for (int ii = 0; ii < dim; ii++)
	  distance_sq_sum += sqr(p1[ii] - p2[ii]);

	return distance_sq_sum;
}	

//Calculates all distances between X and centroids and store it in distance_output
void calc_all_distances(int dim, int n, int k, double *X, double *centroid, double *distance_output)
{
	for (int ii = 0; ii < n; ii++) // for each point
	{
		for (int jj = 0; jj < k; jj++) // for each cluster
		{
			// calculate distance between point and cluster centroid
			distance_output[ii*k + jj] = calc_distance(dim, &X[ii*dim], &centroid[jj*dim]);
		}
	}
}

//Calculates total distance between X and assigned cluster centroid
double calc_total_distance(int dim, int n, int k, double *X, double *centroids, int *cluster_assignment_index)
// NOTE: a point with cluster assignment -1 is ignored
{
  double tot_D = 0;
  
  // for every point
  for (int ii = 0; ii < n; ii++)
    {
      // which cluster is it in?
      int active_cluster = cluster_assignment_index[ii];
      
      // sum distance
      if (active_cluster != -1)
        tot_D += calc_distance(dim, &X[ii*dim], &centroids[active_cluster*dim]);
    }
    
  return tot_D;
}

//Assigns new clusters based on distance_array
void choose_all_clusters_from_distances(int dim, int n, int k, double *distance_array, int *cluster_assignment_index)
{
  // for each point
  for (int ii = 0; ii < n; ii++)
    {
      int best_index = -1;
      double closest_distance = BIG_double;
      
      // for each cluster
      for (int jj = 0; jj < k; jj++)
        {
          // distance between point and cluster centroid
          double cur_distance = distance_array[ii*k + jj];
          if (cur_distance < closest_distance)
            {
              best_index = jj;
              closest_distance = cur_distance;
            }
        }

      // record in array
      cluster_assignment_index[ii] = best_index;
    }
}

//Calculates cluster centroids
void calc_cluster_centroids(int dim, int n, int k, double *X, int *cluster_assignment_index, double *new_cluster_centroid)
{
  int cluster_member_count[MAX_CLUSTERS];

 // initialize cluster centroid coordinate sums to zero
  for (int ii = 0; ii < k; ii++) 
    {
      cluster_member_count[ii] = 0;
      
      for (int jj = 0; jj < dim; jj++)
        new_cluster_centroid[ii*dim + jj] = 0;
   }

 // sum all points
 // for every point
  for (int ii = 0; ii < n; ii++)
    {
     // which cluster is it in?
      int active_cluster = cluster_assignment_index[ii];

     // update count of members in that cluster
      cluster_member_count[active_cluster]++;
      
     // sum point coordinates for finding centroid
      for (int jj = 0; jj < dim; jj++)
        new_cluster_centroid[active_cluster*dim + jj] += X[ii*dim + jj];
    }
    
    
 // now divide each coordinate sum by number of members to find mean/centroid
 // for each cluster
  for (int ii = 0; ii < k; ii++) 
    {          
     // for each dimension
      for (int jj = 0; jj < dim; jj++)
        new_cluster_centroid[ii*dim + jj] /= cluster_member_count[ii];  /// XXXX will divide by zero here for any empty clusters!

    }
}

//Get Number of members in the cluster
void get_cluster_member_count(int n, int k, int *cluster_assignment_index, int *cluster_member_count)
{
 // initialize cluster member counts
  for (int ii = 0; ii < k; ii++) 
    cluster_member_count[ii] = 0;

 // count members of each cluster    
  for (int ii = 0; ii < n; ii++)
    cluster_member_count[cluster_assignment_index[ii]]++;
}

//Updates delta score table
void update_delta_score_table(int dim, int n, int k, double *X, int *cluster_assignment_cur, double *cluster_centroid, int *cluster_member_count, double *point_move_score_table, int cc)
{
 // for every point (both in and not in the cluster)
  for (int ii = 0; ii < n; ii++)
    {
      double dist_sum = 0;
      for (int kk = 0; kk < dim; kk++)
        {
          double axis_dist = X[ii*dim + kk] - cluster_centroid[cc*dim + kk]; 
          dist_sum += sqr(axis_dist);
        }
        
      double mult = ((double)cluster_member_count[cc] / (cluster_member_count[cc] + ((cluster_assignment_cur[ii]==cc) ? -1 : +1)));

      point_move_score_table[ii*dim + cc] = dist_sum * mult;
    }
}
  
//Perform a move  
void  perform_move(int dim, int n, int k, double *X, int *cluster_assignment, double *cluster_centroid, int *cluster_member_count, int move_point, int move_target_cluster)
{
  int cluster_old = cluster_assignment[move_point];
  int cluster_new = move_target_cluster;

 // update cluster assignment array
  cluster_assignment[move_point] = cluster_new;
  
 // update cluster count array
  cluster_member_count[cluster_old]--;
  cluster_member_count[cluster_new]++;
      
 // update centroid array
  for (int ii = 0; ii < dim; ii++)
    {
      cluster_centroid[cluster_old*dim + ii] -= (X[move_point*dim + ii] - cluster_centroid[cluster_old*dim + ii]) / cluster_member_count[cluster_old];
      cluster_centroid[cluster_new*dim + ii] += (X[move_point*dim + ii] - cluster_centroid[cluster_new*dim + ii]) / cluster_member_count[cluster_new];
    }
}  

//Forms a cluster diagram
void cluster_diag(int dim, int n, int k, double *X, int *cluster_assignment_index, double *cluster_centroid)
{
  int cluster_member_count[MAX_CLUSTERS];
  
  get_cluster_member_count(n, k, cluster_assignment_index, cluster_member_count);
   
}

//Copies an array from src to tgt
void copy_assignment_array(int n, int *src, int *tgt)
{
  for (int ii = 0; ii < n; ii++)
    tgt[ii] = src[ii];
}
  
//Calculates number of changes in the assignement of clusters
int assignment_change_count(int n, int a[], int b[])
{
  int change_count = 0;

  for (int ii = 0; ii < n; ii++)
    if (a[ii] != b[ii])
      change_count++;
      
  return change_count;
}

//Standard K Means clustering algorithm
void kmeans(
            int  dim,		                       // dimension of data 
            double *X,                         // pointer to data
            int   n,                           // number of elements
            int   k,                           // number of clusters
            double *cluster_centroid,          // initial cluster centroids
            int   *cluster_assignment_final,   // output
            uint maxIters
           )
  {
    double *dist                    = (double *)malloc(sizeof(double) * n * k);
    int   *cluster_assignment_cur  = (int *)malloc(sizeof(int) * n);
    int   *cluster_assignment_prev = (int *)malloc(sizeof(int) * n);
    double *point_move_score        = (double *)malloc(sizeof(double) * n * k);
    
    
    if (!dist || !cluster_assignment_cur || !cluster_assignment_prev || !point_move_score)
      fail("Error allocating dist arrays");
    
   // initial setup  
    calc_all_distances(dim, n, k, X, cluster_centroid, dist);
    choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
    copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);

   // BATCH UPDATE
    double prev_totD = BIG_double;
    int batch_iteration = 0;
    while (batch_iteration < (int)maxIters)
      {        
        // update cluster centroids
         calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

        // see if we've failed to improve
         double totD = calc_total_distance(dim, n, k, X, cluster_centroid, cluster_assignment_cur);
         if (totD > prev_totD)
           // failed to improve - currently solution worse than previous
           {
            // restore old assignments
             copy_assignment_array(n, cluster_assignment_prev, cluster_assignment_cur);
             
            // recalc centroids
             calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);
                         
            // done with this phase
             break;
           }
           
        // save previous step
         copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);
         
        // move all points to nearest cluster
         calc_all_distances(dim, n, k, X, cluster_centroid, dist);
         choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
         
         int change_count = assignment_change_count(n, cluster_assignment_cur, cluster_assignment_prev);
         
         
        // done with this phase if nothing has changed
         if (change_count == 0)
           {
             break;
           }

         prev_totD = totD;
                        
         batch_iteration++;
      }

	cluster_diag(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

   // write to output array
    copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_final);    
    
    free(dist);
    free(cluster_assignment_cur);
    free(cluster_assignment_prev);
    free(point_move_score);
  }    

/*
 * posterize.
 * 
 */
void posterize(ImageMatrix *m, uint numClusters, uint max_iters)
{
  uint h, w;
  h = m->height;
  w = m->width;

  //allocate memory
  double *src = new double[3*h*w];
  double *centroids = new double[3*numClusters];
  int *cluster = new int[h*w];

  //flatten input ImageMatrix pixMap and store in src
  for(uint i = 0; i < h; ++i)
  {
      for(uint j = 0; j < w; ++j)
      {
          src[3 * (i * w + j)+0] = (double)(uint)m->pixMap[i][j].r;
          src[3 * (i * w + j)+1] = (double)(uint)m->pixMap[i][j].g;
          src[3 * (i * w + j)+2] = (double)(uint)m->pixMap[i][j].b;
      }
  }

  //centroids initialization using Kmeans++ approach
  for(uint i = 0; i < numClusters; ++i)
  {
    if(i == 0)
    {
      int randNum = rand()%(h*w);
      centroids[3*i+0] = src[3*(randNum) + 0];
      centroids[3*i+1] = src[3*(randNum) + 1];
      centroids[3*i+2] = src[3*(randNum) + 2];
    }
    else
    {
      double* dist = new double[h*w];
      uint max = 0;
      for(uint l = 0; l < h; ++l)
      {
        for(uint m = 0; m < w; ++m)
        {
          for(uint k = 0; k < i; ++k)
          {
            if(k == 0)
            {
              dist[l * w + m] = calc_distance(3,&centroids[3*k], &src[3*(l*w + m)]);
            }
            else
            {
              if(dist[l*w +m] > calc_distance(3,&centroids[3*k], &src[3*(l*w + m)]))
              {
                dist[l*w +m] = calc_distance(3,&centroids[3*k], &src[3*(l*w + m)]);
              }
            }
          }
          if(dist[l*w +m] > dist[max])
          {
            max = l*w+m;
          }
        }
      }
      centroids[3*i+0] = src[3*(max) + 0];
      centroids[3*i+1] = src[3*(max) + 1];
      centroids[3*i+2] = src[3*(max) + 2];
    }
  }

  //perform clustering
  kmeans(3, src, h*w, numClusters, centroids, cluster, max_iters);

  //assign the pixels corresponding to centroids of the clusters to the pixels of the clusters
  for(uint y = 0; y < h; y++)
  {
      for(uint x = 0; x < w; x++)
      {
          m->pixMap[y][x].r = (uint)centroids[3*cluster[(y * w + x)]+0];
          m->pixMap[y][x].g = (uint)centroids[3*cluster[(y * w + x)]+1];
          m->pixMap[y][x].b = (uint)centroids[3*cluster[(y * w + x)]+2];
      }
  }

  delete[] src;
  delete[] centroids;
  delete[] cluster;
}