/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    TSP.c
      Author:  Thomas Stuetzle
      Purpose: TSP related procedures, distance computation, neighbour lists
      Check:   README and gpl.txt
      Copyright (C) 2002  Thomas Stuetzle
 */

/***************************************************************************

    Program's name: acotsp

    Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the 
    symmetric TSP 

    Copyright (C) 2004  Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    email: stuetzle no@spam ulb.ac.be
    mail address: Universite libre de Bruxelles
                  IRIDIA, CP 194/6
                  Av. F. Roosevelt 50
                  B-1050 Brussels
                  Belgium

 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <cstring>

#include "inout.h"
#include "thop.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"

#define M_PI 3.14159265358979323846264

struct problem instance;

static double dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

long int  (*distance)(long int, long int);  /* function pointer */

/*    
      FUNCTION: the following four functions implement different ways of 
                computing distances for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
 */

long int round_distance (long int i, long int j) 
/*    
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
 */
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd) + 0.5;

    return (long int) r;
}

long int ceil_distance (long int i, long int j) 
/*    
      FUNCTION: compute ceiling distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
 */
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd);

    return (long int)(ceil (r));
}

long int geo_distance (long int i, long int j) 
/*    
      FUNCTION: compute geometric distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: adapted from concorde code
                for the definition of how to compute this distance see TSPLIB
 */
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    long int dd;
    double x1 = instance.nodeptr[i].x, x2 = instance.nodeptr[j].x,
            y1 = instance.nodeptr[i].y, y2 = instance.nodeptr[j].y;

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (y1);
    min = y1 - deg;
    longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (y2);
    min = y2 - deg;
    longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;

}

long int att_distance (long int i, long int j) 
/*    
      FUNCTION: compute ATT distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
 */
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    long int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}



long int** compute_distances(void)
/*    
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
 */
{
    long int     i, j;
    long int     **matrix;

    if((matrix = (long int **) malloc(sizeof(long int) * instance.n * instance.n + sizeof(long int *) * instance.n)) == NULL){
        fprintf(stderr,"Out of memory, exit.");
        exit(1);
    }
    for ( i = 0 ; i < instance.n ; i++ ) {
        matrix[i] = (long int *)(matrix + instance.n) + i * instance.n;
    }

    int max_distance = 0;
    for ( i = 0 ; i < instance.n - 1 ; i++ ) {
        for ( j = 0  ; j < instance.n - 1 ; j++ ) {
            matrix[i][j] = distance(i, j);
            if ( matrix[i][j] > max_distance ) max_distance = matrix[i][j];
        }
    }

    for ( i = 0; i < instance.n ; i++ ) {
        matrix[i][instance.n - 1] = matrix[instance.n - 1][i] = max_distance * ( instance.n - 1 );
    }
    matrix[0][instance.n - 1] = matrix[instance.n - 1][0] = 0;
    matrix[instance.n - 2][instance.n - 1] = matrix[instance.n - 1][instance.n - 2] = 0;

    //matrix[instance.n - 2][instance.n - 2] = matrix[instance.n - 1][instance.n - 1] = 0;

    return matrix;
}


long int** compute_nn_lists( void )
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each city
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
 */
{
    long int i, node, nn;
    long int *distance_vector;
    long int *help_vector;
    long int **m_nnear;

    TRACE ( printf("\n computing nearest neighbor lists, "); )

    nn = MAX(nn_ls,nn_ants);
    if ( nn >= instance.n )
        nn = instance.n - 1;
    DEBUG ( assert( instance.n > nn ); )

    TRACE ( printf("nn = %ld ... \n",nn); )

    if((m_nnear = (long int **) malloc(sizeof(long int) * instance.n * nn + instance.n * sizeof(long int *))) == NULL){
        exit(EXIT_FAILURE);
    }
    distance_vector = (long int *) calloc(instance.n, sizeof(long int));
    help_vector = (long int *) calloc(instance.n, sizeof(long int));

    for ( node = 0 ; node < instance.n ; node++ ) {  /* compute cnd-sets for all node */
        m_nnear[node] = (long int *)(m_nnear + instance.n) + node * nn;

        for ( i = 0 ; i < instance.n ; i++ ) {  /* Copy distances from nodes to the others */
            distance_vector[i] = instance.distance[node][i];
            help_vector[i] = i;
        }
        distance_vector[node] = LONG_MAX;  /* city is not nearest neighbour */
        sort2(distance_vector, help_vector, 0, instance.n - 1);
        for ( i = 0 ; i < nn ; i++ ) {
            m_nnear[node][i] = help_vector[i];
        }
    }
    free(distance_vector);
    free(help_vector);
    TRACE ( printf("\n    .. done\n"); )

    return m_nnear;
}

long int compute_fitness( long int *t, char *visited, long int t_size, char *p ) 
/*    
      FUNCTION: compute the fitness of the ThOP solution generated from tour t
      INPUT:    pointer to tour t and pointer to packing plan p
      OUTPUT:   fitness of the ThOP solution generated from tour t
 */
{

    int i, j, k, l;

    /* for ( i = 0; i <= instance.n; ++i) printf("%d", t[i]); printf("\n"); */
    
    if(t[0] != 0 || t[t_size-1] != 0 || t[t_size-3] != instance.n-2 || t[t_size-2] != instance.n-1) {
        printf("error: compute_fitness\n"); exit(0);
    } 
    
    double par_a, par_b, par_c, par_sum;
    long int prev_city, curr_city;
    double _total_time;
    long int _total_weight, total_weight, total_profit;    
    int violate_max_time;
    
    const double v = ( instance.max_speed - instance.min_speed ) / instance.capacity_of_knapsack;
    
    long int *distance_accumulated = (long int *) malloc ( instance.n * sizeof(long int));

    long int total_distance = 0;

    for ( i = 0 ; i < t_size-1 ; i++ ) {
        distance_accumulated[t[i]] = total_distance;
        total_distance += instance.distance[t[i]][t[i+1]];        
    }

    double *item_vector = (double *) malloc(instance.m * sizeof(double));
    double *help_vector = (double *) malloc(instance.m * sizeof(double));
    
    long int *profit_accumulated = (long int *) malloc( instance.n * sizeof(long int));
    long int *weight_accumulated = (long int *) malloc( instance.n * sizeof(long int));
    
    long int best_packing_plan_profit = 0;
    char *tmp_packing_plan = (char *) malloc(instance.m * sizeof(char));
    
    long int _try;
        
    for( _try = 0; _try < max_packing_tries; _try++) {
        
        for ( i = 0 ; i < instance.n ; i++ ) {
            profit_accumulated[i] = weight_accumulated[i] = 0;
        }
        
        par_a = ran01( &seed );  /* uniform random number between [0.0, 1.0] */
        par_b = ran01( &seed );  /* uniform random number between [0.0, 1.0] */
        par_c = ran01( &seed );  /* uniform random number between [0.0, 1.0] */

        par_sum = (par_a + par_b + par_c); 
        par_a /= par_sum; par_b /= par_sum; par_c /= par_sum;
                
        l = 0;

        for ( j = 0 ; j < instance.m ; j++ ) {
            tmp_packing_plan[j] = 0;
            if (visited[instance.itemptr[j].id_city] == FALSE) continue;
            item_vector[l] = ( 
                                -1.0 * pow(instance.itemptr[j].profit, par_a)
                            ) 
                            / 
                            ( 
                                pow(instance.itemptr[j].weight, par_b) * pow((distance_accumulated[instance.n - 2] - distance_accumulated[instance.itemptr[j].id_city]), par_c)
                            );
            help_vector[l] = j;
            l++;
        }

        sort2_double(item_vector, help_vector, 0, l-1);

        total_weight = 0, total_profit = 0;            
        
        for ( k = 0 ; k < l ; k++ ) {

            j = help_vector[k];
                        
            if ( total_weight + instance.itemptr[j].weight > instance.capacity_of_knapsack ) continue;

            profit_accumulated[instance.itemptr[j].id_city] += instance.itemptr[j].profit;
            weight_accumulated[instance.itemptr[j].id_city] += instance.itemptr[j].weight; 
            
            violate_max_time = FALSE;
            _total_time = _total_weight = 0;
            prev_city = 0;
            for ( i = 1 ; i < t_size; i++ ) {
                curr_city = t[i];
                if ( weight_accumulated[curr_city] == 0 && curr_city != instance.n - 2) continue;
                _total_time += instance.distance[prev_city][curr_city] / ( instance.max_speed - v * _total_weight );    
                if ( _total_time - EPSILON > instance.max_time ) {
                    violate_max_time = TRUE; break;
                }                
                _total_weight += weight_accumulated[curr_city];
                prev_city = curr_city;
            }

            if ( violate_max_time == FALSE) {
                total_profit += instance.itemptr[j].profit;     
                total_weight += instance.itemptr[j].weight;
                tmp_packing_plan[j] = 1;
            }
            else {
                profit_accumulated[instance.itemptr[j].id_city] -= instance.itemptr[j].profit;
                weight_accumulated[instance.itemptr[j].id_city] -= instance.itemptr[j].weight; 
            }
        }
        
        if ( total_profit > best_packing_plan_profit) {
            best_packing_plan_profit = total_profit;
            for ( j = 0 ; j < instance.m ; j++ ) {
                p[j] = tmp_packing_plan[j];
            }            
        }
    }
    
    free(distance_accumulated);
    free(item_vector);
    free(help_vector);
    free(profit_accumulated);
    free(weight_accumulated);
    free(tmp_packing_plan);

    return instance.UB + 1 - best_packing_plan_profit;
}

// -------------------------------------------------------------
// Build edge information (distances along the tour)
// -------------------------------------------------------------
static void build_route_edge_info(const long int *t, long int t_size, EdgeInfo *edges) {
    for (long i = 0; i < t_size - 1; ++i) {
        long u = t[i], v = t[i+1];
        edges[i].dist = (double)instance.distance[u][v];
    }
}

// -------------------------------------------------------------
// Compute travel time for a given edge and load
// -------------------------------------------------------------
static inline double compute_travel_time_on_edge(double dist, double vmax, double vmin, double nu, long W) {
    double sp = vmax - nu * (double)W;
    if (sp < vmin) sp = vmin;  // enforce lower bound
    return dist / sp;
}

// -------------------------------------------------------------
// Δt when adding weight w from edge idx to end
// -------------------------------------------------------------
static double travel_time_increment_when_adding_weight(const EdgeInfo *edges, long t_size,
                             const long *W_edge, long idx, long w,
                             double vmax, double vmin, double nu) {
    double dt = 0.0;
    for (long e = idx; e < t_size - 1; ++e) {
        dt += compute_travel_time_on_edge(edges[e].dist, vmax, vmin, nu, W_edge[e] + w)
            - compute_travel_time_on_edge(edges[e].dist, vmax, vmin, nu, W_edge[e]);
    }
    return dt;
}

// -------------------------------------------------------------
// Δt when removing weight w from edge idx to end
// -------------------------------------------------------------
static double travel_time_decrement_when_removing_weight(const EdgeInfo *edges, long t_size,
                                const long *W_edge, long idx, long w,
                                double vmax, double vmin, double nu) {
    double dt = 0.0;
    for (long e = idx; e < t_size - 1; ++e) {
        dt += compute_travel_time_on_edge(edges[e].dist, vmax, vmin, nu, W_edge[e] - w)
            - compute_travel_time_on_edge(edges[e].dist, vmax, vmin, nu, W_edge[e]);
    }
    return dt;
}

// -------------------------------------------------------------
// Greedy packing: pick items with best profit / Δtime ratio
// -------------------------------------------------------------
long greedy_packing_by_marginal_profit_per_time(const long *t, long t_size, const char *visited, char *p) {
    const double vmax = instance.max_speed, vmin = instance.min_speed;
    const double nu   = (vmax - vmin) / (double)instance.capacity_of_knapsack;

    // Reset plan
    for (int j = 0; j < instance.m; ++j) p[j] = 0;

    // Prepare edges and load
    EdgeInfo *edges = (EdgeInfo *)malloc((t_size-1) * sizeof(EdgeInfo));
    build_route_edge_info(t, t_size, edges);

    long *W_edge = (long *)calloc(t_size-1, sizeof(long)); 
    double total_time = 0.0;  
    long total_w = 0, total_profit = 0;

    // Map cityId -> position in tour
    long *posInTour = (long *)malloc(instance.n * sizeof(long));
    for (long i = 0; i < instance.n; ++i) posInTour[i] = -1;
    for (long i = 0; i < t_size; ++i) posInTour[t[i]] = i;

    // Greedy iteration
    while (1) {
        int best_j = -1;
        double best_score = -1.0;
        double best_dt = 0.0;

        for (int j = 0; j < instance.m; ++j) {
            if (p[j]) continue;  
            int city = instance.itemptr[j].id_city;
            if (!visited[city]) continue;

            long w = instance.itemptr[j].weight;
            if (total_w + w > instance.capacity_of_knapsack) continue;

            long idx = posInTour[city];
            if (idx < 0 || idx >= t_size - 1) continue; 

            double dt = travel_time_increment_when_adding_weight(edges, t_size, W_edge, idx, w, vmax, vmin, nu);

            if (total_time + dt > instance.max_time + 1e-12) continue;

            double score = (double)instance.itemptr[j].profit / (dt > 0 ? dt : 1e-9);
            if (score > best_score) {
                best_score = score; best_j = j; best_dt = dt;
            }
        }

        if (best_j == -1) break;  // no feasible item left

        // Pick best_j
        int city = instance.itemptr[best_j].id_city;
        long wj  = instance.itemptr[best_j].weight;
        long idx = posInTour[city];

        for (long e = idx; e < t_size - 1; ++e) W_edge[e] += wj;
        total_time += best_dt;
        total_w    += wj;
        total_profit += instance.itemptr[best_j].profit;
        p[best_j] = 1;
    }

    free(edges);
    free(W_edge);
    free(posInTour);
    return total_profit;
}

// -------------------------------------------------------------
// Local search: try improving profit by swapping items
// -------------------------------------------------------------
long local_search_swap_items_to_improve_profit(const long *t, long t_size, const char *visited, char *p, long current_profit) {
    const double vmax = instance.max_speed, vmin = instance.min_speed;
    const double nu   = (vmax - vmin) / (double)instance.capacity_of_knapsack;

    EdgeInfo *edges = (EdgeInfo *)malloc((t_size-1) * sizeof(EdgeInfo));
    build_route_edge_info(t, t_size, edges);

    long *posInTour = (long *)malloc(instance.n * sizeof(long));
    for (long i = 0; i < instance.n; ++i) posInTour[i] = -1;
    for (long i = 0; i < t_size; ++i) posInTour[t[i]] = i;

    // Rebuild W_edge and total_time from current plan
    long *W_edge = (long *)calloc(t_size-1, sizeof(long));
    for (int j = 0; j < instance.m; ++j) if (p[j]) {
        long idx = posInTour[ instance.itemptr[j].id_city ];
        long w   = instance.itemptr[j].weight;
        for (long e = idx; e < t_size - 1; ++e) W_edge[e] += w;
    }
    double total_time = 0.0;
    for (long e = 0; e < t_size - 1; ++e) {
        total_time += compute_travel_time_on_edge(edges[e].dist, vmax, vmin, nu, W_edge[e]);
    }

    long total_w = 0;
    for (int j = 0; j < instance.m; ++j) if (p[j]) total_w += instance.itemptr[j].weight;

    int improved = 1;
    while (improved) {
        improved = 0;

        for (int i = 0; i < instance.m; ++i) if (p[i]) {
            long ci = instance.itemptr[i].id_city, wi = instance.itemptr[i].weight, pi = instance.itemptr[i].profit;
            long idx_i = posInTour[ci];

            double dt_remove = travel_time_decrement_when_removing_weight(edges, t_size, W_edge, idx_i, wi, vmax, vmin, nu);

            for (int j = 0; j < instance.m; ++j) if (!p[j] && visited[instance.itemptr[j].id_city]) {
                long wj = instance.itemptr[j].weight, pj = instance.itemptr[j].profit;
                if (total_w - wi + wj > instance.capacity_of_knapsack) continue;

                long cj = instance.itemptr[j].id_city;
                long idx_j = posInTour[cj];

                // Recompute new travel time with i removed and j added
                double new_time;
                {
                    long *Wtmp = (long *)malloc((t_size-1) * sizeof(long));
                    memcpy(Wtmp, W_edge, (t_size-1) * sizeof(long));
                    for (long e = idx_i; e < t_size - 1; ++e) Wtmp[e] -= wi;
                    for (long e = idx_j; e < t_size - 1; ++e) Wtmp[e] += wj;

                    double time_after = 0.0;
                    for (long e = 0; e < t_size - 1; ++e)
                        time_after += compute_travel_time_on_edge(edges[e].dist, vmax, vmin, nu, Wtmp[e]);
                    new_time = time_after;

                    free(Wtmp);
                }

                long new_profit = current_profit - pi + pj;

                if (new_time <= instance.max_time + 1e-12 && new_profit > current_profit) {
                    // Accept move
                    p[i] = 0; p[j] = 1;
                    for (long e = idx_i; e < t_size - 1; ++e) W_edge[e] -= wi;
                    for (long e = idx_j; e < t_size - 1; ++e) W_edge[e] += wj;
                    total_time = new_time;
                    total_w = total_w - wi + wj;
                    current_profit = new_profit;
                    improved = 1;
                    goto NEXT_ITER; 
                }
            }
        }
        NEXT_ITER: ;
    }

    free(edges);
    free(W_edge);
    free(posInTour);
    return current_profit;
}

// -------------------------------------------------------------
// Hybrid fitness function
// with greedy + local search refinements
// -------------------------------------------------------------
long int compute_fitness_hybrid(long *t, char *visited, long t_size, char *p) {
    long int best_profit = 0;
    char *best_plan = (char *)malloc(instance.m * sizeof(char));
    char *tmp_plan  = (char *)malloc(instance.m * sizeof(char));

    // Try multiple randomized heuristics (as in original paper)
    for (int trial = 0; trial < max_packing_tries; trial++) {
        double par_a = ran01(&seed);
        double par_b = ran01(&seed);
        double par_c = ran01(&seed);
        double par_sum = par_a + par_b + par_c;
        par_a /= par_sum; par_b /= par_sum; par_c /= par_sum;

        // Reset temporary plan
        for (int j = 0; j < instance.m; j++) tmp_plan[j] = 0;

        // Compute item scores (author’s formula)
        double *item_vector = (double *)malloc(instance.m * sizeof(double));
        double *help_vector = (double *)malloc(instance.m * sizeof(double));
        int l = 0;
        for (int j = 0; j < instance.m; j++) {
            if (!visited[instance.itemptr[j].id_city]) continue;
            item_vector[l] = (-1.0 * pow(instance.itemptr[j].profit, par_a)) /
                             (pow(instance.itemptr[j].weight, par_b) *
                              pow(1.0 + instance.itemptr[j].id_city, par_c)); 
            help_vector[l] = j;
            l++;
        }

        sort2_double(item_vector, help_vector, 0, l-1);

        // Greedy selection (original-style)
        long int profit = 0, weight = 0;
        for (int k = 0; k < l; k++) {
            int j = help_vector[k];
            if (weight + instance.itemptr[j].weight > instance.capacity_of_knapsack) continue;
            tmp_plan[j] = 1;
            profit += instance.itemptr[j].profit;
            weight += instance.itemptr[j].weight;
        }

        free(item_vector);
        free(help_vector);

        // Refine with strict greedy + swap local search
        profit = greedy_packing_by_marginal_profit_per_time(t, t_size, visited, tmp_plan);
        profit = local_search_swap_items_to_improve_profit(t, t_size, visited, tmp_plan, profit);

        if (profit > best_profit) {
            best_profit = profit;
            memcpy(best_plan, tmp_plan, instance.m * sizeof(char));
        }
    }

    memcpy(p, best_plan, instance.m * sizeof(char));

    free(best_plan);
    free(tmp_plan);

    // Return fitness (smaller = better, consistent with original)
    return instance.UB + 1 - best_profit;
}