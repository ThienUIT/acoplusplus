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
      File:    TSP.h
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

#ifndef _THOP_H_
#define _THOP_H_


#define RRR            6378.388
#ifndef PI             /* as in stroustrup */
#define PI             3.14159265358979323846
#endif

struct point {
    double x;
    double y;
};

struct item { 
    long int profit;
    long int weight;
    long int id_city;
};

struct problem {
    char          knapsack_data_type[LINE_BUF_LEN];              /* knapsack data type */
    char          edge_weight_type[LINE_BUF_LEN];                /* selfexplanatory */
    long int      optimum;                /* optimal total profit if known, otherwise a bound */
    long int      n;                      /* number of cities */
    long int      m;                      /* number of items */
    long int      capacity_of_knapsack;   /* capacity of knapsack  */
    double        max_time;               /* maximum time  */
    double        min_speed;              /* minimum speed of the thief */
    double        max_speed;              /* maximum speed of the thief */
    long int      n_near;                 /* number of nearest neighbors */
    struct point  *nodeptr;               /* array of structs containing coordinates of nodes */
    struct item   *itemptr;               /* array of structs containing item data */
    long int      **distance;             /* distance matrix: distance[i][j] gives distance between city i und j */
    long int      **nn_list;              /* nearest neighbor list; contains for each node i a sorted list of n_near nearest neighbors */
    long int      UB;                     /* objective value of the optimal solution of the fractional knapsack problem */
};

typedef struct {
    double dist;
} EdgeInfo;

extern struct problem instance;

extern long int (*distance)(long int, long int);  

long int round_distance(long int i, long int j);

long int ceil_distance(long int i, long int j);

long int geo_distance(long int i, long int j);

long int att_distance(long int i, long int j);

long int** compute_distances(void);

long int** compute_nn_lists(void);

long int compute_fitness(long int *t, char *visited, long int t_size, char *p);

void build_route_edge_info(const long int *t, long int t_size, EdgeInfo *edges);

inline double compute_travel_time_on_edge(double dist, double vmax, double vmin, double nu, long W);

double travel_time_increment_when_adding_weight(const EdgeInfo *edges, long t_size, const long *W_edge, long idx, long w, double vmax, double vmin, double nu);

double travel_time_decrement_when_removing_weight(const EdgeInfo *edges, long t_size, const long *W_edge, long idx, long w, double vmax, double vmin, double nu);

long greedy_packing_by_marginal_profit_per_time(const long *t, long t_size, const char *visited, char *p);

long local_search_swap_items_to_improve_profit(const long *t, long t_size, const char *visited, char *p, long current_profit) ;

long int compute_fitness_hybrid(long *t, char *visited, long t_size, char *p);

long greedy_packing_with_lookahead(const long *t, long t_size, const char *visited, char *p);

double simulate_time(const long int *t, long int t_size,
                            long int *weight_accumulated,
                            const char *plan);

int try_drop_or_swap_items(const long *t, long t_size, const char *visited,
                                  char *p, long *profit);

#endif
