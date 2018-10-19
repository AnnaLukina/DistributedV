/*
Author: Anna Lukina <lukina.a.u@gmail.com> 
distributed version of adaptive horizon and neighborhood MPC with PSO
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "pso.h"
#include <omp.h>


#define PI 3.14159265358979323846

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const int Num = 7;
const int NN_min = 6;
const int NN_max = 7;
const double init_head = PI/2;
const double head = PI/2;
int NN = 7; // Number of neighbors including itself
//const double RN = 5.0; // Radius of neighborhood
// Number of time stpes
const int Steps = 70;

// Initial range

const double Init_box = 3.0; // 4.0;

const double InitVmin = 0.25;
const double InitVmax = 0.75;

// Minimum distance for collision freedom
const double Dmin = 1;
const double Dmax = 100.0;

int Ph = 1;
const int MaxPh = 3;
// wing span
const double w = 1.0;
// y_opt for upwash
const double d0 = 1.0;
// x_opt is 2w-lambda
const double lambda = 0.5 - PI / 8.0;
// angle of clear view cone
const double angle = PI / 6.0;

// Gaussian params for upwash
const double u_sigma1 = 5.0;
const double u_sigma2 = 5.0;
const double d_sigma1 = 1.0 / 0.3;
const double d_sigma2 = 1.0 / 0.7;

// bound on acceleration w.r.t velocity
double delta = 1;

int currentConf;
int currentRun;
const int Clones = 1;

typedef struct {
  double dis;
  int idx;
} neigh_dis;

static inline double randd(double min, double max) {
  return ((double)rand() / (double)RAND_MAX) * (max - min) + min;
}

double gaussian() {
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
  double mu = 0.0f;
  double sigma = 0.1f;
  if (call == 1) {
    call = !call;
    return (mu + sigma * (double)X2);
  }

  do {
    U1 = -1 + ((double)rand() / RAND_MAX) * 2;
    U2 = -1 + ((double)rand() / RAND_MAX) * 2;
    W = pow(U1, 2) + pow(U2, 2);
  } while (W >= 1 || W == 0);

  mult = sqrt((-2 * log(W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double)X1);
}

static inline double norm(double x, double y) { return sqrt(x * x + y * y); }

static inline double mvnpdf(double x, double y, double a, double b) {
  return exp(-0.5 * (a * x * x + b * y * y));
}

static inline double dot(double x1, double y1, double x2, double y2) {
  return x1 * x2 + y1 * y2;
}

double v_matching(double *vx, double *vy, int num) {
  double sum = 0.0f;
  for (int i = 0; i < num; i++) {
    for (int j = i + 1; j < num; j++) {
      double diff = norm(vx[i] - vx[j], vy[i] - vy[j]) /
                    (norm(vx[i], vy[i]) + norm(vx[j], vy[j]));
      sum += diff * diff;
    }
  }
  return sum;
}

// minimum distance
double cpa(double x1, double y1, double x2, double y2, double vx1, double vy1,
           double vx2, double vy2) {
  double dvx = vx1 - vy1, dvy = vx2 - vy2;
  double dv2 = dot(dvx, dvy, dvx, dvy);
  double cpatime = 0;
  if (dv2 > 1e-8) {
    double wx = x1 - x2, wy = y1 - y2;
    cpatime = -dot(wx, wy, dvx, dvy) / dv2;
  }
  if (cpatime <= 0 || cpatime > 1)
    return Dmin;//INFINITY;
  else
    return norm(x1 - x2 + cpatime * (vx1 - vx2),
                y1 - y2 + cpatime * (vy1 - vy2));
}

//
void trim(double *x, double *y, double mag) {
  if (norm(*x, *y) <= mag)
    return;
  double theta = atan2(*y, *x);
  *x = mag * cos(theta);
  *y = mag * sin(theta);
  //     printf("trim\n");
}

void init(double x[Steps][Num], double y[Steps][Num], double vx[Steps][Num],
          double vy[Steps][Num], int horizon[Steps][Num],
          int neighbors[Steps][Num]) {
  int isCollision, isDivergent;
  //srand(11);
  do {
    isCollision = 0;
    isDivergent = 0;
    for (int i = 0; i < Steps; i++) {
      for (int j = 0; j < Num; j++) {
        if (i != 0) {
          x[i][j] = y[i][j] = vx[i][j] = vy[i][j] = 0;
        } else {
          x[i][j] = randd(0, Init_box) + 100;
          y[i][j] = randd(0, Init_box) + 100;
          vx[i][j] = randd(InitVmin, InitVmax);
          vy[i][j] = 0;//randd(InitVmin, InitVmax);
      }
      }
    }
    for (int bi = 0; bi < Num; bi++) {
      for (int bj = bi + 1; bj < Num; bj++) {
        if (norm(x[0][bi] - x[0][bj], y[0][bi] - y[0][bj]) < Dmin)
          isCollision = 1;
      }
    }
    // check the heading angle
    //for (int bi = 0; bi < Num; bi++) {
    //    if (atan2(vx[0][bi],vy[0][bi]) != init_head)
    //      isDivergent = 1;
    //}
  } while (isCollision || isDivergent);
  for (int t = 0; t < Steps; t++) {
    for (int j = 0; j < Num; j++) {
      horizon[t][j] = 0;
    }
  }
  for (int t = 0; t < Steps; t++) {
    for (int j = 0; j < Num; j++) {
      neighbors[t][j] = 0;
    }
  }
}

double calculateFitness(double *nextVX, double *nextVY, double *nextX,
                        double *nextY, double *firstVX, double *firstVY,
                        void *params, int forTheRecord, int num) {
  flock_info *info = (flock_info *)params;

  double obstacle = 0.0, benefit = 0.0, ca = 0.0, co = 0.0, heading = 0.0;
  double blocks[num][2];
  double px = 0.0, py = 0.0, k = 0.0, A = 0.0, B = 0.0, C = 0.0, side = 0.0,
         h_dis = 0.0, v_dis = 0.0, sm = 0.0, dot_prod = 0.0, ub_j = 0.0;
  //#pragma omp parallel for
  for (int i = 0; i < num; i++) {
    memset(blocks, 0, sizeof(double) * num * 2);
    A = nextVX[i];
    B = nextVY[i];
    C = -nextVY[i] * nextY[i] - nextVX[i] * nextX[i];
    ub_j = 0;
    // compute center of mass
    double cmX = 0, cmY = 0;
    for (int k = 0; k < num; k++) {
    	cmX += nextX[k];
    	cmY += nextY[k];
    }
    cmX /= num; cmY /= num;
    // check the heading angle
	  if (atan2(nextVX[i],nextVY[i]) != head) {
	      heading += fabs(atan2(nextVX[i],nextVY[i]) - head);
	  }
    //#pragma omp parallel for
    for (int j = 0; j < num; j++) {
      if (j != i) {
        if (forTheRecord == 0) {
          if (cpa(info->cx[i], info->cy[i], info->cx[j], info->cy[j],
                  firstVX[i], firstVY[i], firstVX[j], firstVY[j]) < Dmin) {
            ca = INFINITY;
            break;
          }
        }
        if (nextVX[i] == 0.0) {
          px = nextX[j];
          py = nextY[i];
        } else if (nextVY[i] == 0.0) {
          px = nextX[i];
          py = nextY[j];
        } else {
          k = -nextVX[i] / nextVY[i];
          px = (k * nextX[i] + nextX[j] / k + nextY[j] - nextY[i]) /
               (k + 1.0 / k);
          py = -1.0 / k * (px - nextX[j]) + nextY[j];
        }

        side = A * nextX[j] + B * nextY[j] + C;
        h_dis = norm(px - nextX[i], py - nextY[i]);
        v_dis = fabs(side) / norm(A, B);

        if (side >= 0.0 && (h_dis < w || (h_dis - w) / v_dis < tan(angle))) {
          blocks[j][0] = atan(v_dis / (h_dis + w));
          blocks[j][1] = atan2(v_dis, h_dis - w);
          if (blocks[j][0] < PI / 2.0 - angle)
            blocks[j][0] = PI / 2.0 - angle;
          if (blocks[j][1] > PI / 2.0 + angle)
            blocks[j][1] = PI / 2.0 + angle;
          obstacle += (blocks[j][1] - blocks[j][0]) / (angle);
        }

        sm = erf((h_dis - (w - lambda)) * sqrt(2.0) * 8.0);
        dot_prod = (nextVX[i] * nextVX[j] + nextVY[i] * nextVY[j]) /
                   (norm(nextVX[i], nextVY[i]) * norm(nextVX[j], nextVY[j]));
        if (side > 0.0 && h_dis >= w - lambda)
          ub_j += dot_prod * sm * mvnpdf(h_dis - (2.0 * w - lambda), v_dis - d0,
                                         u_sigma1, u_sigma2);
        else if (side >= 0.0 && h_dis < w - lambda)
          ub_j += sm * mvnpdf(h_dis, v_dis, d_sigma1, d_sigma2);
      }
    }

    if (ub_j < 1.0) {
      benefit += ub_j;
    } else {
      benefit += 1.0;
    }

    // benefit += MIN(ub_j,1);
    // cohesion
    co += norm(cmX - nextX[i], cmY - nextY[i]);
  }
  if (co > Dmax) {
  	co -= Dmax;
  } else {
  	co = 0;
  }
  //printf("PSO:%f\n", pow(v_matching(nextVX, nextVY, num), 2) + pow(obstacle, 2) + pow(num - 1.0 - benefit, 2) + ca);
  return pow(v_matching(nextVX, nextVY, num), 2) + pow(obstacle, 2) + pow(num - 1.0 - benefit, 2) + ca + pow(heading/num, 2);
  //return pow(obstacle, 2) + pow(num - 1.0 - benefit, 2) + ca*ca;
}

double flock_fit(double *va, int dim, void *params, int *ph) {
  flock_info *info = (flock_info *)params;
  int nvars = dim / 2 / Ph;
  int num = NN;
  double curFitness;
  double bestFitness = DBL_MAX;
  double *nextX, *nextY, *nextVX, *nextVY, *prevX, *prevY, *prevVX, *prevVY,
      *firstVX, *firstVY;
  nextX = malloc(num * sizeof(double));
  nextY = malloc(num * sizeof(double));
  nextVX = malloc(num * sizeof(double));
  nextVY = malloc(num * sizeof(double));
  prevX = malloc(num * sizeof(double));
  prevY = malloc(num * sizeof(double));
  prevVX = malloc(num * sizeof(double));
  prevVY = malloc(num * sizeof(double));
  firstVX = malloc(num * sizeof(double));
  firstVY = malloc(num * sizeof(double));
  for (int j = 0; j < Ph; j++) {
    for (int i = 0, k = 0; i < num; i++) {
      int b = info->neigh_idx[i];
      if (j == 0) {
        if (info->is_fixed[b]) {
          nextVX[i] = info->cvx[i] + info->fixed_sol[b] * cos(info->fixed_sol[b + MaxPh * Num]);
          nextVY[i] = info->cvy[i] + info->fixed_sol[b] * sin(info->fixed_sol[b + MaxPh * Num]);
        } else {
          nextVX[i] = info->cvx[i] + va[k] * cos(va[k + Ph * nvars]);
          nextVY[i] = info->cvy[i] + va[k] * sin(va[k + Ph * nvars]);
          k++;
        }
        nextX[i] = info->cx[i] + nextVX[i];
        nextY[i] = info->cy[i] + nextVY[i];
        firstVX[i] = nextVX[i];
        firstVY[i] = nextVY[i];
      } else {
        if (info->is_fixed[b]) {
          nextVX[i] = prevVX[i] + info->fixed_sol[b + j * Num] * cos(info->fixed_sol[b + (j + MaxPh) * Num]);
          nextVY[i] = prevVY[i] + info->fixed_sol[b + j * Num] * sin(info->fixed_sol[b + (j + MaxPh) * Num]);
        } else {
          nextVX[i] = prevVX[i] + va[k + j * nvars] * cos(va[k + (Ph + j) * nvars]);
          nextVY[i] = prevVY[i] + va[k + j * nvars] * sin(va[k + (Ph + j) * nvars]);
          k++;
        }
        nextX[i] = prevX[i] + nextVX[i];
        nextY[i] = prevY[i] + nextVY[i];
      }
      prevVX[i] = nextVX[i];
      prevVY[i] = nextVY[i];
      prevX[i] = nextX[i];
      prevY[i] = nextY[i];
    }
    curFitness = calculateFitness(nextVX, nextVY, nextX, nextY, firstVX,
                                  firstVY, params, 0, num);
    if (curFitness <= bestFitness) {
      *ph = j + 1;
      bestFitness = curFitness;
    }
  }
  free(nextX);
  free(nextY);
  free(nextVX);
  free(nextVY);
  free(prevX);
  free(prevY);
  free(prevVX);
  free(prevVY);
  free(firstVX);
  free(firstVY);
  return bestFitness;
}

void disp(int t, double x[Steps][Num], double y[Steps][Num],
          double vx[Steps][Num], double vy[Steps][Num]) {
  for (int i = 0; i < Num; i++) {
    printf("%.4f\t%.4f\t%.4f\t%.4f\n", x[t][i], y[t][i], vx[t][i], vy[t][i]);
  }
}

char *concat(char *s1, char *s2) {
  char *result = malloc(strlen(s1) + strlen(s2) + 1);
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}

void save(double x[Steps][Num], double y[Steps][Num], double vx[Steps][Num],
          double vy[Steps][Num], int horizon[Steps][Num],
          int neighbors[Steps][Num]) {
  FILE *fp;

  char numstr[21];
  sprintf(numstr, "%i", currentConf);
  char numstr2[21];
  sprintf(numstr2, "%i", currentRun);
  fp = fopen(
      concat(concat(concat(concat("config", numstr), "_"), numstr2), ".txt"),
      "w");

  for (int t = 0; t < Steps; t++) {
    for (int i = 0; i < Num; i++) {
      fprintf(fp, "%f\t", x[t][i]);
    }
    for (int i = 0; i < Num; i++) {
      fprintf(fp, "%f\t", y[t][i]);
    }
    for (int i = 0; i < Num; i++) {
      fprintf(fp, "%f\t", vx[t][i]);
    }
    for (int i = 0; i < Num; i++) {
      fprintf(fp, "%f\t", vy[t][i]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  FILE *fp2;
  fp2 = fopen(
      concat(concat(concat(concat("horizon", numstr), "_"), numstr2), ".txt"),
      "w");
  for (int t = 0; t < Steps; t++) {
    for (int i = 0; i < Num; i++) {
      fprintf(fp2, "%d\t", horizon[t][i]);
    }
    fprintf(fp2, "\n");
  }
  fclose(fp2);

  FILE *fp3;
  fp3 = fopen(
      concat(concat(concat(concat("neighbors", numstr), "_"), numstr2), ".txt"),
      "w");
  for (int t = 0; t < Steps; t++) {
    for (int i = 0; i < Num; i++) {
      fprintf(fp3, "%d\t", neighbors[t][i]);
    }
    fprintf(fp3, "\n");
  }
  fclose(fp3);
}

void saveFitness(int timestep, double fitness) {
  //  printf("fitness at step %i: %.15f\n",timestep, fitness);
  char numstr[21];
  sprintf(numstr, "%i", currentConf);
  char *filename = concat(concat("fit", numstr), ".txt");
  FILE *fp;
  fp =
      fopen(filename,
            "a");
  if (timestep == 0) {
    if (currentRun == 0) {
      fprintf(fp, "%i;", 0);
    } else {
      fprintf(fp, "\n%i;", currentRun);
    }
  }
  fprintf(fp, "%f;", fitness);
  fclose(fp);
}

/*p = a
void example(int **p){
  *p = malloc(.....)
}

int *a;
example(&a);
*a = 4;*/

void get_neigh(int b, double cur_x[Num], double cur_y[Num], double cur_vx[Num],
               double cur_vy[Num], double **neigh_x, double **neigh_y,
               double **neigh_vx, double **neigh_vy, int *neigh,
               int *neigh_idx, int is_fixed[Num]) {
  neigh_dis cur_dis[Num];
  int num_of_neigh = 0;
  for (int i = 0; i < Num; i++) {
    cur_dis[i].dis = norm(cur_x[b] - cur_x[i], cur_y[b] - cur_y[i]);
    cur_dis[i].idx = i;
  }
  neigh_dis temp;
  for (int i = 0; i < Num - 1; i++) {
    for (int j = 0; j < Num - 1 - i; j++) {
      if (cur_dis[j].dis > cur_dis[j + 1].dis) {
        temp = cur_dis[j];
        cur_dis[j] = cur_dis[j + 1];
        cur_dis[j + 1] = temp;
      }
    }
  }
  *neigh_x = (double *)malloc(NN * sizeof(double));
  *neigh_y = (double *)malloc(NN * sizeof(double));
  *neigh_vx = (double *)malloc(NN * sizeof(double));
  *neigh_vy = (double *)malloc(NN * sizeof(double));
  for (int i = 0; i < NN; i++) {
    (*neigh_x)[i] = cur_x[cur_dis[i].idx];
    (*neigh_y)[i] = cur_y[cur_dis[i].idx];
    (*neigh_vx)[i] = cur_vx[cur_dis[i].idx];
    (*neigh_vy)[i] = cur_vy[cur_dis[i].idx];
    neigh_idx[i] = cur_dis[i].idx;
    if (!is_fixed[neigh_idx[i]]) {
      num_of_neigh++;
    }
  }
  *neigh = num_of_neigh;
}

void flock_pso(double x[Steps][Num], double y[Steps][Num],
               double vx[Steps][Num], double vy[Steps][Num],
               int horizon[Steps][Num], int neighbors[Steps][Num]) {
  objFcn obj_fun[Num];
  for (int b = 0; b < Num; b++) {
  		obj_fun[b] = flock_fit;
	}

  double curFitness = 0;
  int globalFitness = 0;
  double cur_fit_i[Num];
  int level = 0;
  double curFitPrev = DBL_MAX;
  double levelDist = 0;

  int t = 0; // time step = 1
  int v_reached = 0; // flag up if reached V-formation
  int checkpoint = 0; // when to run consensus
  int period = 0;//Num - NN; // how often to run consensus 
  flock_info info[Num];
  flock_info all_info;
  double fixed_sol[Num][2 * Num * MaxPh];
  double temp_sol[Num][2 * Num * MaxPh];
  double self_sol[Num][2 * Num * MaxPh];
  int is_fixed[Num][2 * Num * MaxPh];
  int neigh_idx[Num][Num];
  int num_of_neigh[Num];
  double *neigh_x[Num], *neigh_y[Num], *neigh_vx[Num], *neigh_vy[Num];
  pso_result sol[Num];
  pso_options options[Num];
  double NeighCount[Num];
  int nvars[Num];
  int LNN[Num];
  int last_b = 0;
  //double TotalDisagreement = 0;
 
  for (t = 0; t < Steps - 1 && v_reached == 0; t++) {
    // global consensus
    if (t == 0 || t == checkpoint + period) {
      if (t > 0) {
        // average-consensus protocol here
        /*for (int bi = 0; bi < Num; bi++) {
          NeighCount[bi] = 1; 
          memset(self_sol[bi], 0, 2 * Num * MaxPh * sizeof(double));
          int b_ind = info[bi].neigh_idx[0];
          memcpy(self_sol[bi], info[bi].fixed_sol,
                      2 * Num * MaxPh * sizeof(double));
          printf("consensus: bi=%d\t, neigh_idx: %d\n", bi, info[bi].neigh_idx[0]);
          printf("b_ind=%d\t, self_sol: %.4f\n", b_ind, self_sol[bi][b_ind]);
          // for each bird bi collect solutions suggested by neighbirhoods 
          for (int bj = 0; bj < Num; bj++) {
            if (bi != bj) {
              for (int nj = 0; nj < NN; nj++) {
                // find all the neighborhoods the bird bi participated in
                if (info[bj].neigh_idx[nj] == b_ind) {
                  printf("consensus: bj=%d\t, nj=%d\t, other_sol: %.4f\n", bj, nj, info[bj].fixed_sol[b_ind]);
                  self_sol[bi][b_ind] = info[bj].fixed_sol[b_ind];//0.1 * (info[bi].fixed_sol[b_ind] - info[bj].fixed_sol[b_ind]);
                  self_sol[bi][b_ind + MaxPh * Num] = info[bj].fixed_sol[b_ind + MaxPh * Num];//0.1 * (info[bi].fixed_sol[b_ind + MaxPh * Num] - info[bj].fixed_sol[b_ind + MaxPh * Num]);
                  if (info[bj].fixed_sol[b_ind] != 0 && info[bi].fixed_sol[b_ind] != 0) {
                    NeighCount[b_ind]++;
                  }
                  printf("b_ind=%d\t, sol: %.4f\n", b_ind, self_sol[bi][b_ind]);
                }
              }
            }
          }
          //self_sol[bi][b_ind] -= info[bi].fixed_sol[b_ind];
          //self_sol[bi][b_ind + MaxPh * Num] -= info[bi].fixed_sol[b_ind + MaxPh * Num];
          //self_sol[bi][b_ind] /= NeighCount[b_ind];
          //self_sol[bi][b_ind + MaxPh * Num] /= NeighCount[b_ind];
        }*/
        // reconstract the global state
        //#pragma omp parallel for
        for (int i = 0; i < Num; i++) {
          int j = info[i].neigh_idx[0];
          //TotalDisagreement += norm(self_sol[i][j] / NeighCount[j] - info[i].fixed_sol[j], self_sol[i][j + MaxPh * Num] / NeighCount[j] - info[i].fixed_sol[j + MaxPh * Num]);
          //info[i].fixed_sol[j] = self_sol[i][j];
          //info[i].fixed_sol[j + MaxPh * Num] = self_sol[i][j + MaxPh * Num];
          //printf("reconstruct: b_ind=%d\t, sol: %.4f\n", j, info[i].fixed_sol[j]);
          double ax = info[i].fixed_sol[j] * cos(info[i].fixed_sol[j + MaxPh * Num]);
          double ay = info[i].fixed_sol[j] * sin(info[i].fixed_sol[j + MaxPh * Num]);
          //trim(&ax, &ay, delta * norm(info[i].cvx[0], info[i].cvy[0]));
          //printf("all before: i=%d\t, cx: %.4f\n", 0, info[i].cx[0]);
          info[i].cvx[0] = info[i].cvx[0] + ax;
          info[i].cvy[0] = info[i].cvy[0] + ay;
          //trim(&info[i].cvx[0], &info[i].cvx[0], delta * norm(info[i].cvx[0], info[i].cvy[0]));
          info[i].cx[0] = info[i].cx[0] + info[i].cvx[0];
          info[i].cy[0] = info[i].cy[0] + info[i].cvy[0];
          //printf("all after: i=%d\t, cx: %.4f\n", i, info[i].cx[0]);
          all_info.cvx[i] = info[i].cvx[0];
          all_info.cvy[i] = info[i].cvy[0];
          all_info.cx[i] = info[i].cx[0];
          all_info.cy[i] = info[i].cy[0];
          //printf("temp_x:%f\t, temp_vx:%f\n", all_info.cx[i], all_info.cvx[i]);
        }
      } else {
        //flock_info all_info;
        all_info.cx = x[t];
        all_info.cy = y[t];
        all_info.cvx = vx[t];
        all_info.cvy = vy[t];
        all_info.step = t;
      }
      //for (int i = 0; i < Num; i++) {
        //printf("temp_all_x:%f\t, temp_all_vx:%f\n", all_info.cx[i], all_info.cvx[i]);
      //}
      //if (t == 0) {
	      // compute global cost
	      curFitness = calculateFitness(all_info.cvx, all_info.cvy, all_info.cx, all_info.cy, all_info.cvx, all_info.cvy, &all_info, 1, Num);// + TotalDisagreement;
	  //    globalFitness = 0;
	  //} else { 
	      // find the last achived fitness among all birds
	  //    curFitness = info[last_b].cur_fit;
	  //}
      printf("step: %d, level: %d, curFitness: %f, level_dist: %f, neighbors: %d\n", t, level, curFitness, levelDist, NN);
      if (curFitness <= 1e-1) {
        v_reached = 1;
        for (int i = 0; i < Num; i++) {
          vx[t][i] = all_info.cvx[i];
          vy[t][i] = all_info.cvy[i];
          x[t][i] = all_info.cx[i];
          y[t][i] = all_info.cy[i];

          // save configuration of neighbors
          for (int b = 0; b < Num; b++) {
            FILE *fp;
            char numstr0[21];
            sprintf(numstr0, "%i", b);
            char numstr1[21];
            sprintf(numstr1, "%i", currentConf);
            char numstr2[21];
            sprintf(numstr2, "%i", currentRun);
            fp = fopen(
                concat(concat(concat(concat(concat(concat("config", numstr0), "_"), numstr1), "_"), numstr2), ".txt"),
                "a");
            for (int i = 0; i < NN; i++) {
              int bi = 0;
              while (i != info[b].neigh_idx[bi] && bi < NN) {bi++;}
              fprintf(fp, "%f\t", all_info.cx[i]);
            }
            for (int i = 0; i < NN; i++) {
              int bi = 0;
              while (i != info[b].neigh_idx[bi] && bi < NN) {bi++;}
              fprintf(fp, "%f\t", all_info.cy[i]);
            }
            for (int i = 0; i < NN; i++) {
              int bi = 0;
              while (i != info[b].neigh_idx[bi] && bi < NN) {bi++;}
              fprintf(fp, "%f\t", all_info.cvx[i]);
            }
            for (int i = 0; i < NN; i++) {
              int bi = 0;
              while (i != info[b].neigh_idx[bi] && bi < NN) {bi++;}
              fprintf(fp, "%f\t", all_info.cvy[i]);
            }
            fprintf(fp, "\n");
            fclose(fp);
          }

          printf("end_x:%f\t, end_y:%f\n", x[t][i], y[t][i]);
          free(neigh_x[i]);
          free(neigh_y[i]);
          free(neigh_vx[i]);
          free(neigh_vy[i]);
        }
        break;
      }
      // check if global const improved by levelDist
      if (curFitPrev - curFitness > levelDist) {
        // proceede to the next global state
        for (int i = 0; i < Num; i++) {
          vx[t][i] = all_info.cvx[i];
          vy[t][i] = all_info.cvy[i];
          x[t][i] = all_info.cx[i];
          y[t][i] = all_info.cy[i];

          //printf("next_x:%f\t, next_vx:%f\n", x[t][i], vx[t][i]);
        } 
        // decrease neighborhood size proportionally to cost
        NN = MIN(MAX(NN_min, NN - ceil(1 - curFitness / NN)), NN_max);
        level ++;
        curFitPrev = curFitness; //store fitness from the previous consensus round
        levelDist = curFitPrev / (Steps - level); //Lyapunov dynamic distance function
      } else {
        // increase neighborhood size by 1
        NN = MIN(NN_max, NN + 1);
        // roll back to the previous global state
        for (int i = 0; i < Num; i++) {
          vx[t][i] = all_info.cvx[i];
          vy[t][i] = all_info.cvy[i];
          x[t][i] = all_info.cx[i];
          y[t][i] = all_info.cy[i];

          //printf("back_x:%f\t, back_vx:%f\n", x[t][i], vx[t][i]);
        }
      }
      period = 0;//Num - NN;
      checkpoint = t + 1;
      printf("step: %d, level: %d, curFitness: %f, level_dist: %f, neighbors: %d\n", t, level, curFitness, levelDist, NN);

      // local copies of positions and costs for each bird
      //#pragma omp parallel for
      //flock_info info[Num];
      //#pragma omp parallel for
      for (int b = 0; b < Num; b++) {
        if (t == 0) {
        	LNN[b] = NN;
          memset(fixed_sol[b], 0, 2 * Num * MaxPh * sizeof(double));
        } else {
          memcpy(fixed_sol[b], info[b].fixed_sol, 2 * Num * MaxPh * sizeof(double));
        }
        memset(is_fixed[b], 0, Num * sizeof(int));
        num_of_neigh[b] = 0;
        memset(neigh_idx[b], 0, Num * sizeof(int));
        // find neighborhood for each bird
        get_neigh(b, x[t], y[t], vx[t], vy[t], &neigh_x[b], &neigh_y[b], &neigh_vx[b],
                &neigh_vy[b], &num_of_neigh[b], neigh_idx[b], is_fixed[b]);
        info[b].cx = neigh_x[b];
        info[b].cy = neigh_y[b];
        info[b].cvx = neigh_vx[b];
        info[b].cvy = neigh_vy[b];
        info[b].step = t;
        info[b].fixed_sol = fixed_sol[b];
        info[b].is_fixed = is_fixed[b];
        info[b].neigh_idx = neigh_idx[b];
        if (t == 0) {
        	info[b].cur_fit = curFitPrev;
        }
        //for (int i = 0; i < Num; i++) {
        //  printf("neigh_idx: %d\n", info[b].neigh_idx[i]);
        //  printf("%.4f\t%.4f\t%.4f\t%.4f\n", info[b].cx[i], info[b].cy[i], info[b].cvx[i], info[b].cvy[i]);
        //}
      }
    } else {
      // no change in the global state
      for (int b = 0; b < Num; b++) {
        vx[t][b] = all_info.cvx[b];
        vy[t][b] = all_info.cvy[b];
        x[t][b] = all_info.cx[b];
        y[t][b] = all_info.cy[b];
        if (t == 0) {
          memset(fixed_sol[b], 0, 2 * Num * MaxPh * sizeof(double));
        } else {
          memcpy(fixed_sol[b], info[b].fixed_sol, 2 * Num * MaxPh * sizeof(double));
        }
        memset(is_fixed[b], 0, Num * sizeof(int));
        num_of_neigh[b] = 0;
        memset(neigh_idx[b], 0, Num * sizeof(int));

        //printf("same_x:%f\t, same_vx:%f\n", x[t][b], vx[t][b]);
        info[b].step = t;
        info[b].fixed_sol = fixed_sol[b];
        info[b].is_fixed = is_fixed[b];
        info[b].neigh_idx = neigh_idx[b];
        if (t == 0) {
        	info[b].cur_fit = curFitPrev;
        }
      }
    }
    printf("step: %d, level: %d, curFitness: %f, level_dist: %f, neighbors: %d\n", t, level, curFitness, levelDist, NN);

    // local control independently for each bird
    //#pragma omp parallel for
    for (int b = 0; b < Num; b++) { 
      // compute local cost
      /*if (t < checkpoint + period && t >= checkpoint) {
        //#pragma omp parallel for
        for (int i = 0; i < LNN[b]; i++) {
          int b_ind = info[b].neigh_idx[i];
          printf("neigh_idx: %d\n", b_ind);
          printf("consensused: %.4f\t%.4f\t%.4f\t%.4f\n", info[b].cx[i], info[b].cy[i], info[b].cvx[i], info[b].cvy[i]);
          double ax = info[b].fixed_sol[b_ind] * cos(info[b].fixed_sol[b_ind + MaxPh * Num]);
          double ay = info[b].fixed_sol[b_ind] * sin(info[b].fixed_sol[b_ind + MaxPh * Num]);
          trim(&ax, &ay, delta * norm(info[b].cvx[i], info[b].cvy[i]));
          info[b].cvx[i] = info[b].cvx[i] + ax;
          info[b].cvy[i] = info[b].cvy[i] + ay;
          info[b].cx[i] = info[b].cx[i] + info[b].cvx[i];
          info[b].cy[i] = info[b].cy[i] + info[b].cvy[i];
        }
      }
      for (int i = 0; i < LNN[b]; i++) {
        printf("neigh_idx: %d\n", info[b].neigh_idx[i]);
        printf("continued: %.4f\t%.4f\t%.4f\t%.4f\n", info[b].cx[i], info[b].cy[i], info[b].cvx[i], info[b].cvy[i]);
      }*/
      cur_fit_i[b] = curFitness;//calculateFitness(info[b].cvx, info[b].cvy, info[b].cx, info[b].cy, info[b].cvx, info[b].cvy, &info[b], 1, NN);
      // find the worst fitness among your neighbors
      //cur_fit_i = info[b].cur_fit;
      /*for (int bk = 0; bk < NN; bk++) {
      	int ind = info[b].neigh_idx[bk];
      	if (info[ind].cur_fit > cur_fit_i) {
      		cur_fit_i = info[ind].cur_fit;
      	} 
      }*/
      printf("b=%d\t, info_fit: %f\t, cur_fit_i: %f\n", b, info[b].cur_fit, cur_fit_i[b]);
      //if (cur_fit_i < info[b].cur_fit) {
      //  cur_fit_i = info[b].cur_fit;
      //}

      // save configuration of neighbors
      FILE *fp;
      char numstr0[21];
      sprintf(numstr0, "%i", b);
      char numstr1[21];
      sprintf(numstr1, "%i", currentConf);
      char numstr2[21];
      sprintf(numstr2, "%i", currentRun);
      fp = fopen(
          concat(concat(concat(concat(concat(concat("config", numstr0), "_"), numstr1), "_"), numstr2), ".txt"),
          "a");
      for (int i = 0; i < LNN[b]; i++) {
        int bi = 0;
        while (i != info[b].neigh_idx[bi] && bi < LNN[b]) {bi++;}
        fprintf(fp, "%f\t", info[b].cx[i]);
      }
      for (int i = 0; i < LNN[b]; i++) {
        int bi = 0;
        while (i != info[b].neigh_idx[bi] && bi < LNN[b]) {bi++;}
        fprintf(fp, "%f\t", info[b].cy[i]);
      }
      for (int i = 0; i < LNN[b]; i++) {
        int bi = 0;
        while (i != info[b].neigh_idx[bi] && bi < LNN[b]) {bi++;}
        fprintf(fp, "%f\t", info[b].cvx[i]);
      }
      for (int i = 0; i < LNN[b]; i++) {
        int bi = 0;
        while (i != info[b].neigh_idx[bi] && bi < LNN[b]) {bi++;}
        fprintf(fp, "%f\t", info[b].cvy[i]);
      }
      fprintf(fp, "\n");
      fclose(fp);

      // explore horizons for a local solution
      for (int Ph = 1; Ph <= MaxPh; Ph++) {
        nvars[b] = 2 * LNN[b] * Ph;
        settings(&options[b], nvars[b]);
        sol[b].gbest = malloc(options[b].dim * sizeof(double));
        memset(sol[b].gbest, 0, options[b].dim * sizeof(double));
        options[b].lb = malloc(nvars[b] * sizeof(double));
        options[b].ub = malloc(nvars[b] * sizeof(double));
        for (int i = 0; i < nvars[b]; i++) {
          if (i < Ph * LNN[b]) {
            options[b].lb[i] = 0.0;
            options[b].ub[i] = 0.0;
          } else {
            options[b].lb[i] = 0.0;
            options[b].ub[i] = 2.0 * PI;
          }
        }

        for (int i = 0, j = 0; i < LNN[b]; i++) {
          //printf("b=%d, neigh=%d, fixed=%d\n", b, info[b].neigh_idx[i], info[b].is_fixed[info[b].neigh_idx[i]]);
          if (info[b].is_fixed[info[b].neigh_idx[i]]) continue;
          options[b].ub[j] = delta * norm(info[b].cvx[i], info[b].cvy[i]);
          //printf("options_ub: %f\n", options[b].ub[j]);
          j++;
        }
        for (int p = 1; p < Ph; p++) {
          memcpy(options[b].ub + p * LNN[b], options[b].ub,
                LNN[b] * sizeof(double));
        }

        // run PSO for a neighborhood
        pso(obj_fun[b], &info[b], &sol[b], &options[b]);
        //printf("PSO_fit: %f\n", sol[b].error);
        for (int i = 0; i < options[b].dim; i++) {
        //printf("gbest: %f\n", sol[b].gbest[i]);
      }
        // check if PSO found an improvement or search is exhausted
        if (cur_fit_i[b] - sol[b].error > cur_fit_i[b] / (Steps - t) || Ph == MaxPh) {
          /*if (info[b].cur_fit - sol[b].error <= info[b].cur_fit / (Steps - t) && Ph == MaxPh) {
              //use buffered solution
              memset(temp_sol[b], 0, 2 * Num * MaxPh);
              memcpy(temp_sol[b], info[b].fixed_sol, 2 * Num * MaxPh * sizeof(double));
              for (int i = 0; i <  Num; i++) {
                  int n_idx = info[b].neigh_idx[i];
                for (int p = 0; p < MaxPh-1; p++) {
                  printf("temp_sol:%f\t%f\n", temp_sol[b][i + p * Num], temp_sol[b][i + (p + 1) * Num]);
                  info[b].fixed_sol[n_idx + p * Num] = temp_sol[b][n_idx + (p + 1) * Num];
                  info[b].fixed_sol[n_idx + (p + MaxPh) * Num] = temp_sol[b][n_idx + (p + MaxPh + 1) * Num];
                  printf("n_idx=%d\t, b=%d\t, i=%d\t, buffer: %f\t, %f\n", n_idx, b, i, info[b].fixed_sol[n_idx + p * Num], info[b].fixed_sol[n_idx + (p + MaxPh) * Num]);
                }
              } 
            }*/
          if (sol[b].error < cur_fit_i[b] || t == 0) {
            info[b].cur_fit = sol[b].error;
            printf("b=%d\t, info_fit: %f\t, sol_error: %f\n", b, info[b].cur_fit, sol[b].error);
            // store horizon and neighborhood size for the current bird
            horizon[t][b] = Ph;
            neighbors[t][b] = LNN[b];
            //memset(temp_sol, 0, 2 * Num * MaxPh);
            //memcpy(temp_sol, sol[b].gbest,
            //          2 * NN * Ph * sizeof(double));
            
            // accept the new local solution
            memset(fixed_sol[b], 0, 2 * Num * MaxPh * sizeof(double));
            info[b].fixed_sol = fixed_sol[b];
            //#pragma omp parallel for
            for (int i = 0; i < LNN[b]; i++) {
              int n_idx = info[b].neigh_idx[i];
              for (int p = 0; p < Ph; p++) {
                info[b].fixed_sol[n_idx + p * Num] = sol[b].gbest[i + p * LNN[b]];
                info[b].fixed_sol[n_idx + (p + MaxPh) * Num] = sol[b].gbest[i + (p + Ph) * LNN[b]];
                info[b].is_fixed[n_idx] = 1;
                //printf("n_idx=%d\t, b=%d\t, i=%d\t, before: %f\t, %f\n", n_idx, b, i, info[b].fixed_sol[n_idx + p * Num], info[b].fixed_sol[n_idx + (p + MaxPh) * Num]);
              }
            }
            // current bird found its solution = it's fixed, broadcast
            for (int k = 0; k < LNN[b]; k++){
              int bb = info[b].neigh_idx[k];
              for (int j = 0; j < LNN[b]; j++) {
                if (b == info[bb].neigh_idx[j]) {
                  for (int p = 0; p < Ph; p++) {
                    info[bb].fixed_sol[b + p * Num] = info[b].fixed_sol[b + p * Num];
                    info[bb].fixed_sol[b + (p + MaxPh) * Num] = info[b].fixed_sol[b + (p + MaxPh) * Num];
                    //printf("b=%d\t, bb=%d\t, j=%d\t, before: %f\t, %f\n", b, bb, j, info[bb].fixed_sol[b + p * Num], info[bb].fixed_sol[b + (p + MaxPh) * Num]);
                  }
                  info[bb].is_fixed[b] = 1;
                  //printf("%s\n", "fixed");
                  //break;
                }
              }
            }
            // triger global communication
            if (info[b].cur_fit <= 1e-1) {
            	globalFitness = 1;
            }
		    // decrease neighborhood size locally proportionally to the local cost
		    LNN[b] = MIN(MAX(NN_min, LNN[b] - ceil(1 - info[b].cur_fit / LNN[b])), NN_max);
		    printf("LNN%d= %d\n", b, LNN[b]);
		  }
          free(sol[b].gbest);
          free(options[b].lb);
          free(options[b].ub);
          break;
        } else {
          // increase neighborhood size by 1
	      LNN[b] = MIN(NN_max, LNN[b] + 1);
          free(sol[b].gbest);
          free(options[b].lb);
          free(options[b].ub);
        }
      }
    last_b = b;
    }
  }
}

void help() {
  printf("======= HELP ========\n First argument: [int] unique number of "
         "current CONFIGURATION tried\n Second argument: [int] number of RUNS "
         "performed by this current configuration\n ======================\n");
}

int main(int argc, char *argv[]) {

  if (argc > 1) {
    int helpcalled = 0;
    for (int i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-help") == 0) {
        help();
        helpcalled = 1;
      }
    }
    if (helpcalled == 0) {
      if (argc > 2) {
        currentConf = atoi(argv[1]);
        int runs = atoi(argv[2]);
        for (int k = 0; k < runs;
             k++) { // iterate through the amount of runs per configuration
          currentRun = k;
          srand(currentConf + time(NULL));
          double x[Steps][Num], y[Steps][Num], vx[Steps][Num], vy[Steps][Num];
          int horizon[Steps][Num];
          int neighbors[Steps][Num];
          int loading = 0;
          if (!loading) {
            init(x, y, vx, vy, horizon, neighbors);
            // saveConf(0,x,y,vx,vy);
          }
          disp(0, x, y, vx, vy);
          srand(time(NULL));
          clock_t begin, end;
          double time_spent;
          begin = clock();

          for (int b = 0; b < Num; b++) {
            FILE *fp;
            char numstr0[21];
            sprintf(numstr0, "%i", b);
            char numstr1[21];
            sprintf(numstr1, "%i", currentConf);
            char numstr2[21];
            sprintf(numstr2, "%i", currentRun);
            fp = fopen(
                concat(concat(concat(concat(concat(concat("config", numstr0), "_"), numstr1), "_"), numstr2), ".txt"),
                "w");
            fclose(fp);
          }

          flock_pso(x, y, vx, vy, horizon, neighbors);
          // disp(Steps - 1, x, y, vx, vy);
          save(x, y, vx, vy, horizon, neighbors);

          // saveConf(Steps-1,x,y,vx,vy);

          end = clock();
          time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
          printf("duration: %f\n", time_spent);
          saveFitness(0,time_spent);
        }
      } else
        help();
    }
  } else {
    help();
  }
}
