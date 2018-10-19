#ifndef PSO_H_
#define PSO_H_

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
    double error;
    double *gbest; // should contain DIM elements!!
	int ph;
} pso_result;

typedef struct {
    int dim; // problem dimensionality
    double* lb; // lower range limit
    double* ub; // higher range limit

    int SwarmSize; // swarm size (number of particles)

    int StallIterLimit;
    int MaxIter; // maximum number of iterations
    double TolFun;

    double c1; // cognitive coefficient
    double c2; // social coefficient
    double w_max; // max inertia weight value
    double w_min; // min inertia weight value

    double MinFractionNeighbors; // neighborhood size
} pso_options;

typedef struct {
	double *cx, *cy, *cvx, *cvy;
	int step;
	double cur_fit;
    double *fixed_sol;
    int *is_fixed;
    int *neigh_idx;
} flock_info;

typedef double (*objFcn)(double *, int, void *, int *);

void pso(objFcn fun, void *params,
	       pso_result *solution, pso_options *options);

void settings(pso_options *options, int nvars);

#endif /* PSO_H_ */
