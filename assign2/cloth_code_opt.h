#define MAX(a,b) ( (a) > (b) ? (a) : (b))
#define MIN(a,b) ( (a) < (b) ? (a) : (b))

void initMatrix(int n, double mass, double fcon, 
		int delta, double grav, double sep, double rball,
		double offset, double dt, 
		double** x,   double** y,   double** z,
		double** cpx, double** cpy, double** cpz,
		double **fx, double **fy, double **fz,
		double **vx, double **vy, double **vz,
		double **oldfx, double **oldfy, double **oldfz);

void loopcode(int n, double mass, double fcon,
	      int delta, double grav, double sep, double rball,
	      double xball, double yball, double zball,
	      double dt, double *x, double *y, double *z,
	      double *fx, double *fy, double *fz,
	      double *vx, double *vy, double *vz,
	      double *oldfx, double *oldfy, double *oldfz, 
	      double *pe, double *ke,double *te);

double eval_pef(int n, int delta, double grav, double sep,
                double fcon, double *x, double *y, double *z,		
                double *fx, double *fy, double *fz);
