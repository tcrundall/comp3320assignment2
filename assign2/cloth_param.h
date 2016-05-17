// Default values
int n=10, delta=2, maxiter=10000;
double sep=1.0, mass=1.0, fcon=10;
double grav=1.0, dt=0.01;
double xball=0.0, yball=0.0, zball=0.0, rball=3.0, offset=0.0;

// Pointers to cloth data structures
double *x,*y,*z,*fx,*fy,*fz,*vx,*vy,*vz,*oldfx,*oldfy,*oldfz;

//OpenGL related stuff
double *cpx,*cpy,*cpz;
int update=3, rendermode=1, paused=0, loop=0; 
double reset_time = 100;

