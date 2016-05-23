#include "cloth_code.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <emmintrin.h>

void initMatrix(int n, double mass, double fcon, 
                int delta, double grav, double sep, double rball,
                double offset,
                double dt, double **x, double **y, double **z,
                double** cpx, double** cpy, double** cpz,
                double **fx, double **fy, double **fz,
                double **vx, double **vy, double **vz,
                double **oldfx, double **oldfy, double **oldfz)
{
  //printf("In Init\n");
  int i,nx,ny;
  
  // Free any existing
  free(*x); free(*y); free(*z);
  free(*cpx); free(*cpy); free(*cpz);
  
  // allocate arrays to hold locations of nodes
  *x=(double*)malloc(n*n*sizeof(double));
  *y=(double*)malloc(n*n*sizeof(double));
  *z=(double*)malloc(n*n*sizeof(double));
  //This is for opengl stuff
  *cpx=(double*)malloc(n*n*sizeof(double));
  *cpy=(double*)malloc(n*n*sizeof(double));
  *cpz=(double*)malloc(n*n*sizeof(double));

  //initialize coordinates of cloth
  for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny++){
      (*x)[n*nx+ny] = nx*sep-(n-1)*sep*0.5+offset;
      (*z)[n*nx+ny] = rball+1;
      (*y)[n*nx+ny] = ny*sep-(n-1)*sep*0.5+offset;
      (*cpx)[n*nx+ny] = 0;
      (*cpz)[n*nx+ny] = 1;
      (*cpy)[n*nx+ny] = 0;
    }
  }

  // Throw away existing arrays
  free(*fx); free(*fy); free(*fz);
  free(*vx); free(*vy); free(*vz);
  free(*oldfx); free(*oldfy); free(*oldfz);
  // Alloc new
  *fx=(double*)malloc(n*n*sizeof(double));
  *fy=(double*)malloc(n*n*sizeof(double));
  *fz=(double*)malloc(n*n*sizeof(double));
  *vx=(double*)malloc(n*n*sizeof(double));
  *vy=(double*)malloc(n*n*sizeof(double));
  *vz=(double*)malloc(n*n*sizeof(double));
  *oldfx=(double*)malloc(n*n*sizeof(double));
  *oldfy=(double*)malloc(n*n*sizeof(double));
  *oldfz=(double*)malloc(n*n*sizeof(double));
  for (i=0;i<n*n;i++){
    (*vx)[i]=0.0;
    (*vy)[i]=0.0;
    (*vz)[i]=0.0;
    (*fx)[i]=0.0;
    (*fy)[i]=0.0;
    (*fz)[i]=0.0;
  }
}

void loopcode(int n, double mass, double fcon,
              int delta, double grav, double sep, double rball,
              double xball, double yball, double zball,
              double dt, double *x, double *y, double *z,
              double *fx, double *fy, double *fz,
              double *vx, double *vy, double *vz,
              double *oldfx, double *oldfy, double *oldfz, 
              double *pe, double *ke,double *te)
{
  //printf("In loop:\n");
  int i,j;
  double rlen,xdiff,ydiff,zdiff,vmag,vmag_i,damp,vdot;
        
  //update position as per MD simulation
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      x[i*n+j] +=dt*(vx[i*n+j]+dt*fx[i*n+j]*0.5);
      oldfx[i*n+j]=fx[i*n+j];
    }
  }
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
    y[i*n+j] +=dt*(vy[i*n+j]+dt*fy[i*n+j]*0.5);
    oldfy[i*n+j]=fy[i*n+j];
    }
  }
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      z[i*n+j] +=dt*(vz[i*n+j]+dt*fz[i*n+j]*0.5);
      oldfz[i*n+j]=fz[i*n+j];
    }
  }

  //apply constraints - push cloth outside of ball, set to zero velocity
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      xdiff = x[i*n+j]-xball;
      ydiff = y[i*n+j]-yball;
      zdiff = z[i*n+j]-zball;
      vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff);
      if (vmag < rball){
        x[i*n+j]=xball+xdiff*rball/vmag;
        y[i*n+j]=yball+ydiff*rball/vmag;
        z[i*n+j]=zball+zdiff*rball/vmag;
        //vx[i*n+j]=0.0;
        //vy[i*n+j]=0.0;
        //vz[i*n+j]=0.0;

        //Set non-tangential component of velocity to zero:
        vdot = (vx[i*n+j]*xdiff + vy[i*n+j]*ydiff + vz[i*n+j]*zdiff)/vmag;

        vx[i*n+j] -= vdot*xdiff/vmag;
        vy[i*n+j] -= vdot*ydiff/vmag;
        vz[i*n+j] -= vdot*zdiff/vmag;
      }
    }
  }

  *pe=eval_pef(n,delta,grav,sep,fcon,x,y,z,fx,fy,fz);

  //Add a damping factor to eventually set velocity to zero
  damp=0.995;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      vx[i*n+j]=(vx[i*n+j]+dt*(fx[i*n+j]+oldfx[i*n+j])*0.5)*damp;
      vy[i*n+j]=(vy[i*n+j]+dt*(fy[i*n+j]+oldfy[i*n+j])*0.5)*damp;
      vz[i*n+j]=(vz[i*n+j]+dt*(fz[i*n+j]+oldfz[i*n+j])*0.5)*damp;
    }
  }
  *ke=0.0;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      *ke += vx[i*n+j]*vx[i*n+j] + vy[i*n+j]*vy[i*n+j] + vz[i*n+j]*vz[i*n+j];
    }
  }
  *ke=*ke/2.0;
  *te=*pe+*ke;

}

void VectorAddSSE(float* __restrict__ a,
                  float* __restrict__ b,
                  float* __restrict__ c,
                  size_t size)
{
  size_t i;
  for (i=0; i<(size/4)*4; i+=4)
  {
    __m128 sse_a = _mm_load_ps(&a[i]);
    __m128 sse_b = _mm_load_ps(&b[i]);

    __m128 sse_c = _mm_add_ps(sse_a, sse_b);

    _mm_store_ps(&c[i], sse_c);
  }
}


double eval_pef(int n, int delta, double grav, double sep,
                double fcon, double *x, double *y, double *z,                
                double *fx, double *fy, double *fz)
{
  //printf("In eval\n");
  double pe,rlen,xdiff,ydiff,zdiff,vmag,vmag_i;
  int nx, ny, dx, dy;
  double rlenArr[2];

  pe = 0.0;
  //loop over particles
  for (ny=0;ny<n;ny++){
    //printf("L.1: %d\n", ny);
    for (nx=0;nx<n;nx++){
      //printf("L.2: %d\n", nx);
      fx[ny*n+nx]=0.0;
      fy[ny*n+nx]=0.0;
      fz[ny*n+nx]=-grav;
      //loop over displacements

      double pe_sum[2] = {0.0, 0.0};
      double fx_sum[2] = {0.0, 0.0};
      double fy_sum[2] = {0.0, 0.0};
      double fz_sum[2] = {0.0, 0.0};
      double dummy[2]  = {0.0, 0.0};

      // Section A
      for (dy=MAX(ny-delta,0);dy<ny;dy++){
        //printf("-- A.1: %d\n", dy);

        // Handles the edge column if unaligned
        if (MAX(nx-delta,0)%2 == 1) 
        {
          dx = MAX(nx-delta,0);
          //printf("Here5\n");
          //printf("-- -- A.2': %d\n", dx);
          // compute reference distance
          rlen=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)))*sep;
          
          // compute actual distance
          xdiff = x[dy*n+dx]-x[ny*n+nx];
          ydiff = y[dy*n+dx]-y[ny*n+nx];
          zdiff = z[dy*n+dx]-z[ny*n+nx];
          vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff); 

          //potential energy and force
          pe += fcon*(vmag-rlen);
          fx[ny*n+nx]+=fcon*xdiff*(vmag-rlen)/vmag;
          fy[ny*n+nx]+=fcon*ydiff*(vmag-rlen)/vmag;
          fz[ny*n+nx]+=fcon*zdiff*(vmag-rlen)/vmag;
        }

        for (dx=MAX(nx-delta,0)+MAX(nx-delta,0)%2; //forces even starting index
            dx<(MIN(nx+delta+1,n)/2) * 2;dx+=2) 
        {
          //printf("-- -- A.2: %d\n", dx);
          // exclude self interaction
          // compute reference distance
          
          rlenArr[0] = sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)))*sep;
          rlenArr[1] = sqrt((double)((nx+1-dx)*(nx+1-dx)+(ny-dy)*(ny-dy)))*sep;
          
          //printf("Here1\n");
          __m128d sse_rlen = _mm_load_pd(&rlenArr[0]);

          __m128d sse_nx = _mm_set1_pd(x[ny*n+nx]);
          __m128d sse_ny = _mm_set1_pd(y[ny*n+nx]);
          __m128d sse_nz = _mm_set1_pd(z[ny*n+nx]);

          __m128d sse_dx = _mm_load_pd(&x[dy*n+dx]);
          __m128d sse_dy = _mm_load_pd(&y[dy*n+dx]);
          __m128d sse_dz = _mm_load_pd(&z[dy*n+dx]);

          __m128d sse_xdiff = _mm_sub_pd(sse_dx, sse_nx);
          __m128d sse_ydiff = _mm_sub_pd(sse_dy, sse_ny);
          __m128d sse_zdiff = _mm_sub_pd(sse_dz, sse_nz);

          //printf("Here2\n");
          __m128d sse_sum = _mm_set1_pd(0.0);

          __m128d sse_mul = _mm_mul_pd(sse_xdiff, sse_xdiff);
          sse_sum = _mm_add_pd(sse_mul, sse_sum);

          sse_mul = _mm_mul_pd(sse_ydiff, sse_ydiff);
          sse_sum = _mm_add_pd(sse_mul, sse_sum);

          sse_mul = _mm_mul_pd(sse_zdiff, sse_zdiff);
          sse_sum = _mm_add_pd(sse_mul, sse_sum);

          __m128d sse_vmag = _mm_sqrt_pd(sse_sum);

          __m128d sse_fcon = _mm_set1_pd(fcon);
          __m128d sse_PE = _mm_sub_pd(sse_vmag, sse_rlen);
          sse_PE = _mm_mul_pd(sse_fcon, sse_PE);
          _mm_store_pd(pe_sum, sse_PE);

          //printf("Here3\n");
          //printf("PE contriubtion is: %f\n", pe_sum[0]);
          pe += pe_sum[0];
          pe += pe_sum[1];

          //printf("Here3.25\n");
          __m128d sse_fx = _mm_mul_pd(sse_PE, sse_rlen);
          sse_fx = _mm_div_pd(sse_fx, sse_vmag);
          //printf("Here3.33\n");
          _mm_store_pd(fx_sum, sse_fx);
          fx[ny*n+nx] += fx_sum[0] + fx_sum[1];
          
          //printf("Here3.5\n");
          __m128d sse_fy = _mm_mul_pd(sse_PE, sse_rlen);
          sse_fy = _mm_div_pd(sse_fy, sse_vmag);
          _mm_store_pd(fy_sum, sse_fy);
          fy[ny*n+nx] += fy_sum[0] + fy_sum[1];

          //printf("Here3.75\n");
          //printf("fz contribution is: %f\n", fz_sum[0]);
          __m128d sse_fz = _mm_mul_pd(sse_PE, sse_rlen);
          sse_fz = _mm_div_pd(sse_fz, sse_vmag);
          _mm_store_pd(fz_sum, sse_fz);
          fy[ny*n+nx] += fz_sum[0] + fz_sum[1];

          //printf("Here4\n");
          //scanf("%s");
        }

        // handles the last column if unaligned
        if (MIN(nx+delta+1,n)%2 == 0)
        {
          //printf("Here5\n");
          //printf("-- -- A.2': %d\n", dx);
          // compute reference distance
          rlen=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)))*sep;
          
          // compute actual distance
          xdiff = x[dy*n+dx]-x[ny*n+nx];
          ydiff = y[dy*n+dx]-y[ny*n+nx];
          zdiff = z[dy*n+dx]-z[ny*n+nx];
          vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff); 

          //potential energy and force
          pe += fcon*(vmag-rlen);
          fx[ny*n+nx]+=fcon*xdiff*(vmag-rlen)/vmag;
          fy[ny*n+nx]+=fcon*ydiff*(vmag-rlen)/vmag;
          fz[ny*n+nx]+=fcon*zdiff*(vmag-rlen)/vmag; 
        }
      }

      // Section B
      for (dy=ny+1;dy<MIN(ny+delta+1,n);dy++){
        //printf("-- B.1: %d\n", dy);
        for (dx=MAX(nx-delta,0);dx<MIN(nx+delta+1,n);dx++){
          //printf("-- -- B.2: %d\n", dx);
          // exclude self interaction
          // compute reference distance
          rlen=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)))*sep;
          // compute actual distance
          xdiff = x[dy*n+dx]-x[ny*n+nx];
          ydiff = y[dy*n+dx]-y[ny*n+nx];
          zdiff = z[dy*n+dx]-z[ny*n+nx];
          vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff); //changed to 1/vmag_i
          vmag_i = 1/vmag;
          //potential energy and force
          pe += fcon*(vmag-rlen);
          fx[ny*n+nx]+=fcon*xdiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
          fy[ny*n+nx]+=fcon*ydiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
          fz[ny*n+nx]+=fcon*zdiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
        }
      }
      // Section C
      dy=ny;
      for (dx=MAX(nx-delta,0);dx<nx;dx++){
        // compute reference distance
        rlen=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)))*sep;
        // compute actual distance
        xdiff = x[dy*n+dx]-x[ny*n+nx];
        ydiff = y[dy*n+dx]-y[ny*n+nx];
        zdiff = z[dy*n+dx]-z[ny*n+nx];
        vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff); //changed to 1/vmag_i
        vmag_i = 1/vmag;
        //potential energy and force
        pe += fcon*(vmag-rlen);
        fx[ny*n+nx]+=fcon*xdiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
        fy[ny*n+nx]+=fcon*ydiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
        fz[ny*n+nx]+=fcon*zdiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
      }
      // Section D
      for (dx=nx+1;dx<MIN(nx+delta+1,n);dx++){
        // compute reference distance
        rlen=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)))*sep;
        // compute actual distance
        xdiff = x[dy*n+dx]-x[ny*n+nx];
        ydiff = y[dy*n+dx]-y[ny*n+nx];
        zdiff = z[dy*n+dx]-z[ny*n+nx];
        vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff); //changed to 1/vmag_i
        vmag_i = 1/vmag;
        //potential energy and force
        pe += fcon*(vmag-rlen);
        fx[ny*n+nx]+=fcon*xdiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
        fy[ny*n+nx]+=fcon*ydiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
        fz[ny*n+nx]+=fcon*zdiff*(vmag-rlen)*vmag_i; //changed '/' to '*'
      }
    }
  }
  return pe;
} 
