#include "cloth_code.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void initMatrix(int n, double mass, double fcon, 
                int delta, double grav, double sep, double rball,
                double offset,
                double dt, double **x, double **y, double **z,
                double** cpx, double** cpy, double** cpz,
                double **fx, double **fy, double **fz,
                double **vx, double **vy, double **vz,
                double **oldfx, double **oldfy, double **oldfz)
{
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

double eval_pef(int n, int delta, double grav, double sep,
                double fcon, double *x, double *y, double *z,                
                double *fx, double *fy, double *fz)
{
  double pe,rlen,xdiff,ydiff,zdiff,vmag,vmag_i;
  int nx, ny, dx, dy;

  pe = 0.0;
  //loop over particles
  for (ny=0;ny<n;ny++){
    for (nx=0;nx<n;nx++){
      fx[ny*n+nx]=0.0;
      fy[ny*n+nx]=0.0;
      fz[ny*n+nx]=-grav;
      //loop over displacements
      // Section A
      for (dy=MAX(ny-delta,0);dy<ny;dy++){
        for (dx=MAX(nx-delta,0);dx<MIN(nx+delta+1,n);dx++){
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
      // Section B
      for (dy=ny+1;dy<MIN(ny+delta+1,n);dy++){
        for (dx=MAX(nx-delta,0);dx<MIN(nx+delta+1,n);dx++){
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
