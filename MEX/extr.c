/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

/************************************************************************/
/*                                                                      */
/* INITIALIZATION OF EXTREMA STRUCTURE                                  */
/*                                                                      */
/************************************************************************/

extrema_t init_extr(int n) {
    extrema_t ex;
    ex.x_min=(double *)malloc(n*sizeof(double));
    ex.x_max=(double *)malloc(n*sizeof(double));
    ex.y_min=(double *)malloc(n*sizeof(double));
    ex.y_max=(double *)malloc(n*sizeof(double));
    ex.x_tmp=(double *)malloc(n*sizeof(double));
    ex.y_tmp=(double *)malloc(n*sizeof(double));
    ex.i_min=(int *)malloc(n*sizeof(int));
    ex.i_max=(int *)malloc(n*sizeof(int));
    return ex;
}


/************************************************************************/
/*                                                                      */
/* DETECTION OF LOCAL EXTREMA                                           */
/*                                                                      */
/************************************************************************/

void extr(double x[],double y[],int n,extrema_t *ex) {
    int cour, iter, count;
    double prev, next;
    ex->n_min=0;
    ex->n_max=0;

    double dmin = 0, dmax = 0;

  /* search for extrema */
    if (ex->is_circular) {
        if (y[0]<=y[n-1] && y[0]<=y[1]) /* local minimum */ {
            ex->x_min[NBSYM]=x[0];
            ex->y_min[NBSYM]=y[0];
            ex->i_min[NBSYM]=0;
            ex->n_min++;
        }
        if (y[0]>=y[n-1] && y[0]>=y[1]) /* local maximum */ {
            ex->x_max[NBSYM]=x[0];
            ex->y_max[NBSYM]=y[0];
            ex->i_max[NBSYM]=0;
            ex->n_max++;
        }
    }
    for(cour=1;cour<(n-1);cour++) {
        if (y[cour]<=y[cour-1] && y[cour]<=y[cour+1]) /* local minimum */ {
            ex->x_min[ex->n_min+NBSYM]=x[cour];
            ex->y_min[ex->n_min+NBSYM]=y[cour];
            ex->i_min[ex->n_min+NBSYM]=cour;

            if (ex->n_min) {
              dmin += x[cour] - ex->x_min[ex->n_min+NBSYM-1];
            }

            ex->n_min++;
        }
        if (y[cour]>=y[cour-1] && y[cour]>=y[cour+1]) /* local maximum */ {
            ex->x_max[ex->n_max+NBSYM]=x[cour];
            ex->y_max[ex->n_max+NBSYM]=y[cour];
            ex->i_max[ex->n_max+NBSYM]=cour;

            if (ex->n_max) {
              dmax += x[cour] - ex->x_max[ex->n_max+NBSYM-1];
            }

            ex->n_max++;
        }
    }

    if (ex->is_circular) {
        if (y[cour]<=y[cour-1] && y[cour]<=y[0]) /* local minimum */ {
            ex->x_min[ex->n_min+NBSYM]=x[cour];
            ex->y_min[ex->n_min+NBSYM]=y[cour];
            ex->i_min[ex->n_min+NBSYM]=cour;

            dmin += x[cour] - ex->x_min[ex->n_min+NBSYM-1];

            ex->n_min++;
        }
        if (y[cour]>=y[cour-1] && y[cour]>=y[0]) /* local maximum */ {
            ex->x_max[ex->n_max+NBSYM]=x[cour];
            ex->y_max[ex->n_max+NBSYM]=y[cour];
            ex->i_max[ex->n_max+NBSYM]=cour;

            dmax += x[cour] - ex->x_max[ex->n_max+NBSYM-1];

            ex->n_max++;
        }
    }

    /* Fill large gaps with intermediate extrema to avoid explosions */
    dmin = COEFF*(dmin / ex->n_min);
    dmax = COEFF*(dmax / ex->n_max);

    ex->x_tmp[0] = ex->x_min[NBSYM];
    ex->y_tmp[0] = ex->y_min[NBSYM];

    count = 1;
    prev = ex->x_tmp[0];
    for(iter=NBSYM+1;iter<ex->n_min+NBSYM;iter++) {
      if ((ex->x_min[iter] - prev) > 2*dmin) {
        next = ex->x_min[iter];
        for(cour=ex->i_min[iter-1];cour<ex->i_min[iter];cour++) {
          if ((x[cour] - prev) > dmin && (next - x[cour]) > dmin) {
            ex->x_tmp[count] = x[cour];
            ex->y_tmp[count] = y[cour];

            prev = x[cour];
            count++;
          }
        }
      }
      ex->x_tmp[count] = ex->x_min[iter];
      ex->y_tmp[count] = ex->y_min[iter];

      prev = ex->x_tmp[count];

      count++;
    }

    ex->n_min = count;
    memcpy(ex->x_min+NBSYM, ex->x_tmp, count*sizeof(double));
    memcpy(ex->y_min+NBSYM, ex->y_tmp, count*sizeof(double));

    ex->x_tmp[0] = ex->x_max[NBSYM];
    ex->y_tmp[0] = ex->y_max[NBSYM];

    count = 1;
    prev = ex->x_tmp[0];
    for(iter=NBSYM+1;iter<ex->n_max+NBSYM;iter++) {
      if ((ex->x_max[iter] - prev) > 2*dmax) {
        next = ex->x_max[iter];
        for(cour=ex->i_max[iter-1];cour<ex->i_max[iter];cour++) {
          if ((x[cour] - prev) > dmax && (next - x[cour]) > dmax) {
            ex->x_tmp[count] = x[cour];
            ex->y_tmp[count] = y[cour];

            prev = x[cour];
            count++;
          }
        }
      }
      ex->x_tmp[count] = ex->x_max[iter];
      ex->y_tmp[count] = ex->y_max[iter];

      prev = ex->x_tmp[count];

      count++;
    }

    ex->n_max = count;
    memcpy(ex->x_max+NBSYM, ex->x_tmp, count*sizeof(double));
    memcpy(ex->y_max+NBSYM, ex->y_tmp, count*sizeof(double));
}

/************************************************************************/
/*                                                                      */
/* EXTRAPOLATION OF EXTREMA TO LIMIT BORDER EFFECTS                     */
/*                                                                      */
/************************************************************************/

void boundary_conditions(double x[],double y[],int n,extrema_t *ex) {
    int cour,nbsym;
    double dx, tmp_x;
    
    nbsym = NBSYM;

    /* reduce the number of symmetrized points if there is not enough extrema */
    //while(ex->n_min < nbsym+1 && ex->n_max < nbsym+1) nbsym--;
    nbsym = SMALLER(SMALLER(ex->n_min, ex->n_max), nbsym);

    if (ex->is_circular) {

      if (nbsym < NBSYM) {
          for(cour=0;cour<ex->n_max;cour++) {
              ex->x_max[nbsym+cour] = ex->x_max[NBSYM+cour];
              ex->y_max[nbsym+cour] = ex->y_max[NBSYM+cour];
          }
          for(cour=0;cour<ex->n_min;cour++) {
              ex->x_min[nbsym+cour] = ex->x_min[NBSYM+cour];
              ex->y_min[nbsym+cour] = ex->y_min[NBSYM+cour];
          }
      }

      dx = ((x[1] - x[0]) + (x[n-1] - x[n-2])) / 2;
      tmp_x = x[0] - dx - x[n-1];

      for(cour=0;cour<nbsym;cour++) {
          ex->x_max[cour] = tmp_x + ex->x_max[ex->n_max+cour];
          ex->y_max[cour] = ex->y_max[ex->n_max+cour];
          ex->x_min[cour] = tmp_x + ex->x_min[ex->n_min+cour];
          ex->y_min[cour] = ex->y_min[ex->n_min+cour];
      }
      (ex->n_min) += nbsym;
      (ex->n_max) += nbsym;

      tmp_x = x[n-1] + dx - x[0];
      for(cour=0;cour<nbsym;cour++) {
          ex->x_max[ex->n_max+cour] = tmp_x + ex->x_max[nbsym+cour];
          ex->y_max[ex->n_max+cour] = ex->y_max[nbsym+cour];
          ex->x_min[ex->n_min+cour] = tmp_x + ex->x_min[nbsym+cour];
          ex->y_min[ex->n_min+cour] = ex->y_min[nbsym+cour];
      }
      (ex->n_min) += nbsym;
      (ex->n_max) += nbsym;

    } else {

      if (nbsym < NBSYM) {
          for(cour=0;cour<ex->n_max;cour++) {
              ex->x_max[nbsym+cour] = ex->x_max[NBSYM+cour];
              ex->y_max[nbsym+cour] = ex->y_max[NBSYM+cour];
          }
          for(cour=0;cour<ex->n_min;cour++) {
              ex->x_min[nbsym+cour] = ex->x_min[NBSYM+cour];
              ex->y_min[nbsym+cour] = ex->y_min[NBSYM+cour];
          }
      }
      
      /* select the symmetrized points and the axis of symmetry at the beginning of the signal*/
      if (ex->x_max[nbsym] < ex->x_min[nbsym]) { /* first = max */
          if (y[0] > ex->y_min[nbsym]) { /* the edge is not a min */
              if (2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                      ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
                      ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                      ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
                  }
              } else { /* symmetrized parts are long enough */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[cour] = 2*ex->x_max[nbsym]-ex->x_max[2*nbsym-cour];
                      ex->y_max[cour] = ex->y_max[2*nbsym-cour];
                      ex->x_min[cour] = 2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1-cour];
                      ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
                  }
              }
          } else { /* edge is a min -> sym with respect to the edge*/
              for(cour=0;cour<nbsym;cour++) {
                  ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                  ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
              }
              for(cour=0;cour<nbsym-1;cour++) {
                  ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-2-cour];
                  ex->y_min[cour] = ex->y_min[2*nbsym-2-cour];
              }
              ex->x_min[nbsym-1] = x[0];
              ex->y_min[nbsym-1] = y[0];
          }
      } else { /* first = min */
          if (y[0] < ex->y_max[nbsym]) { /* the edge is not a max */
              if (2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                      ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
                      ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                      ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
                  }
              } else { /* symmetrized parts are long enough */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[cour] = 2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1-cour];
                      ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
                      ex->x_min[cour] = 2*ex->x_min[nbsym]-ex->x_min[2*nbsym-cour];
                      ex->y_min[cour] = ex->y_min[2*nbsym-cour];
                  }
              }
          } else { /* edge is a max -> sym with respect to the edge*/
              for(cour=0;cour<nbsym;cour++) {
                  ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                  ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
              }
              for(cour=0;cour<nbsym-1;cour++) {
                  ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-2-cour];
                  ex->y_max[cour] = ex->y_max[2*nbsym-2-cour];
              }
              ex->x_max[nbsym-1] = x[0];
              ex->y_max[nbsym-1] = y[0];
          }
      }
      
      (ex->n_min) += nbsym-1;
      (ex->n_max) += nbsym-1;
      
      /* select the symmetrized points and the axis of symmetry at the end of the signal*/
      if (ex->x_max[ex->n_max] < ex->x_min[ex->n_min]) { /* last is a min */
          if (y[n-1] < ex->y_max[ex->n_max]) { /* the edge is not a max */
              if (2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                      ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
                      ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                      ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
                  }
              } else { /* symmetrized parts are long enough */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[ex->n_max+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-cour];
                      ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
                      ex->x_min[ex->n_min+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_min[ex->n_min-1-cour];
                      ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-1-cour];
                  }
              }
          } else { /* edge is a max -> sym with respect to the edge*/
              for(cour=0;cour<nbsym;cour++) {
                  ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                  ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
              }
              for(cour=0;cour<nbsym-1;cour++) {
                  ex->x_max[ex->n_max+2+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                  ex->y_max[ex->n_max+2+cour] = ex->y_max[ex->n_max-cour];
              }
              ex->x_max[ex->n_max+1] = x[n-1];
              ex->y_max[ex->n_max+1] = y[n-1];
          }
      } else {  /* last is a max */
          if (y[n-1] > ex->y_min[ex->n_min]) { /* the edge is not a min */
              if (2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                      ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
                      ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                      ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
                  }
              } else { /* symmetrized parts are long enough */
                  for(cour=0;cour<nbsym;cour++) {
                      ex->x_max[ex->n_max+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_max[ex->n_max-1-cour];
                      ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-1-cour];
                      ex->x_min[ex->n_min+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-cour];
                      ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
                  }
              }
          } else { /* edge is a min -> sym with respect to the edge*/
              for(cour=0;cour<nbsym;cour++) {
                  ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                  ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
              }
              for(cour=0;cour<nbsym-1;cour++) {
                  ex->x_min[ex->n_min+2+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                  ex->y_min[ex->n_min+2+cour] = ex->y_min[ex->n_min-cour];
              }
              ex->x_min[ex->n_min+1] = x[n-1];
              ex->y_min[ex->n_min+1] = y[n-1];
          }
      }
      
      (ex->n_min) = ex->n_min + nbsym + 1;
      (ex->n_max) = ex->n_max + nbsym + 1;
    }
}


/************************************************************************/
/*                                                                      */
/* FREE ALLOCATED MEMORY                                                */
/*                                                                      */
/************************************************************************/

void free_extr(extrema_t ex) {
    free(ex.x_max);
    free(ex.x_min);
    free(ex.y_max);
    free(ex.y_min);
    free(ex.i_max);
    free(ex.i_min);
    free(ex.x_tmp);
    free(ex.y_tmp);
}
