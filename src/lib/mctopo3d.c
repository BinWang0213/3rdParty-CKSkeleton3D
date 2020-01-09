/*
Copyright ESIEE (2009) 

m.couprie@esiee.fr

This software is an image processing library whose purpose is to be
used primarily for research and teaching.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software. You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
/* $Id: mctopo3d.c,v 1.12 2006/12/04 14:34:06 michel Exp $ */
/* 
Librairie mctopo3D : 

Calcul des nombres topologiques en 3D

Version calculant les nombres de connexité T et Tb directement
d'après la definition de G. Bertrand [Ber94].

[Ber94] G. Bertrand, "Simple points, topological numbers and geodesic
neighborhoods in cubic grids", Pattern Recognition Letters, 
Vol. 15, pp. 1003-1011, 1994.

Michel Couprie 1998-2006

Update nov. 2006 : modif geodesic_neighborhood pour compatibilité 64 bits
*/

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <mclifo.h>
#include <mcutil.h>
#include <mccodimage.h>
#include <mctopo3d.h>

/* globales privees */
static Lifo * LIFO_topo3d1 = NULL;
static Lifo * LIFO_topo3d2 = NULL;
static voxel cube_topo3d[27];
static voxel cubec_topo3d[27];
  
/* ========================================== */
void init_topo3d()
/* ========================================== */
#undef F_NAME
#define F_NAME "init_topo3d"
{
  LIFO_topo3d1 = CreeLifoVide(27);
  LIFO_topo3d2 = CreeLifoVide(27);
  if ((LIFO_topo3d1 == NULL) || (LIFO_topo3d2 == NULL))
  {   
    fprintf(stderr, "mccube() : CreeLifoVide failed\n");
    exit(0);
  }
  construitcube(cube_topo3d);
  construitcube(cubec_topo3d);
} /* init_topo3d() */

/* ========================================== */
void termine_topo3d()
/* ========================================== */
{
  LifoTermine(LIFO_topo3d1);
  LifoTermine(LIFO_topo3d2);
} /* termine_topo3d() */

/* ========================================== */
void construitcube(voxel * cube)
/* ========================================== */
{
  uint8_t x,y,z,u,v,w;
  pvoxel p;
  for (z = 0; z < 3; z++)
    for (y = 0; y < 3; y++)
      for (x = 0; x < 3; x++)
      {
      	p = &(cube[encode(x,y,z)]);
        p->x = x;
        p->y = y;
        p->z = z;
        p->n = encode(x,y,z);
      	p->val = 0;
      	p->lab = 0;
      	p->lab2 = 0;
      	p->n6v = 0;
      	p->n12v = 0;
      	p->n8v = 0;
      	p->n18v = 0;
      	p->n26v = 0;

        if ((x == 1) && (y == 1) && (z == 1))  p->type = centre;
        else
          if (abs(1-x)+abs(1-y)+abs(1-z) == 1) p->type = face;
          else
            if (abs(1-x)+abs(1-y)+abs(1-z) <= 2) p->type = arete;
            else
              p->type = coin;

      	for (w = 0; w < 3; w++)
       	  for (v = 0; v < 3; v++)
      	    for (u = 0; u < 3; u++)
      	    {
              if (abs(u-x)+abs(v-y)+abs(w-z) == 1)
	      {
                p->v6[p->n6v++] = &(cube[encode(u,v,w)]);
                p->v18[p->n18v++] = &(cube[encode(u,v,w)]);
                p->v26[p->n26v++] = &(cube[encode(u,v,w)]);
	      }
              else
              if (max(abs(u-x), max(abs(v-y), abs(w-z))) == 1)
              {
                if (abs(u-x)+abs(v-y)+abs(w-z) <= 2)
		{
                  p->v12[p->n12v++] = &(cube[encode(u,v,w)]);
                  p->v18[p->n18v++] = &(cube[encode(u,v,w)]);
                  p->v26[p->n26v++] = &(cube[encode(u,v,w)]);
		}
                else 
		{
                  p->v8[p->n8v++] = &(cube[encode(u,v,w)]);
                  p->v26[p->n26v++] = &(cube[encode(u,v,w)]);
		}
              }
      	    }  /* for w v u */
      } /* for z y x */
} /* construitcube() */

/* ========================================== */
uint32_t encodecube()
/* ========================================== */
{
  uint8_t n;
  pvoxel p;
  uint32_t i = 0;

  for (n = 0; n < 27; n++)
  {
    p = &(cube_topo3d[n]);
    if (p->val) i = i | (1<<n);        
  } /* for n */
  return i;
} /* encodecube() */

/* ========================================== */
void geodesic_neighborhood(voxel * cube, uint8_t connex, uint8_t s)
/* ========================================== */
#undef F_NAME
#define F_NAME ""
/* 
  met a 1 le champ lab des points appartenant au voisinage geodesique d'ordre s du point central,
  met a 0 les autres
*/
{
  uint8_t n;
  pvoxel p, pp, pc;
  Lifo * LIFOtmp;
  
  if ((LIFO_topo3d1 == NULL) || (LIFO_topo3d2 == NULL))
  { 
    fprintf(stderr, "geodesic_neighborhood: LIFO_topo3d1 and LIFO_topo3d2 must be allocated\n"); 
    exit(0); 
  }

  if (s < 1) 
  { 
    fprintf(stderr, "geodesic_neighborhood: order %d not allowed (must be > 0)\n", s); 
    exit(0); 
  }
  
  /* met a 0 le champ lab des points du cube */
  for (n = 0; n < 27; n++) cube[n].lab = 0;

  /* met a 1 le champ lab des voisins de valeur 1 du point central pc */
  pc = &(cube[13]);
  for (n = 0; n < pc->n6v; n++)
  {
    p = pc->v6[n];
    if (p->val == 1)
    { p->lab = 1; LifoPush(LIFO_topo3d1, (int32_t)(p-pc)); }
  }
  if (connex > 6)
    for (n = 0; n < pc->n12v; n++)
    {
      p = pc->v12[n];
      if (p->val == 1)
      { p->lab = 1; LifoPush(LIFO_topo3d1, (int32_t)(p-pc)); }
    }
  if (connex > 18)
    for (n = 0; n < pc->n8v; n++)
    {
      p = pc->v8[n];
      if (p->val == 1)
      { p->lab = 1; LifoPush(LIFO_topo3d1, (int32_t)(p-pc)); }
    }
  s--;

  while (s > 0)
  {
    while (!LifoVide(LIFO_topo3d1))
    {
      p = pc + LifoPop(LIFO_topo3d1);
      /* met a 1 le champ lab des voisins de valeur 1 du point p (sauf pc) */
      for (n = 0; n < p->n6v; n++)
      {
        pp = p->v6[n];
        if ((pp != pc) && (pp->val == 1) && (pp->lab == 0))
        { pp->lab = 1; LifoPush(LIFO_topo3d2, (int32_t)(pp-pc)); }
      }
      if (connex > 6)
        for (n = 0; n < p->n12v; n++)
        {
          pp = p->v12[n];
          if ((pp != pc) && (pp->val == 1) && (pp->lab == 0))
          { pp->lab = 1; LifoPush(LIFO_topo3d2, (int32_t)(pp-pc)); }
        }
      if (connex > 18)
        for (n = 0; n < p->n8v; n++)
        {
          pp = p->v8[n];
          if ((pp != pc) && (pp->val == 1) && (pp->lab == 0))
          { pp->lab = 1; LifoPush(LIFO_topo3d2, (int32_t)(pp-pc)); }
        }
    } /* while (!LifoVide(LIFO_topo3d1)) */
    s--;
    LIFOtmp = LIFO_topo3d1;
    LIFO_topo3d1 = LIFO_topo3d2;
    LIFO_topo3d2 = LIFOtmp;
  } /* while (s > 0) */

  LifoFlush(LIFO_topo3d1);
  
} /* geodesic_neighborhood() */

/* ========================================== */
void G6(voxel * cube)
/* ========================================== */
{
  geodesic_neighborhood(cube, 6, 2);	
} /* G6() */

/* ========================================== */
void G6p(voxel * cube)
/* ========================================== */
{
  geodesic_neighborhood(cube, 6, 3);	
} /* G6p() */

/* ========================================== */
void G18(voxel * cube)
/* ========================================== */
{
  geodesic_neighborhood(cube, 18, 2);	
} /* G18() */

/* ========================================== */
void G26(voxel * cube)
/* ========================================== */
{
  geodesic_neighborhood(cube, 26, 1);	
} /* G26() */

/* ========================================== */
uint8_t nbcomp(voxel * cube, uint8_t connex)
/* ========================================== */
/*
  retourne le nombre de composantes connexes de l'objet marque par un lab=1 
*/
{
  uint8_t ncc;
  uint8_t n,v;
  pvoxel p,pp;

  ncc = 0;
  for (n = 0; n < 27; n++) cube[n].lab2 = 0;
  for (n = 0; n < 27; n++)
  {
    p = &(cube[n]);
    if ((p->lab == 1) && (p->lab2 == 0))
    {
      ncc++;
      p->lab2 = ncc;
      LifoPush(LIFO_topo3d1, (int32_t)p);
      while (!LifoVide(LIFO_topo3d1))
      {
        p = (pvoxel)LifoPop(LIFO_topo3d1);
        for (v = 0; v < p->n6v; v++)
        {
          pp = p->v6[v];
          if ((pp->lab == 1) && (pp->lab2 == 0))
	  {
            pp->lab2 = ncc;
            LifoPush(LIFO_topo3d1, (int32_t)pp);
	  }
        } /* for v */
        if (connex > 6)
          for (v = 0; v < p->n12v; v++)
          {
            pp = p->v12[v];
            if ((pp->lab == 1) && (pp->lab2 == 0)) 
	    {
              pp->lab2 = ncc;
              LifoPush(LIFO_topo3d1, (int32_t)pp);
	    }
          } /* for v */
        if (connex > 18)
          for (v = 0; v < p->n8v; v++)
          {
            pp = p->v8[v];
            if ((pp->lab == 1) && (pp->lab2 == 0))
	    {
              pp->lab2 = ncc;
              LifoPush(LIFO_topo3d1, (int32_t)pp);
	    }
          } /* for v */
      } /* while (!LifoVide(LIFO_topo3d1)) */
    } /* if */
  } /* for n */

  return ncc;
} /* nbcomp() */

/* ========================================== */
uint8_t nbvois6(voxel * cube)
/* ========================================== */
/*
  retourne le nombre de 6-voisins du point central appartenant a l'objet
*/
{
  uint8_t v, nbvois;
  pvoxel p,pp;

  nbvois = 0;
  p = &(cube[13]);
  for (v = 0; v < p->n6v; v++)
  {
    pp = p->v6[v];
    if (pp->val) nbvois++;
  }

  return nbvois;
} /* nbvois6() */

/* ========================================== */
uint8_t nbvois26(voxel * cube)
/* ========================================== */
/*
  retourne le nombre de 26-voisins du point central appartenant a l'objet
*/
{
  uint8_t v, nbvois;
  pvoxel p,pp;

  nbvois = 0;
  p = &(cube[13]);
  for (v = 0; v < p->n26v; v++)
  {
    pp = p->v26[v];
    if (pp->val) nbvois++;
  }

  return nbvois;
} /* nbvois26() */

/* ========================================== */
int32_t nbvoisc6(
  uint8_t *B,            /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ========================================== */
/*
  retourne le nombre de 6-voisins du point central de niveau nul
*/
{
  int32_t nbvois = 0;
  if ((i%rs!=rs-1) && !B[i+1])    nbvois++;
  if (((i%ps)>=rs) && !B[i-rs])   nbvois++;
  if ((i%rs!=0) && !B[i-1])       nbvois++;
  if (((i%ps)<ps-rs) && !B[i+rs]) nbvois++;
  if ((i>=ps) && !B[i-ps])        nbvois++;
  if ((i<N-ps) && !B[i+ps])      nbvois++;
  return nbvois;
} /* nbvoisc6() */

/* ========================================== */
int32_t nbvoisc18(
  uint8_t *B,            /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ========================================== */
/*
  retourne le nombre de 18-voisins du point central de niveau nul
*/
{
  int32_t nbvois = 0;
  if (((i<N-ps)&&(i%rs!=rs-1)) && !B[ps+i+1]) nbvois++;
  if (((i<N-ps)&&(i%ps>=rs)) && !B[ps+i-rs]) nbvois++;
  if (((i<N-ps)&&(i%rs!=0)) && !B[ps+i-1]) nbvois++;
  if (((i<N-ps)&&(i%ps<ps-rs)) && !B[ps+i+rs]) nbvois++;
  if (((i<N-ps)) && !B[ps+i]) nbvois++;
  if (((i%rs!=rs-1)) && !B[i+1]) nbvois++;
  if (((i%rs!=rs-1)&&(i%ps>=rs)) && !B[i+1-rs]) nbvois++;
  if (((i%ps>=rs)) && !B[i-rs]) nbvois++;
  if (((i%ps>=rs)&&(i%rs!=0)) && !B[i-rs-1]) nbvois++;
  if (((i%rs!=0)) && !B[i-1]) nbvois++;
  if (((i%rs!=0)&&(i%ps<ps-rs)) && !B[i-1+rs]) nbvois++;
  if (((i%ps<ps-rs)) && !B[i+rs]) nbvois++;
  if (((i%ps<ps-rs)&&(i%rs!=rs-1)) && !B[i+rs+1]) nbvois++;
  if (((i>=ps)&&(i%rs!=rs-1)) && !B[-ps+i+1]) nbvois++;
  if (((i>=ps)&&(i%ps>=rs)) && !B[-ps+i-rs]) nbvois++;
  if (((i>=ps)&&(i%rs!=0)) && !B[-ps+i-1]) nbvois++;
  if (((i>=ps)&&(i%ps<ps-rs)) && !B[-ps+i+rs]) nbvois++;
  if (((i>=ps)) && !B[-ps+i]) nbvois++;
  return nbvois;
} /* nbvoisc18() */

/* ========================================== */
int32_t nbvoisc26(
  uint8_t *B,            /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ========================================== */
/*
  retourne le nombre de 26-voisins du point central de niveau nul
*/
{
  int32_t nbvois = 0;
  if (((i<N-ps)&&(i%rs!=rs-1)) && !B[ps+i+1]) nbvois++;
  if (((i<N-ps)&&(i%rs!=rs-1)&&(i%ps>=rs)) && !B[ps+i+1-rs]) nbvois++;
  if (((i<N-ps)&&(i%ps>=rs)) && !B[ps+i-rs]) nbvois++;
  if (((i<N-ps)&&(i%ps>=rs)&&(i%rs!=0)) && !B[ps+i-rs-1]) nbvois++;
  if (((i<N-ps)&&(i%rs!=0)) && !B[ps+i-1]) nbvois++;
  if (((i<N-ps)&&(i%rs!=0)&&(i%ps<ps-rs)) && !B[ps+i-1+rs]) nbvois++;
  if (((i<N-ps)&&(i%ps<ps-rs)) && !B[ps+i+rs]) nbvois++;
  if (((i<N-ps)&&(i%ps<ps-rs)&&(i%rs!=rs-1)) && !B[ps+i+rs+1]) nbvois++;
  if (((i<N-ps)) && !B[ps+i]) nbvois++;
  if (((i%rs!=rs-1)) && !B[i+1]) nbvois++;
  if (((i%rs!=rs-1)&&(i%ps>=rs)) && !B[i+1-rs]) nbvois++;
  if (((i%ps>=rs)) && !B[i-rs]) nbvois++;
  if (((i%ps>=rs)&&(i%rs!=0)) && !B[i-rs-1]) nbvois++;
  if (((i%rs!=0)) && !B[i-1]) nbvois++;
  if (((i%rs!=0)&&(i%ps<ps-rs)) && !B[i-1+rs]) nbvois++;
  if (((i%ps<ps-rs)) && !B[i+rs]) nbvois++;
  if (((i%ps<ps-rs)&&(i%rs!=rs-1)) && !B[i+rs+1]) nbvois++;
  if (((i>=ps)&&(i%rs!=rs-1)) && !B[-ps+i+1]) nbvois++;
  if (((i>=ps)&&(i%rs!=rs-1)&&(i%ps>=rs)) && !B[-ps+i+1-rs]) nbvois++;
  if (((i>=ps)&&(i%ps>=rs)) && !B[-ps+i-rs]) nbvois++;
  if (((i>=ps)&&(i%ps>=rs)&&(i%rs!=0)) && !B[-ps+i-rs-1]) nbvois++;
  if (((i>=ps)&&(i%rs!=0)) && !B[-ps+i-1]) nbvois++;
  if (((i>=ps)&&(i%rs!=0)&&(i%ps<ps-rs)) && !B[-ps+i-1+rs]) nbvois++;
  if (((i>=ps)&&(i%ps<ps-rs)) && !B[-ps+i+rs]) nbvois++;
  if (((i>=ps)&&(i%ps<ps-rs)&&(i%rs!=rs-1)) && !B[-ps+i+rs+1]) nbvois++;
  if (((i>=ps)) && !B[-ps+i]) nbvois++;
  return nbvois;
} /* nbvoisc26() */

/* ========================================== */
int32_t nbvoiso6(
  uint8_t *B,            /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ========================================== */
/*
  retourne le nombre de 6-voisins du point central de niveau NON nul
*/
{
  int32_t nbvois = 0;
  if ((i%rs!=rs-1) && B[i+1])    nbvois++;
  if (((i%ps)>=rs) && B[i-rs])   nbvois++;
  if ((i%rs!=0) && B[i-1])       nbvois++;
  if (((i%ps)<ps-rs) && B[i+rs]) nbvois++;
  if ((i>=ps) && B[i-ps])        nbvois++;
  if ((i<N-ps) && B[i+ps])      nbvois++;
  return nbvois;
} /* nbvoiso6() */

/* ========================================== */
int32_t nbvoiso18(
  uint8_t *B,            /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ========================================== */
/*
  retourne le nombre de 18-voisins du point central de niveau NON nul
*/
{
  int32_t nbvois = 0;
  if (((i<N-ps)&&(i%rs!=rs-1)) && B[ps+i+1]) nbvois++;
  if (((i<N-ps)&&(i%ps>=rs)) && B[ps+i-rs]) nbvois++;
  if (((i<N-ps)&&(i%rs!=0)) && B[ps+i-1]) nbvois++;
  if (((i<N-ps)&&(i%ps<ps-rs)) && B[ps+i+rs]) nbvois++;
  if (((i<N-ps)) && B[ps+i]) nbvois++;
  if (((i%rs!=rs-1)) && B[i+1]) nbvois++;
  if (((i%rs!=rs-1)&&(i%ps>=rs)) && B[i+1-rs]) nbvois++;
  if (((i%ps>=rs)) && B[i-rs]) nbvois++;
  if (((i%ps>=rs)&&(i%rs!=0)) && B[i-rs-1]) nbvois++;
  if (((i%rs!=0)) && B[i-1]) nbvois++;
  if (((i%rs!=0)&&(i%ps<ps-rs)) && B[i-1+rs]) nbvois++;
  if (((i%ps<ps-rs)) && B[i+rs]) nbvois++;
  if (((i%ps<ps-rs)&&(i%rs!=rs-1)) && B[i+rs+1]) nbvois++;
  if (((i>=ps)&&(i%rs!=rs-1)) && B[-ps+i+1]) nbvois++;
  if (((i>=ps)&&(i%ps>=rs)) && B[-ps+i-rs]) nbvois++;
  if (((i>=ps)&&(i%rs!=0)) && B[-ps+i-1]) nbvois++;
  if (((i>=ps)&&(i%ps<ps-rs)) && B[-ps+i+rs]) nbvois++;
  if (((i>=ps)) && B[-ps+i]) nbvois++;
  return nbvois;
} /* nbvoiso18() */

/* ========================================== */
int32_t nbvoiso26(
  uint8_t *B,            /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ========================================== */
/*
  retourne le nombre de 26-voisins du point central de niveau NON nul
*/
{
  int32_t nbvois = 0;
  if (((i<N-ps)&&(i%rs!=rs-1)) && B[ps+i+1]) nbvois++;
  if (((i<N-ps)&&(i%rs!=rs-1)&&(i%ps>=rs)) && B[ps+i+1-rs]) nbvois++;
  if (((i<N-ps)&&(i%ps>=rs)) && B[ps+i-rs]) nbvois++;
  if (((i<N-ps)&&(i%ps>=rs)&&(i%rs!=0)) && B[ps+i-rs-1]) nbvois++;
  if (((i<N-ps)&&(i%rs!=0)) && B[ps+i-1]) nbvois++;
  if (((i<N-ps)&&(i%rs!=0)&&(i%ps<ps-rs)) && B[ps+i-1+rs]) nbvois++;
  if (((i<N-ps)&&(i%ps<ps-rs)) && B[ps+i+rs]) nbvois++;
  if (((i<N-ps)&&(i%ps<ps-rs)&&(i%rs!=rs-1)) && B[ps+i+rs+1]) nbvois++;
  if (((i<N-ps)) && B[ps+i]) nbvois++;
  if (((i%rs!=rs-1)) && B[i+1]) nbvois++;
  if (((i%rs!=rs-1)&&(i%ps>=rs)) && B[i+1-rs]) nbvois++;
  if (((i%ps>=rs)) && B[i-rs]) nbvois++;
  if (((i%ps>=rs)&&(i%rs!=0)) && B[i-rs-1]) nbvois++;
  if (((i%rs!=0)) && B[i-1]) nbvois++;
  if (((i%rs!=0)&&(i%ps<ps-rs)) && B[i-1+rs]) nbvois++;
  if (((i%ps<ps-rs)) && B[i+rs]) nbvois++;
  if (((i%ps<ps-rs)&&(i%rs!=rs-1)) && B[i+rs+1]) nbvois++;
  if (((i>=ps)&&(i%rs!=rs-1)) && B[-ps+i+1]) nbvois++;
  if (((i>=ps)&&(i%rs!=rs-1)&&(i%ps>=rs)) && B[-ps+i+1-rs]) nbvois++;
  if (((i>=ps)&&(i%ps>=rs)) && B[-ps+i-rs]) nbvois++;
  if (((i>=ps)&&(i%ps>=rs)&&(i%rs!=0)) && B[-ps+i-rs-1]) nbvois++;
  if (((i>=ps)&&(i%rs!=0)) && B[-ps+i-1]) nbvois++;
  if (((i>=ps)&&(i%rs!=0)&&(i%ps<ps-rs)) && B[-ps+i-1+rs]) nbvois++;
  if (((i>=ps)&&(i%ps<ps-rs)) && B[-ps+i+rs]) nbvois++;
  if (((i>=ps)&&(i%ps<ps-rs)&&(i%rs!=rs-1)) && B[-ps+i+rs+1]) nbvois++;
  if (((i>=ps)) && B[-ps+i]) nbvois++;
  return nbvois;
} /* nbvoiso26() */

/* ========================================== */
uint8_t T6(voxel * cube)
/* ========================================== */
{
  G6(cube);
  return nbcomp(cube, 6);
} /* T6() */

/* ========================================== */
uint8_t T6p(voxel * cube)
/* ========================================== */
{
  G6p(cube);
  return nbcomp(cube, 6);
} /* T6p() */

/* ========================================== */
uint8_t T18(voxel * cube)
/* ========================================== */
{
  G18(cube);
  return nbcomp(cube, 18);
} /* T18() */

/* ========================================== */
uint8_t T26(voxel * cube)
/* ========================================== */
{
  G26(cube);
  return nbcomp(cube, 26);
} /* T26() */

/* ========================================== */
static uint8_t simple(voxel * cube, voxel * cubec, uint8_t connex)
/* ========================================== */
#undef F_NAME
#define F_NAME ""
{
  switch (connex)
  {
    case 6: return (uint8_t)((T6(cube) == 1) && (T26(cubec) == 1));
    case 18: return (uint8_t)((T18(cube) == 1) && (T6p(cubec) == 1));
    case 26: return (uint8_t)((T26(cube) == 1) && (T6(cubec) == 1));
    default: 
      fprintf(stderr, "simple: mauvaise connexite : %d\n", connex); 
      exit(0); 
  } /* switch (connex) */
} /* simple() */

/* ==================================== */
int32_t preparecubes(
  uint8_t *B,            /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*
  Transfere le voisinage de i pour l'image 3d img dans les 
  structures cube_topo3d (vois. original) et cubec_topo3d (complementaire).
  ATTENTION: i ne doit pas etre un point de bord (test a faire avant).
 */
{

  /* plan "HAUT" (+ps) */
  if (B[ps+i+1])    cube_topo3d[23].val = 1; else cube_topo3d[23].val = 0;
  if (B[ps+i+1-rs]) cube_topo3d[20].val = 1; else cube_topo3d[20].val = 0;
  if (B[ps+i-rs])   cube_topo3d[19].val = 1; else cube_topo3d[19].val = 0;
  if (B[ps+i-rs-1]) cube_topo3d[18].val = 1; else cube_topo3d[18].val = 0;
  if (B[ps+i-1])    cube_topo3d[21].val = 1; else cube_topo3d[21].val = 0;
  if (B[ps+i-1+rs]) cube_topo3d[24].val = 1; else cube_topo3d[24].val = 0;
  if (B[ps+i+rs])   cube_topo3d[25].val = 1; else cube_topo3d[25].val = 0;
  if (B[ps+i+rs+1]) cube_topo3d[26].val = 1; else cube_topo3d[26].val = 0;
  if (B[ps+i])      cube_topo3d[22].val = 1; else cube_topo3d[22].val = 0;
  /* plan "COURANT" () */
  if (B[i+1])       cube_topo3d[14].val = 1; else cube_topo3d[14].val = 0;
  if (B[i+1-rs])    cube_topo3d[11].val = 1; else cube_topo3d[11].val = 0;
  if (B[i-rs])      cube_topo3d[10].val = 1; else cube_topo3d[10].val = 0;
  if (B[i-rs-1])    cube_topo3d[9].val = 1; else cube_topo3d[9].val = 0;
  if (B[i-1])       cube_topo3d[12].val = 1; else cube_topo3d[12].val = 0;
  if (B[i-1+rs])    cube_topo3d[15].val = 1; else cube_topo3d[15].val = 0;
  if (B[i+rs])      cube_topo3d[16].val = 1; else cube_topo3d[16].val = 0;
  if (B[i+rs+1])    cube_topo3d[17].val = 1; else cube_topo3d[17].val = 0;
  if (B[i])         cube_topo3d[13].val = 1; else cube_topo3d[13].val = 0;
  /* plan "BAS" (-ps) */
  if (B[-ps+i+1])    cube_topo3d[5].val = 1; else cube_topo3d[5].val = 0;
  if (B[-ps+i+1-rs]) cube_topo3d[2].val = 1; else cube_topo3d[2].val = 0;
  if (B[-ps+i-rs])   cube_topo3d[1].val = 1; else cube_topo3d[1].val = 0;
  if (B[-ps+i-rs-1]) cube_topo3d[0].val = 1; else cube_topo3d[0].val = 0;
  if (B[-ps+i-1])    cube_topo3d[3].val = 1; else cube_topo3d[3].val = 0;
  if (B[-ps+i-1+rs]) cube_topo3d[6].val = 1; else cube_topo3d[6].val = 0;
  if (B[-ps+i+rs])   cube_topo3d[7].val = 1; else cube_topo3d[7].val = 0;
  if (B[-ps+i+rs+1]) cube_topo3d[8].val = 1; else cube_topo3d[8].val = 0;
  if (B[-ps+i])      cube_topo3d[4].val = 1; else cube_topo3d[4].val = 0;
  
  for (i = 0; i < 27; i++) 
    if (cube_topo3d[i].val == 1) cubec_topo3d[i].val = 0; else cubec_topo3d[i].val = 1;
} /* preparecubes() */

/* ==================================== */
int32_t preparecubesh(
  uint8_t *img,          /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t h,                       /* seuil */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*
  Transfere le voisinage de i pour l'image 3d img seuillee au niveau h dans les 
  structures cube_topo3d (vois. original) et cubec_topo3d (complementaire).
  ATTENTION: i ne doit pas etre un point de bord (test a faire avant).
 */
{
  /* plan "ARRIERE" (+ps) */
  if (img[ps+i+1]>=h)    cube_topo3d[17].val = 1; else cube_topo3d[17].val = 0;
  if (img[ps+i+1-rs]>=h) cube_topo3d[26].val = 1; else cube_topo3d[26].val = 0;
  if (img[ps+i-rs]>=h)   cube_topo3d[25].val = 1; else cube_topo3d[25].val = 0;
  if (img[ps+i-rs-1]>=h) cube_topo3d[24].val = 1; else cube_topo3d[24].val = 0;
  if (img[ps+i-1]>=h)    cube_topo3d[15].val = 1; else cube_topo3d[15].val = 0;
  if (img[ps+i-1+rs]>=h) cube_topo3d[6].val = 1; else cube_topo3d[6].val = 0;
  if (img[ps+i+rs]>=h)   cube_topo3d[7].val = 1; else cube_topo3d[7].val = 0;
  if (img[ps+i+rs+1]>=h) cube_topo3d[8].val = 1; else cube_topo3d[8].val = 0;
  if (img[ps+i]>=h)      cube_topo3d[16].val = 1; else cube_topo3d[16].val = 0;
  /* plan "COURANT" () */
  if (img[i+1]>=h)       cube_topo3d[14].val = 1; else cube_topo3d[14].val = 0;
  if (img[i+1-rs]>=h)    cube_topo3d[23].val = 1; else cube_topo3d[23].val = 0;
  if (img[i-rs]>=h)      cube_topo3d[22].val = 1; else cube_topo3d[22].val = 0;
  if (img[i-rs-1]>=h)    cube_topo3d[21].val = 1; else cube_topo3d[21].val = 0;
  if (img[i-1]>=h)       cube_topo3d[12].val = 1; else cube_topo3d[12].val = 0;
  if (img[i-1+rs]>=h)    cube_topo3d[3].val = 1; else cube_topo3d[3].val = 0;
  if (img[i+rs]>=h)      cube_topo3d[4].val = 1; else cube_topo3d[4].val = 0;
  if (img[i+rs+1]>=h)    cube_topo3d[5].val = 1; else cube_topo3d[5].val = 0;
  if (img[i]>=h)         cube_topo3d[13].val = 1; else cube_topo3d[13].val = 0;
  /* plan "AVANT" (-ps) */
  if (img[-ps+i+1]>=h)    cube_topo3d[11].val = 1; else cube_topo3d[11].val = 0;
  if (img[-ps+i+1-rs]>=h) cube_topo3d[20].val = 1; else cube_topo3d[20].val = 0;
  if (img[-ps+i-rs]>=h)   cube_topo3d[19].val = 1; else cube_topo3d[19].val = 0;
  if (img[-ps+i-rs-1]>=h) cube_topo3d[18].val = 1; else cube_topo3d[18].val = 0;
  if (img[-ps+i-1]>=h)    cube_topo3d[9].val = 1; else cube_topo3d[9].val = 0;
  if (img[-ps+i-1+rs]>=h) cube_topo3d[0].val = 1; else cube_topo3d[0].val = 0;
  if (img[-ps+i+rs]>=h)   cube_topo3d[1].val = 1; else cube_topo3d[1].val = 0;
  if (img[-ps+i+rs+1]>=h) cube_topo3d[2].val = 1; else cube_topo3d[2].val = 0;
  if (img[-ps+i]>=h)      cube_topo3d[10].val = 1; else cube_topo3d[10].val = 0;
  
  for (i = 0; i < 27; i++) 
    if (cube_topo3d[i].val == 1) cubec_topo3d[i].val = 0; else cubec_topo3d[i].val = 1;
} /* preparecubesh() */

/* ==================================== */
int32_t preparecubesh_l(
  uint32_t *img,          /* pointeur base image */
  int32_t i,                       /* index du point */
  int32_t h,                      /* seuil */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*
  Transfere le voisinage de i pour l'image 3d img seuillee au niveau h dans les 
  structures cube_topo3d (vois. original) et cubec_topo3d (complementaire).
  ATTENTION: i ne doit pas etre un point de bord (test a faire avant).
 */
{
  /* plan "ARRIERE" (+ps) */
  if (img[ps+i+1]>=h)    cube_topo3d[17].val = 1; else cube_topo3d[17].val = 0;
  if (img[ps+i+1-rs]>=h) cube_topo3d[26].val = 1; else cube_topo3d[26].val = 0;
  if (img[ps+i-rs]>=h)   cube_topo3d[25].val = 1; else cube_topo3d[25].val = 0;
  if (img[ps+i-rs-1]>=h) cube_topo3d[24].val = 1; else cube_topo3d[24].val = 0;
  if (img[ps+i-1]>=h)    cube_topo3d[15].val = 1; else cube_topo3d[15].val = 0;
  if (img[ps+i-1+rs]>=h) cube_topo3d[6].val = 1; else cube_topo3d[6].val = 0;
  if (img[ps+i+rs]>=h)   cube_topo3d[7].val = 1; else cube_topo3d[7].val = 0;
  if (img[ps+i+rs+1]>=h) cube_topo3d[8].val = 1; else cube_topo3d[8].val = 0;
  if (img[ps+i]>=h)      cube_topo3d[16].val = 1; else cube_topo3d[16].val = 0;
  /* plan "COURANT" () */
  if (img[i+1]>=h)       cube_topo3d[14].val = 1; else cube_topo3d[14].val = 0;
  if (img[i+1-rs]>=h)    cube_topo3d[23].val = 1; else cube_topo3d[23].val = 0;
  if (img[i-rs]>=h)      cube_topo3d[22].val = 1; else cube_topo3d[22].val = 0;
  if (img[i-rs-1]>=h)    cube_topo3d[21].val = 1; else cube_topo3d[21].val = 0;
  if (img[i-1]>=h)       cube_topo3d[12].val = 1; else cube_topo3d[12].val = 0;
  if (img[i-1+rs]>=h)    cube_topo3d[3].val = 1; else cube_topo3d[3].val = 0;
  if (img[i+rs]>=h)      cube_topo3d[4].val = 1; else cube_topo3d[4].val = 0;
  if (img[i+rs+1]>=h)    cube_topo3d[5].val = 1; else cube_topo3d[5].val = 0;
  if (img[i]>=h)         cube_topo3d[13].val = 1; else cube_topo3d[13].val = 0;
  /* plan "AVANT" (-ps) */
  if (img[-ps+i+1]>=h)    cube_topo3d[11].val = 1; else cube_topo3d[11].val = 0;
  if (img[-ps+i+1-rs]>=h) cube_topo3d[20].val = 1; else cube_topo3d[20].val = 0;
  if (img[-ps+i-rs]>=h)   cube_topo3d[19].val = 1; else cube_topo3d[19].val = 0;
  if (img[-ps+i-rs-1]>=h) cube_topo3d[18].val = 1; else cube_topo3d[18].val = 0;
  if (img[-ps+i-1]>=h)    cube_topo3d[9].val = 1; else cube_topo3d[9].val = 0;
  if (img[-ps+i-1+rs]>=h) cube_topo3d[0].val = 1; else cube_topo3d[0].val = 0;
  if (img[-ps+i+rs]>=h)   cube_topo3d[1].val = 1; else cube_topo3d[1].val = 0;
  if (img[-ps+i+rs+1]>=h) cube_topo3d[2].val = 1; else cube_topo3d[2].val = 0;
  if (img[-ps+i]>=h)      cube_topo3d[10].val = 1; else cube_topo3d[10].val = 0;
  
  for (i = 0; i < 27; i++) 
    if (cube_topo3d[i].val == 1) cubec_topo3d[i].val = 0; else cubec_topo3d[i].val = 1;
} /* preparecubesh_l() */

/* ******************************************************************************* */
/* ******************************************************************************* */
/*                               PRIMITIVES 3D BINAIRES                            */
/* ******************************************************************************* */
/* ******************************************************************************* */

/* ==================================== */
int32_t top6(                   /* pour un objet en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N,                       /* taille image */
  int32_t *t,
  int32_t *tb)                     /* resultats */
/* ==================================== */
/*
  ATTENTION: p ne doit pas etre un point de bord (test a faire avant).
*/
{
  preparecubes(img, p, rs, ps, N);
  *t = T6(cube_topo3d);
  *tb = T26(cubec_topo3d);
} /* top6() */

/* ==================================== */
int32_t top18(                   /* pour un objet en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N,                       /* taille image */
  int32_t *t,
  int32_t *tb)                     /* resultats */
/* ==================================== */
/*
  ATTENTION: p ne doit pas etre un point de bord (test a faire avant).
*/
{
  preparecubes(img, p, rs, ps, N);
  *t = T18(cube_topo3d);
  *tb = T6p(cubec_topo3d);
} /* top18() */

/* ==================================== */
int32_t top26(                   /* pour un objet en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N,                       /* taille image */
  int32_t *t,
  int32_t *tb)                     /* resultats */
/* ==================================== */
/*
  ATTENTION: p ne doit pas etre un point de bord (test a faire avant).
*/
{
  preparecubes(img, p, rs, ps, N);
  *t = T26(cube_topo3d);
  *tb = T6(cubec_topo3d);
} /* top26() */

/* ==================================== */
int32_t simple6(                   /* pour un objet en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
#undef F_NAME
#define F_NAME "simple6"
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubes(img, p, rs, ps, N);
  return ((T6(cube_topo3d) == 1) && (T26(cubec_topo3d) == 1));
} /* simple6() */

/* ==================================== */
int32_t simple18(                  /* pour un objet en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
#undef F_NAME
#define F_NAME "simple18"
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubes(img, p, rs, ps, N);
  return ((T18(cube_topo3d) == 1) && (T6p(cubec_topo3d) == 1));
} /* simple18() */

/* ==================================== */
int32_t simple26(                  /* pour un objet en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
#undef F_NAME
#define F_NAME "simple26"
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubes(img, p, rs, ps, N);
  return ((T26(cube_topo3d) == 1) && (T6(cubec_topo3d) == 1));
} /* simple26() */

/* ==================================== */
int32_t simple6h(                   /* pour un objet en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t h,                       /* seuil */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
#undef F_NAME
#define F_NAME "simple6h"
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, h, rs, ps, N);
  return ((T6(cube_topo3d) == 1) && (T26(cubec_topo3d) == 1));
} /* simple6h() */

/* ==================================== */
int32_t simple18h(                  /* pour un objet en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t h,                       /* seuil */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
#undef F_NAME
#define F_NAME "simple18h"
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, h, rs, ps, N);
  return ((T18(cube_topo3d) == 1) && (T6p(cubec_topo3d) == 1));
} /* simple18h() */

/* ==================================== */
int32_t simple26h(                  /* pour un objet en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t h,                       /* seuil */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
#undef F_NAME
#define F_NAME "simple26h"
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, h, rs, ps, N);
  return ((T26(cube_topo3d) == 1) && (T6(cubec_topo3d) == 1));
} /* simple26h() */

/* ==================================== */
int32_t tbar6h(               /* pour un objet en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t h,                       /* seuil */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return -1;
  preparecubesh(img, p, h, rs, ps, N);
  return T26(cubec_topo3d);
} /* tbar6h() */

/* ==================================== */
int32_t tbar26h(              /* pour un objet en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t h,                       /* seuil */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return -1;
  preparecubesh(img, p, h, rs, ps, N);
  return T6(cubec_topo3d);
} /* tbar26h() */

/* ========================================== */
uint8_t P_simple(voxel * cube, voxel * cubep, voxel * cubec, uint8_t connex)
/* ========================================== */
#undef F_NAME
#define F_NAME "P_simple"
/*
  cube contient X
  cubep contient P
  cubec (auxiliaire) n'a pas besoin d'etre initialise
  d'apres: "Some topological properties of surfaces in Z3", G. Bertrand & R. Malgouyres
           Theoreme 6
*/
{
  uint8_t n;
  uint8_t v;
  pvoxel x;  /* point central de cube */
  pvoxel y;  /* point de cube */
  pvoxel xc; /* point central de cubec */
  pvoxel yc; /* point de cubec */
  pvoxel xp; /* point central de cubep */
  pvoxel yp; /* point de cubep */

  for (n = 0; n < 27; n++) if (cube[n].val == 1) cubec[n].val = 0; else cubec[n].val = 1;

  switch (connex) /* teste la condition 2 (theoreme 6) */
  {
    case 6:  
      if (T26(cubec) != 1) return 0; 
      break;
    case 18: 
      if (T6p(cubec) != 1) return 0; 
      break;
    case 26: 
      if (T6(cubec) != 1) return 0; 
      break;
    default: 
      fprintf(stderr, "P_simple: mauvaise connexite : %d\n", connex); 
      exit(0); 
  } /* switch (connex) */
  
  x = &(cube[13]);
  xc = &(cubec[13]);
  xp = &(cubep[13]);
  switch (connex) /* teste la condition 4 (theoreme 6) */
  {
    case 6: 
      for (n = 0; n < x->n26v; n++)
      {
        yp = xp->v26[n];
        if (yp->val)
        {
          yc = xc->v26[n];
          v = yc->val;
          yc->val = 1;
          if (T26(cubec) != 1) return 0;
          yc->val = v;
        } /* if (yp->val) */
      } /* for (n = 0; n < x->n26v; n++) */
      break;
    case 18: 
      for (n = 0; n < x->n6v; n++)
      {
        yp = xp->v6[n];
        if (yp->val)
        {
          yc = xc->v6[n];
          v = yc->val;
          yc->val = 1;
          if (T6p(cubec) != 1) return 0;
          yc->val = v;
        } /* if (yp->val) */
      } /* for (n = 0; n < x->n6v; n++) */
      break;
    case 26: 
      for (n = 0; n < x->n6v; n++)
      {
        yp = xp->v6[n];
        if (yp->val)
        {
          yc = xc->v6[n];
          v = yc->val;
          yc->val = 1;
          if (T6(cubec) != 1) return 0;
          yc->val = v;
        } /* if (yp->val) */
      } /* for (n = 0; n < x->n6v; n++) */
      break;
    default: 
      fprintf(stderr, "P_simple: mauvaise connexite : %d\n", connex); 
      exit(0); 
  } /* switch (connex) */
  
  for (n = 0; n < 27; n++) /* calcule et range dans cubec l'ensemble R = X - P  */
  {
    y = &(cube[n]);
    yp = &(cubep[n]);
    yc = &(cubec[n]);
    if (y->val && !yp->val) yc->val = 1; else yc->val = 0;
  } /* for (n = 0; n < 27; n++) */

  switch (connex) /* teste la condition 1 (theoreme 6) */
  {
    case 6:  
      if (T6(cubec) != 1) return 0;
      break;
    case 18: 
      if (T18(cubec) != 1) return 0;
      break;
    case 26: 
      if (T26(cubec) != 1) return 0;
      break;
    default: 
      fprintf(stderr, "P_simple: mauvaise connexite : %d\n", connex); 
      exit(0); 
  } /* switch (connex) */
  
  switch (connex) /* teste la condition 3 (theoreme 6) */
  {
    case 6: 
      for (n = 0; n < x->n6v; n++)
      {
        yp = xp->v6[n];
        if (yp->val)
        {
          yc = xc->v6[n];
          v = yc->val;
          yc->val = 1;
          if (T6(cubec) != 1) return 0;
          yc->val = v;
        } /* if (yp->val) */
      } /* for (n = 0; n < x->n6v; n++) */
      break;
    case 18: 
      for (n = 0; n < x->n18v; n++)
      {
        yp = xp->v18[n];
        if (yp->val)
        {
          yc = xc->v18[n];
          v = yc->val;
          yc->val = 1;
          if (T18(cubec) != 1) return 0;
          yc->val = v;
        } /* if (yp->val) */
      } /* for (n = 0; n < x->n18v; n++) */
      break;
    case 26: 
      for (n = 0; n < x->n26v; n++)
      {
        yp = xp->v26[n];
        if (yp->val)
        {
          yc = xc->v26[n];
          v = yc->val;
          yc->val = 1;
          if (T26(cubec) != 1) return 0;
          yc->val = v;
        } /* if (yp->val) */
      } /* for (n = 0; n < x->n26v; n++) */
      break;
    default: 
      fprintf(stderr, "P_simple: mauvaise connexite : %d\n", connex); 
      exit(0); 
  } /* switch (connex) */
  return 1;
} /* P_simple() */





/* ******************************************************************************* */
/* ******************************************************************************* */
/*                                  PRIMITIVES 3D NDG                              */
/* ******************************************************************************* */
/* ******************************************************************************* */

/* ==================================== */
int32_t pdestr6(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return ((T26(cube_topo3d) == 1) && (T6(cubec_topo3d) == 1));
} /* pdestr6() */

/* ==================================== */
int32_t pdestr18(                  /* pour des minima en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return ((T6p(cube_topo3d) == 1) && (T18(cubec_topo3d) == 1));
} /* pdestr18() */

/* ==================================== */
int32_t pdestr26(                  /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return ((T6(cube_topo3d) == 1) && (T26(cubec_topo3d) == 1));
} /* pdestr26() */

/* ==================================== */
int32_t plevdestr6(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return (T6(cubec_topo3d) == 1);
} /* plevdestr6() */

/* ==================================== */
int32_t plevdestr18(                  /* pour des minima en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return (T18(cubec_topo3d) == 1);
} /* plevdestr18() */

/* ==================================== */
int32_t plevdestr26(                  /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return (T26(cubec_topo3d) == 1);
} /* plevdestr26() */

/* ==================================== */
int32_t pconstr6(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return ((T26(cube_topo3d) == 1) && (T6(cubec_topo3d) == 1));
} /* pconstr6() */

/* ==================================== */
int32_t pconstr18(                  /* pour des minima en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return ((T6p(cube_topo3d) == 1) && (T18(cubec_topo3d) == 1));
} /* pconstr18() */

/* ==================================== */
int32_t pconstr26(                  /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return ((T6(cube_topo3d) == 1) && (T26(cubec_topo3d) == 1));
} /* pconstr26() */

/* ==================================== */
int32_t plevconstr6(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return (T26(cube_topo3d) == 1);
} /* plevconstr6() */

/* ==================================== */
int32_t plevconstr18(                  /* pour des minima en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return (T6p(cube_topo3d) == 1);
} /* plevconstr18() */

/* ==================================== */
int32_t plevconstr26(                  /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return (T6(cube_topo3d) == 1);
} /* plevconstr26() */

/* ==================================== */
int32_t peak6(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return (T26(cube_topo3d) == 0);
} /* peak6() */

/* ==================================== */
int32_t peak26(                    /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return (T6(cube_topo3d) == 0);
} /* peak26() */

/* ==================================== */
int32_t well6(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return (T6(cubec_topo3d) == 0);
} /* well6() */

/* ==================================== */
int32_t well26(                    /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return (T26(cubec_topo3d) == 0);
} /* well26() */

/* ==================================== */
uint8_t alpha26m(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* retourne le sup des valeurs < img[x] dans le 26-voisinage de x, */
/* ou img[x] si pas de telles valeurs */
/* ==================================== */
{
	register uint8_t val = *(img+p);
	register int32_t q;
	register uint8_t v;
	register int32_t alpha = NDG_MIN - 1;
        register int32_t k;

        for (k = 0; k < 26; k += 1)
        {
          q = voisin26(p, k, rs, ps, N);
          if ((q != -1) && ((v=img[q]) < val) && ((int32_t)v > alpha)) alpha = (int32_t)v;
	}
        if (alpha == NDG_MIN - 1) 
          return val;
        else
          return (uint8_t)alpha;
} /* alpha26m() */

/* ==================================== */
uint32_t alpha26m_l(
  uint32_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* retourne le sup des valeurs < img[x] dans le 26-voisinage de x, */
/* ou img[x] si pas de telles valeurs */
/* ==================================== */
{
	register uint32_t val = *(img+p);
	register int32_t q;
	register uint32_t v;
	register int32_t alpha = NDG_MIN - 1;
        register int32_t k;

        for (k = 0; k < 26; k += 1)
        {
          q = voisin26(p, k, rs, ps, N);
          if ((q != -1) && ((v=img[q]) < val) && ((int32_t)v > alpha)) alpha = (int32_t)v;
	}
        if (alpha == NDG_MIN - 1) 
          return val;
        else
          return (uint32_t)alpha;
} /* alpha26m_l() */

/* ==================================== */
uint8_t alpha6m(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* retourne le sup des valeurs < img[x] dans le 6-voisinage de x, */
/* ou img[x] si pas de telles valeurs */
/* ==================================== */
{
	register uint8_t val = *(img+p);
	register int32_t q;
	register uint8_t v;
	register int32_t alpha = NDG_MIN - 1;
        register int32_t k;

        for (k = 0; k <= 10; k += 2)
        {
          q = voisin6(p, k, rs, ps, N);
          if ((q != -1) && ((v=img[q]) < val) && ((int32_t)v > alpha)) alpha = (int32_t)v;
	}
        if (alpha == NDG_MIN - 1) 
          return val;
        else
          return (uint8_t)alpha;
} /* alpha6m() */

/* ==================================== */
uint8_t alpha26p(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* retourne le inf des valeurs > img[x] dans le 26-voisinage de x, */
/* ou img[x] si pas de telles valeurs */
/* ==================================== */
{
	register uint8_t val = *(img+p);
	register int32_t q;
	register uint8_t v;
	register int32_t alpha = NDG_MAX + 1;
        register int32_t k;

        for (k = 0; k < 26; k += 1)
        {
          q = voisin26(p, k, rs, ps, N);
          if ((q != -1) && ((v=img[q]) > val) && ((int32_t)v < alpha)) alpha = (int32_t)v;
	}
        if (alpha == NDG_MAX + 1) 
          return val;
        else
          return (uint8_t)alpha;
} /* alpha26p() */

/* ==================================== */
uint8_t alpha6p(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* retourne le inf des valeurs > img[x] dans le 6-voisinage de x, */
/* ou img[x] si pas de telles valeurs */
/* ==================================== */
{
	register uint8_t val = *(img+p);
	register int32_t q;
	register uint8_t v;
	register int32_t alpha = NDG_MAX + 1;
        register int32_t k;

        for (k = 0; k <= 10; k += 2)
        {
          q = voisin6(p, k, rs, ps, N);
          if ((q != -1) && ((v=img[q]) > val) && ((int32_t)v < alpha)) alpha = (int32_t)v;
	}
        if (alpha == NDG_MAX + 1) 
          return val;
        else
          return (uint8_t)alpha;
} /* alpha6p() */

/* ==================================== */
uint8_t delta6m( 
/* retourne la valeur max. a laquelle p est destructible - minima 6-connexes */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{	
  uint8_t ret, sav = img[p];
  while (pdestr6(img, p, rs, ps, N)) img[p] = alpha26m(img, p, rs, ps, N);
  ret = img[p];
  img[p] = sav;
  return ret;
} /* delta6m() */

/* ==================================== */
uint8_t delta26m( 
/* retourne la valeur max. a laquelle p est destructible - minima 26-connexes */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{	
  uint8_t ret, sav = img[p];
  while (pdestr26(img, p, rs, ps, N)) img[p] = alpha26m(img, p, rs, ps, N);
  ret = img[p];
  img[p] = sav;
  return ret;
} /* delta26m() */

/* ==================================== */
uint8_t delta6p( 
/* retourne la valeur min. a laquelle p est constructible - minima 6-connexes */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{	
  uint8_t ret, sav = img[p];
  while (pconstr6(img, p, rs, ps, N)) img[p] = alpha26p(img, p, rs, ps, N);
  ret = img[p];
  img[p] = sav;
  return ret;
} /* delta6p() */

/* ==================================== */
uint8_t delta26p( 
/* retourne la valeur min. a laquelle p est constructible - minima 26-connexes */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{	
  uint8_t ret, sav = img[p];
  while (pconstr26(img, p, rs, ps, N)) img[p] = alpha26p(img, p, rs, ps, N);
  ret = img[p];
  img[p] = sav;
  return ret;
} /* delta26p() */

/* ==================================== */
int32_t separant6(  /* teste si un point est separant - minima 6-connexes
	         ie- s'il est separant pour une coupe <= img[p] */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  int32_t k, q;

  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;

  preparecubesh(img, p, img[p], rs, ps, N);
  if (T6(cubec_topo3d) >= 2) return 1;

  for (k = 0; k < 26; k += 1)
  {
    q = voisin26(p, k, rs, ps, N);
    if ((q != -1) && (img[q] <= img[p]))
    {
      preparecubesh(img, p, img[q], rs, ps, N);
      if (T6(cubec_topo3d) >= 2) return 1;
    }
  }	
  return 0;
} /* separant6() */

/* ==================================== */
int32_t hseparant6(  /* teste si un point est hseparant - minima 6-connexes
	         ie- s'il est separant pour la coupe h */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t h,                       /* parametre */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  int32_t k, q;

  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;

  preparecubesh(img, p, h, rs, ps, N);
  if (T6(cubec_topo3d) >= 2) return 1;
  return 0;
} /* hseparant6() */

/* ==================================== */
int32_t hfseparant6(  /* teste si un point est hfseparant - minima 6-connexes
	         ie- s'il est separant pour une coupe c telle que h < c <= img[p] */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t h,                       /* parametre */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  int32_t k, q;

  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;

  preparecubesh(img, p, img[p], rs, ps, N);
  if (T6(cubec_topo3d) >= 2) return 1;

  for (k = 0; k < 26; k += 1)
  {
    q = voisin26(p, k, rs, ps, N);
    if ((q != -1) && (img[q] > h) && (img[q] <= img[p]))
    {
      preparecubesh(img, p, img[q], rs, ps, N);
      if (T6(cubec_topo3d) >= 2) return 1;
    }
  }	
  return 0;
} /* hfseparant6() */

/* ==================================== */
int32_t filsombre6(                /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*
   pour la coupe K (>img[p]), le point doit etre
   - soit un point isole (T == 1; Tb == 0)
   - soit une extremite de courbe (T == 1; Tb == 1 et card(K) = 1)
   - soit un point de courbe (T == 1; Tb == 2 et card(K) = 2)
   - soit un point de croisement de courbe (T == 1; Tb == n et card(K) = n)

*/
{
  int32_t T, Tb, Nb;
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  T = T26(cube_topo3d);
  if (T != 1) return 0;
  Tb = T6(cubec_topo3d);
  if (Tb == 0) return 1;
  Nb = nbvois6(cubec_topo3d);
  if (Tb > 0) return (Nb == Tb);
  return 0;
} /* filsombre6() */

/* ==================================== */
int32_t filsombre26(               /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*
   pour la coupe K (>img[p]), le point doit etre
   - soit un point isole (T == 1; Tb == 0)
   - soit une extremite de courbe (T == 1; Tb == 1 et card(K) = 1)
   - soit un point de courbe (T == 1; Tb == 2 et card(K) = 2)
   - soit un point de croisement de courbe (T == 1; Tb == n et card(K) = n)

*/
{
  int32_t T, Tb, Nb;
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  T = T6(cube_topo3d);
  if (T != 1) return 0;
  Tb = T26(cubec_topo3d);
  if (Tb == 0) return 1;
  Nb = nbvois26(cubec_topo3d);
  if (Tb > 0) return (Nb == Tb);
  return 0;
} /* filsombre26() */

/* ==================================== */
int32_t filclair6(                /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*
   pour la coupe K (>=img[p]), le point doit etre
   - soit un point isole (Tb == 1; T == 0)
   - soit une extremite de courbe (Tb == 1; T == 1 et card(K) = 1)
   - soit un point de courbe (Tb == 1; T == 2 et card(K) = 2)
   - soit un point de croisement de courbe (Tb == 1; T == n et card(K) = n)
*/
{
  int32_t T, Tb, Nb;
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  Tb = T6(cubec_topo3d);
  if (Tb != 1) return 0;
  T = T26(cube_topo3d);
  if (T == 0) return 1;
  Nb = nbvois26(cube_topo3d);
  if (T > 0) return (Nb == T);
  return 0;
} /* filclair6() */

/* ==================================== */
int32_t filclair26(                /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*
   pour la coupe K (>=img[p]), le point doit etre
   - soit un point isole (Tb == 1; T == 0)
   - soit une extremite de courbe (Tb == 1; T == 1 et card(K) = 1)
   - soit un point de courbe (Tb == 1; T == 2 et card(K) = 2)
   - soit un point de croisement de courbe (Tb == 1; T == n et card(K) = n)
*/
{
  int32_t T, Tb, Nb;
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  Tb = T26(cubec_topo3d);
  if (Tb != 1) return 0;
  T = T6(cube_topo3d);
  if (T == 0) return 1;
  Nb = nbvois6(cube_topo3d);
  if (T > 0) return (Nb == T);
  return 0;
} /* filclair26() */

/* ==================================== */
int32_t t6mm(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return T6(cubec_topo3d);
} /* t6mm() */

/* ==================================== */
int32_t t6m(                   /* pour des minima en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return T6(cubec_topo3d);
} /* t6m() */

/* ==================================== */
int32_t t26mm(                   /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return T26(cubec_topo3d);
} /* t26mm() */

/* ==================================== */
int32_t t26m(                   /* pour des minima en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return T26(cubec_topo3d);
} /* t26m() */

/* ==================================== */
int32_t t6pp(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return T6(cube_topo3d);
} /* t6pp() */

/* ==================================== */
int32_t t6p(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return T6(cube_topo3d);
} /* t6p() */

/* ==================================== */
int32_t t26pp(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p]+1, rs, ps, N);
  return T26(cube_topo3d);
} /* t26pp() */

/* ==================================== */
int32_t t26p(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh(img, p, img[p], rs, ps, N);
  return T26(cube_topo3d);
} /* t26p() */

/* ==================================== */
int32_t t26pp_l(
  uint32_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh_l(img, p, img[p]+1, rs, ps, N);
  return T26(cube_topo3d);
} /* t26pp_l() */

/* ==================================== */
int32_t t6pp_l(
  uint32_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    return 0;
  preparecubesh_l(img, p, img[p]+1, rs, ps, N);
  return T6(cube_topo3d);
} /* t6pp_l() */

/* ==================================== */
void nbtopoh3d26_l( /* pour les minima en 26-connexite */ 
  uint32_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  uint32_t h,
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N,                       /* taille image */
  int32_t *t6p,
  int32_t *t26mm)
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    {
      printf("ERREUR: nbtopoh3d26_l: point de bord\n");
      exit(0);
    }
  preparecubesh_l(img, p, h, rs, ps, N);
  *t6p = T6(cube_topo3d);
  *t26mm = T26(cubec_topo3d);
} /* nbtopoh3d26_l() */

/* ==================================== */
void nbtopoh3d6_l( /* pour les minima en 6-connexite */ 
  uint32_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  uint32_t h,
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N,                       /* taille image */
  int32_t *t26p,
  int32_t *t6mm)
/* ==================================== */
{
  if ((p < ps) || (p >= N-ps) ||         /* premier ou dernier plan */
      (p%ps < rs) || (p%ps >= ps-rs) ||  /* premiere ou derniere colonne */
      (p%rs == 0) || (p%rs == rs-1))     /* premiere ou derniere ligne */
    {
      printf("ERREUR: nbtopoh3d6_l: point de bord\n");
      exit(0);
    }
  preparecubesh_l(img, p, h, rs, ps, N);
  *t26p = T26(cube_topo3d);
  *t6mm = T6(cubec_topo3d);
} /* nbtopoh3d6_l() */

/* ==================================== */
int32_t bordext6(uint8_t *F, int32_t x, int32_t rs, int32_t ps, int32_t N)
/* ==================================== */
/* teste si x a un 6-voisin a 0 */
{
  int32_t k, y;
  for (k = 0; k <= 10; k += 2) /* parcourt les voisins en 6-connexite */
  {
    y = voisin6(x, k, rs, ps, N);
    if ((y != -1) && (F[y] == 0)) return 1;
  } /* for k */      
  return 0;
} /* bordext6() */

/* ==================================== */
int32_t bordext26(uint8_t *F, int32_t x, int32_t rs, int32_t ps, int32_t N)
/* ==================================== */
/* teste si x a un 26-voisin a 0 */
{
  int32_t k, y;
  for (k = 0; k < 26; k += 1) /* parcourt les voisins en 26-connexite */
  {
    y = voisin26(x, k, rs, ps, N);
    if ((y != -1) && (F[y] == 0)) return 1;
  } /* for k */      
  return 0;
} /* bordext26() */

/* ==================================== */
int32_t curve6( /* point de courbe en 6-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*  ATTENTION: i ne doit pas etre un point de bord (test a faire avant). */
{
  if (img[p] == 0) return 0;
  preparecubes(img, p, rs, ps, N);
  if ((T6(cube_topo3d) == 2) && (nbvoiso6(img, p, rs, ps, N) == 2)) return 1;
  return 0;
} /* curve6() */

/* ==================================== */
int32_t curve18( /* point de courbe en 18-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*  ATTENTION: i ne doit pas etre un point de bord (test a faire avant). */
{
  if (img[p] == 0) return 0;
  preparecubes(img, p, rs, ps, N);
  if ((T18(cube_topo3d) == 2) && (nbvoiso18(img, p, rs, ps, N) == 2)) return 1;
  return 0;
} /* curve18() */

/* ==================================== */
int32_t curve26( /* point de courbe en 26-connexite */
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */
/* ==================================== */
/*  ATTENTION: i ne doit pas etre un point de bord (test a faire avant). */
{
  if (img[p] == 0) return 0;
  preparecubes(img, p, rs, ps, N);
  if ((T26(cube_topo3d) == 2) && (nbvoiso26(img, p, rs, ps, N) == 2)) return 1;
  return 0;
} /* curve26() */
