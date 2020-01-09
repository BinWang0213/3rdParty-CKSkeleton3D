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
/* $Id: lskelpar3d.c,v 1.6 2006/12/04 14:34:06 michel Exp $ */
/* 
   3D Parallel thinning algorithms based on critical kernels

   Michel Couprie, 2006
*/

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <mccodimage.h>
#include <mcimage.h>
#include <mctopo.h>
#include <mctopo3d.h>
#include <mcutil.h>
#include <lskelpar3d.h>

#define S_OBJECT      1
#define S_SIMPLE      2
#define S_2M_CRUCIAL  4
#define S_1M_CRUCIAL  8
#define S_0M_CRUCIAL 16
#define S_CURVE      32
#define S_SURF       64
#define S_SELECTED  128

#define IS_OBJECT(f)     (f&S_OBJECT)
#define IS_SIMPLE(f)     (f&S_SIMPLE)
#define IS_2M_CRUCIAL(f) (f&S_2M_CRUCIAL)
#define IS_1M_CRUCIAL(f) (f&S_1M_CRUCIAL)
#define IS_0M_CRUCIAL(f) (f&S_0M_CRUCIAL)
#define IS_CURVE(f)      (f&S_CURVE)
#define IS_SURF(f)       (f&S_SURF)
#define IS_SELECTED(f)   (f&S_SELECTED)

#define SET_OBJECT(f)     (f|=S_OBJECT)
#define SET_SIMPLE(f)     (f|=S_SIMPLE)
#define SET_2M_CRUCIAL(f) (f|=S_2M_CRUCIAL)
#define SET_1M_CRUCIAL(f) (f|=S_1M_CRUCIAL)
#define SET_0M_CRUCIAL(f) (f|=S_0M_CRUCIAL)
#define SET_CURVE(f)      (f|=S_CURVE)
#define SET_SURF(f)       (f|=S_SURF)
#define SET_SELECTED(f)   (f|=S_SELECTED)

#define UNSET_OBJECT(f)     (f&=~S_OBJECT)
#define UNSET_SIMPLE(f)     (f&=~S_SIMPLE)
#define UNSET_2M_CRUCIAL(f) (f&=~S_2M_CRUCIAL)
#define UNSET_1M_CRUCIAL(f) (f&=~S_1M_CRUCIAL)
#define UNSET_0M_CRUCIAL(f) (f&=~S_0M_CRUCIAL)
#define UNSET_SELECTED(f)   (f&=~S_SELECTED)

//#define VERBOSE
//#define DEBUG
#ifdef DEBUG
int32_t trace = 1;
#endif

/* ==================================== */
static void extract_vois(
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N,                       /* taille image */
  uint8_t *vois)    
/* 
  retourne dans "vois" les valeurs des 27 voisins de p, dans l'ordre suivant: 

               12      11      10       
               13       8       9
               14      15      16

		3	2	1			
		4      26	0
		5	6	7

               21      20      19
               22      17      18
               23      24      25

  le point p ne doit pas être un point de bord de l'image
*/
/* ==================================== */
{
#undef F_NAME
#define F_NAME "extract_vois"
  register uint8_t * ptr = img+p;
  if ((p%rs==rs-1) || (p%ps<rs) || (p%rs==0) || (p%ps>=ps-rs) || 
      (p < ps) || (p >= N-ps)) /* point de bord */
  {
    printf("%s: ERREUR: point de bord\n", F_NAME);
    exit(0);
  }
  vois[ 0] = *(ptr+1);
  vois[ 1] = *(ptr+1-rs);
  vois[ 2] = *(ptr-rs);
  vois[ 3] = *(ptr-rs-1);
  vois[ 4] = *(ptr-1);
  vois[ 5] = *(ptr-1+rs);
  vois[ 6] = *(ptr+rs);
  vois[ 7] = *(ptr+rs+1);

  vois[ 8] = *(ptr-ps);
  vois[ 9] = *(ptr-ps+1);
  vois[10] = *(ptr-ps+1-rs);
  vois[11] = *(ptr-ps-rs);
  vois[12] = *(ptr-ps-rs-1);
  vois[13] = *(ptr-ps-1);
  vois[14] = *(ptr-ps-1+rs);
  vois[15] = *(ptr-ps+rs);
  vois[16] = *(ptr-ps+rs+1);

  vois[17] = *(ptr+ps);
  vois[18] = *(ptr+ps+1);
  vois[19] = *(ptr+ps+1-rs);
  vois[20] = *(ptr+ps-rs);
  vois[21] = *(ptr+ps-rs-1);
  vois[22] = *(ptr+ps-1);
  vois[23] = *(ptr+ps-1+rs);
  vois[24] = *(ptr+ps+rs);
  vois[25] = *(ptr+ps+rs+1);

  vois[26] = *(ptr);
} /* extract_vois() */

/* ==================================== */
static void insert_vois(
  uint8_t *vois,			
  uint8_t *img,          /* pointeur base image */
  int32_t p,                       /* index du point */
  int32_t rs,                      /* taille rangee */
  int32_t ps,                      /* taille plan */
  int32_t N)                       /* taille image */    
/* 
  recopie vois dans le voisinage de p
  le point p ne doit pas être un point de bord de l'image
*/
/* ==================================== */
{
#undef F_NAME
#define F_NAME "insert_vois"
  register uint8_t * ptr = img+p;
  if ((p%rs==rs-1) || (p%ps<rs) || (p%rs==0) || (p%ps>=ps-rs) || 
      (p < ps) || (p >= N-ps)) /* point de bord */
  {
    printf("%s: ERREUR: point de bord\n", F_NAME);
    exit(0);
  }
  *(ptr+1) = vois[ 0];
  *(ptr+1-rs) = vois[ 1];
  *(ptr-rs) = vois[ 2];
  *(ptr-rs-1) = vois[ 3];
  *(ptr-1) = vois[ 4];
  *(ptr-1+rs) = vois[ 5];
  *(ptr+rs) = vois[ 6];
  *(ptr+rs+1) = vois[ 7];

  *(ptr-ps) = vois[ 8];
  *(ptr-ps+1) = vois[ 9];
  *(ptr-ps+1-rs) = vois[10];
  *(ptr-ps-rs) = vois[11];
  *(ptr-ps-rs-1) = vois[12];
  *(ptr-ps-1) = vois[13];
  *(ptr-ps-1+rs) = vois[14];
  *(ptr-ps+rs) = vois[15];
  *(ptr-ps+rs+1) = vois[16];

  *(ptr+ps) = vois[17];
  *(ptr+ps+1) = vois[18];
  *(ptr+ps+1-rs) = vois[19];
  *(ptr+ps-rs) = vois[20];
  *(ptr+ps-rs-1) = vois[21];
  *(ptr+ps-1) = vois[22];
  *(ptr+ps-1+rs) = vois[23];
  *(ptr+ps+rs) = vois[24];
  *(ptr+ps+rs+1) = vois[25];

  *(ptr) = vois[26];
} /* insert_vois() */

/* ==================================== */
static void print_vois(uint8_t *vois)    
/* 
   affiche vois (debug)
*/
/* ==================================== */
{
  printf("%2d %2d %2d     %2d %2d %2d     %2d %2d %2d\n", 
	 vois[12],vois[11],vois[10],vois[3],vois[2],vois[1],vois[21],vois[20],vois[19]);
  printf("%2d %2d %2d     %2d %2d %2d     %2d %2d %2d\n", 
	 vois[13],vois[8],vois[9],vois[4],vois[26],vois[0],vois[22],vois[17],vois[18]);
  printf("%2d %2d %2d     %2d %2d %2d     %2d %2d %2d\n\n", 
	 vois[14],vois[15],vois[16],vois[5],vois[6],vois[7],vois[23],vois[24],vois[25]);
} /* print_vois() */

/* ==================================== */
static void isometrieXZ_vois(uint8_t *vois) 
// effectue une isométrie du voisinage "vois" par échange des axes X et Z (+ symétries)
// cette isométrie est de plus une involution
/* ==================================== */
{
  uint8_t v[26];
  int32_t i;
  v[ 0] = vois[17];  v[ 1] = vois[20];  v[ 2] = vois[ 2];  v[ 3] = vois[11];
  v[ 4] = vois[ 8];  v[ 5] = vois[15];  v[ 6] = vois[ 6];  v[ 7] = vois[24];
  v[ 8] = vois[ 4];  v[ 9] = vois[22];  v[10] = vois[21];  v[11] = vois[ 3];
  v[12] = vois[12];  v[13] = vois[13];  v[14] = vois[14];  v[15] = vois[ 5];
  v[16] = vois[23];  v[17] = vois[ 0];  v[18] = vois[18];  v[19] = vois[19];
  v[20] = vois[ 1];  v[21] = vois[10];  v[22] = vois[ 9];  v[23] = vois[16];
  v[24] = vois[ 7];  v[25] = vois[25];
  for (i = 0; i < 26; i++) vois[i] = v[i];
} /* isometrieXZ_vois() */

/* ==================================== */
static void isometrieYZ_vois(uint8_t *vois)
// effectue une isométrie du voisinage "vois" par échange des axes Y et Z (+ symétries)  
// cette isométrie est de plus une involution
/* ==================================== */
{
  uint8_t v[26];
  int32_t i;
  v[ 0] = vois[ 0];  v[ 1] = vois[18];  v[ 2] = vois[17];  v[ 3] = vois[22];
  v[ 4] = vois[ 4];  v[ 5] = vois[13];  v[ 6] = vois[ 8];  v[ 7] = vois[ 9];
  v[ 8] = vois[ 6];  v[ 9] = vois[ 7];  v[10] = vois[25];  v[11] = vois[24];
  v[12] = vois[23];  v[13] = vois[ 5];  v[14] = vois[14];  v[15] = vois[15];
  v[16] = vois[16];  v[17] = vois[ 2];  v[18] = vois[ 1];  v[19] = vois[19];
  v[20] = vois[20];  v[21] = vois[21];  v[22] = vois[ 3];  v[23] = vois[12];
  v[24] = vois[11];  v[25] = vois[10];
  for (i = 0; i < 26; i++) vois[i] = v[i];
} /* isometrieYZ_vois() */

/* ==================================== */
int32_t match_end(uint8_t *v)
/* ==================================== */
/*
               12      11      10
               13       8       9
               14      15      16

		3	2	1
		4      26	0
		5	6	7

               21      20      19
               22      17      18
               23      24      25

Teste si au moins un des points 12, 11, 13, 8, 3, 2, 4 est objet et tous les autres fond
(aussi avec les isométries). 
*/
{
  uint8_t t;
#ifdef DEBUG
  if (trace)
  {  
    printf("match_end\n");
    print_vois(v);
  }
#endif

  if ((v[19] || v[20] || v[18] || v[17] || v[2] || v[1] || v[0]) &&
      !v[3 ] && !v[4 ] && !v[21] && !v[22] && 
      !v[25] && !v[24] && !v[7 ] && !v[6 ] && !v[23] && !v[5] && 
      !v[12] && !v[11] && !v[10] && !v[13] && 
      !v[8 ] && !v[9 ] && !v[14] && !v[15] && !v[16]) return 1;

  if ((v[21] || v[20] || v[22] || v[17] || v[2] || v[3] || v[4]) &&
      !v[1 ] && !v[0 ] && !v[19] && !v[18] && 
      !v[23] && !v[24] && !v[5 ] && !v[6 ] && !v[25] && !v[7] && 
      !v[12] && !v[11] && !v[10] && !v[13] && 
      !v[8 ] && !v[9 ] && !v[14] && !v[15] && !v[16]) return 1;

  if ((v[17] || v[22] || v[24] || v[23] || v[4] || v[5] || v[6]) &&
      !v[2 ] && !v[3 ] && !v[20] && !v[21] && 
      !v[18] && !v[25] && !v[0 ] && !v[7 ] && !v[19] && !v[1] && 
      !v[12] && !v[11] && !v[10] && !v[13] && 
      !v[8 ] && !v[9 ] && !v[14] && !v[15] && !v[16]) return 1;

  if ((v[17] || v[18] || v[24] || v[25] || v[0] || v[7] || v[6]) &&
      !v[2 ] && !v[1 ] && !v[20] && !v[19] && 
      !v[22] && !v[23] && !v[4 ] && !v[5 ] && !v[21] && !v[3] && 
      !v[12] && !v[11] && !v[10] && !v[13] && 
      !v[8 ] && !v[9 ] && !v[14] && !v[15] && !v[16]) return 1;

  if ((v[9 ] || v[8] || v[16] || v[15] || v[0] || v[7] || v[6]) &&
      !v[10] && !v[11] && !v[1 ] && !v[2 ] && 
      !v[13] && !v[14] && !v[4 ] && !v[5 ] && !v[12] && !v[3] && 
      !v[21] && !v[20] && !v[19] && !v[22] && 
      !v[17] && !v[18] && !v[23] && !v[24] && !v[25]) return 1;

  if ((v[13] || v[8] || v[14] || v[15] || v[4] || v[5] || v[6]) &&
      !v[12] && !v[11] && !v[3 ] && !v[2 ] && 
      !v[9 ] && !v[16] && !v[0 ] && !v[7 ] && !v[10] && !v[1] && 
      !v[21] && !v[20] && !v[19] && !v[22] && 
      !v[17] && !v[18] && !v[23] && !v[24] && !v[25]) return 1;

  if ((v[11] || v[10] || v[8] || v[9] || v[2] || v[1] || v[0]) &&
      !v[12] && !v[13] && !v[3 ] && !v[4 ] && 
      !v[15] && !v[16] && !v[6 ] && !v[7 ] && !v[14] && !v[5] && 
      !v[21] && !v[20] && !v[19] && !v[22] && 
      !v[17] && !v[18] && !v[23] && !v[24] && !v[25]) return 1;

  if ((v[12] || v[11] || v[13] || v[8] || v[3 ] || v[2 ] || v[4 ]) &&
      !v[10] && !v[9 ] && !v[1 ] && !v[0 ] && 
      !v[14] && !v[15] && !v[5 ] && !v[6 ] && !v[16] && !v[7] && 
      !v[21] && !v[20] && !v[19] && !v[22] && 
      !v[17] && !v[18] && !v[23] && !v[24] && !v[25]) return 1;

  return 0;
} // match_end()

/* ==================================== */
int32_t match_vois2(uint8_t *v)
/* ==================================== */
/*
               12      11      10       
               13       8       9
               14      15      16

		3	2	1			
		4      26	0
		5	6	7
Teste si les conditions suivantes sont réunies:
1: v[8] et v[26] doivent être dans l'objet et simples
2: for i = 0 to 7 do w[i] = v[i] || v[i+9] ; w[0...7] doit être non 2D-simple
Si le test réussit, les points 8, 26 sont marqués 2M_CRUCIAL
*/
{
  uint8_t t;
#ifdef DEBUG
  if (trace)
  {  
    printf("match_vois2\n");
    print_vois(v);
  }
#endif
  if (!IS_SIMPLE(v[8]) || !IS_SIMPLE(v[26])) return 0;
  if (v[0] || v[9]) t = 1; else t = 0;
  if (v[1] || v[10]) t |= 2;
  if (v[2] || v[11]) t |= 4;
  if (v[3] || v[12]) t |= 8;
  if (v[4] || v[13]) t |= 16;
  if (v[5] || v[14]) t |= 32;
  if (v[6] || v[15]) t |= 64;
  if (v[7] || v[16]) t |= 128;
  if ((t4b(t) == 1) && (t8(t) == 1)) return 0; // simple 2D
  SET_2M_CRUCIAL(v[8]);
  SET_2M_CRUCIAL(v[26]);
#ifdef DEBUG
  if (trace)
    printf("match !\n");
#endif
  return 1;
} // match_vois2()

/* ==================================== */
int32_t match_vois2s(uint8_t *v)
/* ==================================== */
/*
               12      11      10       
               13       8       9
               14      15      16

		3	2	1			
		4      26	0
		5	6	7
Pour les conditions de courbe et de surface.
Teste si les deux conditions suivantes sont réunies:
1: v[8] et v[26] doivent être simples
2: for i = 0 to 7 do w[i] = v[i] || v[i+9] ; w[0...7] doit être non 2D-simple
Si le test réussit, alors les points 8, 26 sont marqués 2M_CRUCIAL, de plus:
  Si t4b(w[0...7]) == 0 alors les points 8, 26 sont marqués SURF
  Sinon, si t8(w[0...7]) > 1 alors les points 8, 26 sont marqués CURVE
*/
{
  uint8_t t;
#ifdef DEBUG
  if (trace)
  {  
    printf("match_vois2s\n");
    print_vois(v);
  }
#endif
  if (!IS_SIMPLE(v[8]) || !IS_SIMPLE(v[26])) return 0;
  if (v[0] || v[9]) t = 1; else t = 0;
  if (v[1] || v[10]) t |= 2;
  if (v[2] || v[11]) t |= 4;
  if (v[3] || v[12]) t |= 8;
  if (v[4] || v[13]) t |= 16;
  if (v[5] || v[14]) t |= 32;
  if (v[6] || v[15]) t |= 64;
  if (v[7] || v[16]) t |= 128;
  if ((t4b(t) == 1) && (t8(t) == 1)) return 0; // simple 2D
  SET_2M_CRUCIAL(v[8]);
  SET_2M_CRUCIAL(v[26]);
  if (t4b(t) == 0) { SET_SURF(v[8]); SET_SURF(v[26]); }
  else if (t8(t) > 1) 
  { 
    SET_CURVE(v[8]); 
    SET_CURVE(v[26]); 
  }
#ifdef DEBUG
  if (trace)
    printf("match !\n");
#endif
  return 1;
} // match_vois2s()

/* ==================================== */
int32_t match_vois1(uint8_t *v)
/* ==================================== */
// A A  P1 P2  B B
// A A  P3 P4  B B
// avec pour localisations possibles :
// 12 11   3  2   21 20 
// 13  8   4 26   22 17
// et :
// 11 10    2 1   20 19
//  8  9   26 0   17 18
//
// Teste si les trois conditions suivantes sont réunies:
// 1: (P1 et P4) ou (P2 et P3)
// 2: tous les points Pi non nuls doivent être simples et non marqués 2M_CRUCIAL
// 3: A et B sont tous nuls ou [au moins un A non nul et au moins un B non nul]
// Si le test réussit, les points Pi non nuls sont marques 1M_CRUCIAL
{
  int32_t ret = 0;
#ifdef DEBUG
  if (trace)
  {  
    printf("match_vois1\n");
    print_vois(v);
  }
#endif
  if (!((v[2] && v[4]) || (v[3] && v[26]))) goto next1;
  if ((IS_OBJECT(v[2])  && (!IS_SIMPLE(v[2])  || IS_2M_CRUCIAL(v[2]))) ||
      (IS_OBJECT(v[3])  && (!IS_SIMPLE(v[3])  || IS_2M_CRUCIAL(v[3]))) ||
      (IS_OBJECT(v[4])  && (!IS_SIMPLE(v[4])  || IS_2M_CRUCIAL(v[4]))) ||
      (IS_OBJECT(v[26]) && (!IS_SIMPLE(v[26]) || IS_2M_CRUCIAL(v[26])))) goto next1;
  if ((v[12] || v[11] || v[13] || v[8] || v[21] || v[20] || v[22] || v[17]) &&
      ((!v[12] && !v[11] && !v[13] && !v[8]) || 
       (!v[21] && !v[20] && !v[22] && !v[17]))) goto next1;
  if (v[2])  SET_1M_CRUCIAL(v[2]);
  if (v[3])  SET_1M_CRUCIAL(v[3]);
  if (v[4])  SET_1M_CRUCIAL(v[4]);
  if (v[26]) SET_1M_CRUCIAL(v[26]);
  ret = 1;
 next1:
  if (!((v[2] && v[0]) || (v[1] && v[26]))) goto next2;
  if ((IS_OBJECT(v[2])  && (!IS_SIMPLE(v[2])  || IS_2M_CRUCIAL(v[2]))) ||
      (IS_OBJECT(v[1])  && (!IS_SIMPLE(v[1])  || IS_2M_CRUCIAL(v[1]))) ||
      (IS_OBJECT(v[0])  && (!IS_SIMPLE(v[0])  || IS_2M_CRUCIAL(v[0]))) ||
      (IS_OBJECT(v[26]) && (!IS_SIMPLE(v[26]) || IS_2M_CRUCIAL(v[26])))) goto next2;
  if ((v[10] || v[11] || v[9] || v[8] || v[19] || v[20] || v[18] || v[17]) &&
      ((!v[10] && !v[11] && !v[9] && !v[8]) || 
       (!v[19] && !v[20] && !v[18] && !v[17]))) goto next2;
  if (v[2])  SET_1M_CRUCIAL(v[2]);
  if (v[1])  SET_1M_CRUCIAL(v[1]);
  if (v[0])  SET_1M_CRUCIAL(v[0]);
  if (v[26]) SET_1M_CRUCIAL(v[26]);
  ret = 1;
 next2:
#ifdef DEBUG
  if (trace && ret)
    printf("match !\n");
#endif
  return ret;
} // match_vois1()

/* ==================================== */
int32_t match_vois1s(uint8_t *v)
/* ==================================== */
// A A  P1 P2  B B
// A A  P3 P4  B B
// avec pour localisations possibles :
// 12 11   3  2   21 20 
// 13  8   4 26   22 17
// et :
// 11 10    2 1   20 19
//  8  9   26 0   17 18
//
// Pour la condition de courbe. 
// Teste si les trois conditions suivantes sont réunies:
// 1: (P1 et P4) ou (P2 et P3)
// 2: tous les points Pi non nuls doivent être simples et non 2M_CRUCIAL
// 3: au moins un A non nul et au moins un B non nul
// Si le test réussit, les points Pi non nuls sont marques CURVE
{
  int32_t ret = 0;
#ifdef DEBUG
  if (trace)
  {  
    printf("match_vois1\n");
    print_vois(v);
  }
#endif
  if (!((v[2] && v[4]) || (v[3] && v[26]))) goto next1;
  if ((IS_OBJECT(v[2])  && (!IS_SIMPLE(v[2]) || IS_2M_CRUCIAL(v[2]))) ||
      (IS_OBJECT(v[3])  && (!IS_SIMPLE(v[3]) || IS_2M_CRUCIAL(v[3]))) ||
      (IS_OBJECT(v[4])  && (!IS_SIMPLE(v[4]) || IS_2M_CRUCIAL(v[4]))) ||
      (IS_OBJECT(v[26]) && (!IS_SIMPLE(v[26]) || IS_2M_CRUCIAL(v[26])))) goto next1;
  if ((!v[12] && !v[11] && !v[13] && !v[8]) || 
      (!v[21] && !v[20] && !v[22] && !v[17])) goto next1;
  if (v[2])  SET_CURVE(v[2]);
  if (v[3])  SET_CURVE(v[3]);
  if (v[4])  SET_CURVE(v[4]);
  if (v[26]) SET_CURVE(v[26]);
  ret = 1;
 next1:
  if (!((v[2] && v[0]) || (v[1] && v[26]))) goto next2;
  if ((IS_OBJECT(v[2])  && (!IS_SIMPLE(v[2]) || IS_2M_CRUCIAL(v[2]))) ||
      (IS_OBJECT(v[1])  && (!IS_SIMPLE(v[1]) || IS_2M_CRUCIAL(v[1]))) ||
      (IS_OBJECT(v[0])  && (!IS_SIMPLE(v[0]) || IS_2M_CRUCIAL(v[0]))) ||
      (IS_OBJECT(v[26]) && (!IS_SIMPLE(v[26]) || IS_2M_CRUCIAL(v[26])))) goto next2;
  if ((!v[10] && !v[11] && !v[9] && !v[8]) || 
      (!v[19] && !v[20] && !v[18] && !v[17])) goto next2;
  if (v[2])  SET_CURVE(v[2]);
  if (v[1])  SET_CURVE(v[1]);
  if (v[0])  SET_CURVE(v[0]);
  if (v[26]) SET_CURVE(v[26]);
  ret = 1;
 next2:
#ifdef DEBUG
  if (trace)
    printf("match !\n");
#endif
  return ret;
} // match_vois1s()

/* ==================================== */
int32_t match_vois0(uint8_t *v)
/* ==================================== */
/*
               12      11
               13       8

		3	2
		4      26

Teste si les conditions suivantes sont réunies:
1: au moins un des ensembles {12,26}, {11,4}, {13,2}, {8,3} est inclus dans l'objet, et
2: les points non nuls sont tous simples, non marqués 2M_CRUCIAL et non marqués 1M_CRUCIAL
Si le test réussit, les points non nuls sont marqués 0M_CRUCIAL
*/
{
  uint8_t t;
#ifdef DEBUG
  if (trace)
  {  
    printf("match_vois0\n");
    print_vois(v);
  }
#endif
  if (!((v[12]&&v[26]) || (v[11]&&v[4]) || (v[13]&&v[2]) || (v[8]&&v[3]) )) return 0;

  if (v[12] && (!IS_SIMPLE(v[12]) || IS_2M_CRUCIAL(v[12]) || IS_1M_CRUCIAL(v[12]))) return 0;
  if (v[26] && (!IS_SIMPLE(v[26]) || IS_2M_CRUCIAL(v[26]) || IS_1M_CRUCIAL(v[26]))) return 0;
  if (v[11] && (!IS_SIMPLE(v[11]) || IS_2M_CRUCIAL(v[11]) || IS_1M_CRUCIAL(v[11]))) return 0;
  if (v[ 4] && (!IS_SIMPLE(v[ 4]) || IS_2M_CRUCIAL(v[ 4]) || IS_1M_CRUCIAL(v[ 4]))) return 0;
  if (v[13] && (!IS_SIMPLE(v[13]) || IS_2M_CRUCIAL(v[13]) || IS_1M_CRUCIAL(v[13]))) return 0;
  if (v[ 2] && (!IS_SIMPLE(v[ 2]) || IS_2M_CRUCIAL(v[ 2]) || IS_1M_CRUCIAL(v[ 2]))) return 0;
  if (v[ 8] && (!IS_SIMPLE(v[ 8]) || IS_2M_CRUCIAL(v[ 8]) || IS_1M_CRUCIAL(v[ 8]))) return 0;
  if (v[ 3] && (!IS_SIMPLE(v[ 3]) || IS_2M_CRUCIAL(v[ 3]) || IS_1M_CRUCIAL(v[ 3]))) return 0;

  if (v[12]) SET_0M_CRUCIAL(v[12]);
  if (v[26]) SET_0M_CRUCIAL(v[26]);
  if (v[11]) SET_0M_CRUCIAL(v[11]);
  if (v[ 4]) SET_0M_CRUCIAL(v[ 4]);
  if (v[13]) SET_0M_CRUCIAL(v[13]);
  if (v[ 2]) SET_0M_CRUCIAL(v[ 2]);
  if (v[ 8]) SET_0M_CRUCIAL(v[ 8]);
  if (v[ 3]) SET_0M_CRUCIAL(v[ 3]);
#ifdef DEBUG
  if (trace)
    printf("match !\n");
#endif
  return 1;
} // match_vois0()

/* ==================================== */
int32_t match2(uint8_t *v)
/* ==================================== */
{
  int32_t i, ret = 0;
  if (match_vois2(v)) ret = 1;
  isometrieXZ_vois(v);
  if (match_vois2(v)) ret = 1;
  isometrieXZ_vois(v);
  isometrieYZ_vois(v);
  if (match_vois2(v)) ret = 1;
  isometrieYZ_vois(v);
  return ret;
} /* match2() */

/* ==================================== */
int32_t match2s(uint8_t *v)
/* ==================================== */
{
  int32_t i, ret = 0;
  if (match_vois2s(v)) ret = 1;
  isometrieXZ_vois(v);
  if (match_vois2s(v)) ret = 1;
  isometrieXZ_vois(v);
  isometrieYZ_vois(v);
  if (match_vois2s(v)) ret = 1;
  isometrieYZ_vois(v);
  return ret;
} /* match2s() */

/* ==================================== */
int32_t match1(uint8_t *v)
/* ==================================== */
{
  int32_t i, ret = 0;
  if (match_vois1(v)) ret = 1;
  isometrieXZ_vois(v);
  if (match_vois1(v)) ret = 1;
  isometrieXZ_vois(v);
  isometrieYZ_vois(v);
  if (match_vois1(v)) ret = 1;
  isometrieYZ_vois(v);
  return ret;
} /* match1() */

/* ==================================== */
int32_t match1s(uint8_t *v)
/* ==================================== */
{
  int32_t i, ret = 0;
  if (match_vois1s(v)) ret = 1;
  isometrieXZ_vois(v);
  if (match_vois1s(v)) ret = 1;
  isometrieXZ_vois(v);
  isometrieYZ_vois(v);
  if (match_vois1s(v)) ret = 1;
  isometrieYZ_vois(v);
  return ret;
} /* match1s() */

/* ==================================== */
int32_t match0(uint8_t *v)
/* ==================================== */
{
  int32_t i, ret = 0;
  if (match_vois0(v)) ret = 1;
  return ret;
} /* match0() */

/* ==================================== */
int32_t lskelMK3a(struct xvimage *image, 
	     int32_t nsteps,
	     struct xvimage *inhibit)
/* ==================================== */
/*
Squelette symétrique ultime
Algo MK3 données: S
Répéter jusqu'à stabilité
  P := voxels simples pour S
  R := voxels de P à préserver (match2, match1)
  T :=  [S  \  P]  \cup  R
  S := T \cup [S \ (T \oplus \Gamma_26*)]
*/
#undef F_NAME
#define F_NAME "lskelMK3a"
{ 
  int32_t i, j, x, y, z;
  int32_t rs = rowsize(image);     /* taille ligne */
  int32_t cs = colsize(image);     /* taille colonne */
  int32_t ds = depth(image);       /* nb plans */
  int32_t ps = rs * cs;            /* taille plan */
  int32_t N = ps * ds;             /* taille image */
  uint8_t *S = UCHARDATA(image);      /* l'image de depart */
  struct xvimage *t = copyimage(image); 
  uint8_t *T = UCHARDATA(t);
  struct xvimage *r = copyimage(image); 
  uint8_t *R = UCHARDATA(r);
  int32_t step, nonstab;
  uint8_t v[27];

  if (inhibit != NULL)
  {
    fprintf(stderr, "%s: inhibit image not implemented\n", F_NAME);
    return 0;
  }

  if (nsteps == -1) nsteps = 1000000000;

  for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;

  init_topo3d();

  /* ================================================ */
  /*               DEBUT ALGO                         */
  /* ================================================ */

  step = 0;
  nonstab = 1;
  while (nonstab && (step < nsteps))
  {
    nonstab = 0;
    step++;
#ifdef VERBOSE
    printf("step %d\n", step);
#endif

    // PREMIERE SOUS-ITERATION : MARQUE LES POINTS SIMPLES
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]) && simple26(S, i, rs, ps, N))
	SET_SIMPLE(S[i]);
#ifdef DEBUG
writeimage(image,"_S");
#endif
    // DEUXIEME SOUS-ITERATION : MARQUE LES CLIQUES CRUCIALES CORRESPONDANT AUX 2-FACES
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
#ifdef DEBUG
memset(R, 0, N);
for (i = 0; i < N; i++) if (IS_2M_CRUCIAL(S[i])) R[i] = 255;
writeimage(r,"_M2");
#endif

    // TROISIEME SOUS-ITERATION : MARQUE LES CLIQUES CRUCIALES CORRESPONDANT AUX 1-FACES
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
#ifdef DEBUG
printf("%d %d %d\n", i % rs, (i % ps) / rs, i / ps);
#endif
	extract_vois(S, i, rs, ps, N, v);
	if (match1(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
#ifdef DEBUG
memset(R, 0, N);
for (i = 0; i < N; i++) if (IS_1M_CRUCIAL(S[i])) R[i] = 255;
writeimage(r,"_M1");
#endif

    memset(T, 0, N);
    for (i = 0; i < N; i++) // T := [S \ P] \cup  R, où R représente les pts marqués
      if ((S[i] && !IS_SIMPLE(S[i])) || IS_2M_CRUCIAL(S[i]) || IS_1M_CRUCIAL(S[i]))
	T[i] = 1;
#ifdef DEBUG
writeimage(t,"_T");
#endif

    memset(R, 0, N);
    for (i = 0; i < N; i++)
      if (nbvoiso26(T, i, rs, ps, N) >= 1) R[i] = 1; // calcule R = Dilat(T)
    for (i = 0; i < N; i++)
      if (T[i] || (S[i] && !R[i])) T[i] = 1; else T[i] = 0; // T := T \cup [S \ R]

    for (i = 0; i < N; i++)
      if (S[i] && !T[i]) 
      {
	S[i] = 0; 
	nonstab = 1; 
      }
    for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;
  }

#ifdef VERBOSE1
    printf("number of steps: %d\n", step);
#endif

  for (i = 0; i < N; i++) if (S[i]) S[i] = 255; // normalize values

  freeimage(t);
  freeimage(r);
  termine_topo3d();
  return(1);
} /* lskelMK3a() */

/* ==================================== */
int32_t lskelEK3(struct xvimage *image, 
	     int32_t nsteps,
	     struct xvimage *inhibit)
/* ==================================== */
/*
Squelette symétrique curviligne basé sur les extrémités
Algo EK3 données: S
Répéter jusqu'à stabilité
  E := points extrémité de S
  P := voxels simples pour S et pas dans E
  C2 := voxels 2M-cruciaux (match2)
  C1 := voxels 1M-cruciaux (match1)
  C0 := voxels 0M-cruciaux (match0)
  P := P  \  [C2 \cup C1 \cup C0]
  S := S \ P
*/
#undef F_NAME
#define F_NAME "lskelEK3"
{ 
  int32_t i, j, x, y, z;
  int32_t rs = rowsize(image);     /* taille ligne */
  int32_t cs = colsize(image);     /* taille colonne */
  int32_t ds = depth(image);       /* nb plans */
  int32_t ps = rs * cs;            /* taille plan */
  int32_t N = ps * ds;             /* taille image */
  uint8_t *S = UCHARDATA(image);      /* l'image de depart */
  struct xvimage *t = copyimage(image); 
  uint8_t *T = UCHARDATA(t);
  uint8_t *I;
  int32_t step, nonstab;
  int32_t top, topb;
  uint8_t v[27];

  if (inhibit == NULL) 
  {
    inhibit = copyimage(image); 
    razimage(inhibit);
  }
  I = UCHARDATA(inhibit);

  if (nsteps == -1) nsteps = 1000000000;

  for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;

  init_topo3d();

  /* ================================================ */
  /*               DEBUT ALGO                         */
  /* ================================================ */

  step = 0;
  nonstab = 1;
  while (nonstab && (step < nsteps))
  {
    nonstab = 0;
    step++;
#ifdef VERBOSE
    printf("step %d\n", step);
#endif

    // AJOUTE LES EXTREMITES DANS I
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match_end(v)) I[i] = 1;
      }

    // MARQUE LES POINTS SIMPLES NON DANS I
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]) && !I[i] && simple26(S, i, rs, ps, N))
	SET_SIMPLE(S[i]);

    // MARQUE LES POINTS 2M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 1M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match1(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 0M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match0(v))
	  insert_vois(v, S, i, rs, ps, N);
      }

    memset(T, 0, N);
    for (i = 0; i < N; i++) // T := [S \ P] \cup M, où M représente les pts marqués
      if ((S[i] && !IS_SIMPLE(S[i])) || IS_2M_CRUCIAL(S[i]) || IS_1M_CRUCIAL(S[i]) || IS_0M_CRUCIAL(S[i]))
	T[i] = 1;
#ifdef DEBUG
writeimage(t,"_T");
#endif

    for (i = 0; i < N; i++)
      if (S[i] && !T[i]) 
      {
	S[i] = 0; 
	nonstab = 1; 
      }
    for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;
  }

#ifdef VERBOSE1
    printf("number of steps: %d\n", step);
#endif

  for (i = 0; i < N; i++) if (S[i]) S[i] = 255; // normalize values

  freeimage(t);
  termine_topo3d();
  return(1);
} /* lskelEK3() */

/* ==================================== */
int32_t lskelCK3a(struct xvimage *image, 
	     int32_t nsteps,
	     struct xvimage *inhibit)
/* ==================================== */
/*
Squelette symétrique curviligne
Algo CK3a données: S
Répéter jusqu'à stabilité
  C := points de courbe de S
  I := I \cup C
  P := voxels simples pour S et pas dans I
  C2 := voxels 2M-cruciaux (match2)
  C1 := voxels 1M-cruciaux (match1)
  C0 := voxels 0M-cruciaux (match0)
  P := P  \  [C2 \cup C1 \cup C0]
  S := S \ P
*/
#undef F_NAME
#define F_NAME "lskelCK3a"
{ 
  int32_t i, j, x, y, z;
  int32_t rs = rowsize(image);     /* taille ligne */
  int32_t cs = colsize(image);     /* taille colonne */
  int32_t ds = depth(image);       /* nb plans */
  int32_t ps = rs * cs;            /* taille plan */
  int32_t N = ps * ds;             /* taille image */
  uint8_t *S = UCHARDATA(image);      /* l'image de depart */
  struct xvimage *t = copyimage(image); 
  uint8_t *T = UCHARDATA(t);
  uint8_t *I;
  int32_t step, nonstab;
  int32_t top, topb;
  uint8_t v[27];

  if (inhibit == NULL) 
  {
    inhibit = copyimage(image); 
    razimage(inhibit);
  }
  I = UCHARDATA(inhibit);

  if (nsteps == -1) nsteps = 1000000000;

  for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;

  init_topo3d();

  /* ================================================ */
  /*               DEBUT ALGO                         */
  /* ================================================ */

  step = 0;
  nonstab = 1;
  while (nonstab && (step < nsteps))
  {
    nonstab = 0;
    step++;
#ifdef VERBOSE
    printf("step %d\n", step);
#endif

    // MARQUE LES POINTS SIMPLES NON DANS I
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]) && !I[i] && simple26(S, i, rs, ps, N))
	SET_SIMPLE(S[i]);
    // DEUXIEME SOUS-ITERATION : MARQUE LES POINTS DE COURBE (2)
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2s(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // TROISIEME SOUS-ITERATION : MARQUE LES POINTS DE COURBE (1)
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match1s(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS DE COURBE (3)
    for (i = 0; i < N; i++)
    {
      if (IS_OBJECT(S[i]) && !IS_SIMPLE(S[i]))
      {    
	top26(S, i, rs, ps, N, &top, &topb);
	if (top > 1) SET_CURVE(S[i]);
      }
    }
    // DEMARQUE PTS DE COURBE ET LES MEMORISE DANS I
    for (i = 0; i < N; i++)
    { 
      UNSET_2M_CRUCIAL(S[i]);
      if (IS_CURVE(S[i])) { UNSET_SIMPLE(S[i]); I[i] = 1; }
    }
    // MARQUE LES POINTS 2M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 1M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match1(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 0M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match0(v))
	  insert_vois(v, S, i, rs, ps, N);
      }

    memset(T, 0, N);
    for (i = 0; i < N; i++) // T := [S \ P] \cup M, où M représente les pts marqués
      if ((S[i] && !IS_SIMPLE(S[i])) || IS_2M_CRUCIAL(S[i]) || IS_1M_CRUCIAL(S[i]) || IS_0M_CRUCIAL(S[i]))
	T[i] = 1;
#ifdef DEBUG
writeimage(t,"_T");
#endif

    for (i = 0; i < N; i++)
      if (S[i] && !T[i]) 
      {
	S[i] = 0; 
	nonstab = 1; 
      }
    for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;
  }

#ifdef VERBOSE1
    printf("number of steps: %d\n", step);
#endif

  for (i = 0; i < N; i++) if (S[i]) S[i] = 255; // normalize values

  freeimage(t);
  termine_topo3d();
  return(1);
} /* lskelCK3a() */

/* ==================================== */
int32_t lskelCK3b(struct xvimage *image, 
	     int32_t nsteps,
	     struct xvimage *inhibit)
/* ==================================== */
/*
Squelette symétrique curviligne
Variante avec reconstruction des points extrémités
Algo CK3b données: S
Répéter jusqu'à stabilité
  C := points de courbe de S
  E := points extrémités de S
  C := C union [E inter gamma(C)] 
  I := I \cup C
  P := voxels simples pour S et pas dans I
  C2 := voxels 2M-cruciaux (match2)
  C1 := voxels 1M-cruciaux (match1)
  C0 := voxels 0M-cruciaux (match0)
  P := P  \  [C2 \cup C1 \cup C0]
  S := S \ P
*/
#undef F_NAME
#define F_NAME "lskelCK3b"
{ 
  int32_t i, j, k, x, y, z;
  int32_t rs = rowsize(image);     /* taille ligne */
  int32_t cs = colsize(image);     /* taille colonne */
  int32_t ds = depth(image);       /* nb plans */
  int32_t ps = rs * cs;            /* taille plan */
  int32_t N = ps * ds;             /* taille image */
  uint8_t *S = UCHARDATA(image);      /* l'image de depart */
  struct xvimage *t = copyimage(image); 
  uint8_t *T = UCHARDATA(t);
  struct xvimage *e = copyimage(image); 
  uint8_t *E = UCHARDATA(e);
  uint8_t *I;
  int32_t step, nonstab;
  int32_t top, topb;
  uint8_t v[27];

  if (inhibit == NULL) 
  {
    inhibit = copyimage(image); 
    razimage(inhibit);
  }
  I = UCHARDATA(inhibit);

  if (nsteps == -1) nsteps = 1000000000;

  for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;

  init_topo3d();

  /* ================================================ */
  /*               DEBUT ALGO                         */
  /* ================================================ */

  step = 0;
  nonstab = 1;
  while (nonstab && (step < nsteps))
  {
    nonstab = 0;
    step++;
#ifdef VERBOSE
    printf("step %d\n", step);
#endif

    // AJOUTE LES EXTREMITES DANS E
    memset(E, 0, N);
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match_end(v)) E[i] = 1;
      }
    // MARQUE LES POINTS SIMPLES NON DANS I
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]) && !I[i] && simple26(S, i, rs, ps, N))
	SET_SIMPLE(S[i]);
    // DEUXIEME SOUS-ITERATION : MARQUE LES POINTS DE COURBE (2)
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2s(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // TROISIEME SOUS-ITERATION : MARQUE LES POINTS DE COURBE (1)
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match1s(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS DE COURBE (3)
    for (i = 0; i < N; i++)
    {
      if (IS_OBJECT(S[i]) && !IS_SIMPLE(S[i]))
      {    
	top26(S, i, rs, ps, N, &top, &topb);
	if (top > 1) SET_CURVE(S[i]);
      }
    }

    // DEMARQUE PTS DE COURBE ET LES MEMORISE DANS I
    // AJOUTE AUX POINTS DE COURBE LEURS VOISINS QUI SONT DANS E
    for (i = 0; i < N; i++)
    { 
      UNSET_2M_CRUCIAL(S[i]);
      if (IS_CURVE(S[i])) 
      {
        for (k = 0; k < 26; k += 1)        /* parcourt les voisins en 26-connexite */
        {
          j = voisin26(i, k, rs, ps, N);
          if ((j != -1) && E[j])
	  {
	    UNSET_SIMPLE(S[j]); 
	    I[j] = 1; 
	  }
	}
	UNSET_SIMPLE(S[i]); 
	I[i] = 1; 
      }
    }
    // MARQUE LES POINTS 2M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 1M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match1(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 0M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match0(v))
	  insert_vois(v, S, i, rs, ps, N);
      }

    memset(T, 0, N);
    for (i = 0; i < N; i++) // T := [S \ P] \cup M, où M représente les pts marqués
      if ((S[i] && !IS_SIMPLE(S[i])) || IS_2M_CRUCIAL(S[i]) || IS_1M_CRUCIAL(S[i]) || IS_0M_CRUCIAL(S[i]))
	T[i] = 1;
#ifdef DEBUG
writeimage(t,"_T");
#endif

    for (i = 0; i < N; i++)
      if (S[i] && !T[i]) 
      {
	S[i] = 0; 
	nonstab = 1; 
      }
    for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;
  }

#ifdef VERBOSE1
    printf("number of steps: %d\n", step);
#endif

  for (i = 0; i < N; i++) if (S[i]) S[i] = 255; // normalize values

  freeimage(t);
  freeimage(e);
  termine_topo3d();
  return(1);
} /* lskelCK3b() */

/* ==================================== */
int32_t lskelCK3(struct xvimage *image, 
	     int32_t nsteps,
	     struct xvimage *inhibit)
/* ==================================== */
/*
Squelette symétrique curviligne, variante
Les points "candidats" à devenir des points de courbes sont les
points qui ne sont pas voisins d'un point isthme 2D ni d'un point interieur
Algo CK3 données: S
Répéter jusqu'à stabilité
  C := points de courbe de S
  P := voxels simples pour S et pas dans C
  C2 := voxels 2M-cruciaux (match2)
  C1 := voxels 1M-cruciaux (match1)
  C0 := voxels 0M-cruciaux (match0)
  P := P  \  [C2 \cup C1 \cup C0]
  S := S \ P
*/
#undef F_NAME
#define F_NAME "lskelCK3"
{ 
  int32_t i, j, k, x, y, z;
  int32_t rs = rowsize(image);     /* taille ligne */
  int32_t cs = colsize(image);     /* taille colonne */
  int32_t ds = depth(image);       /* nb plans */
  int32_t ps = rs * cs;            /* taille plan */
  int32_t N = ps * ds;             /* taille image */
  uint8_t *S = UCHARDATA(image);      /* l'image de depart */
  struct xvimage *t = copyimage(image); 
  uint8_t *T = UCHARDATA(t);
  uint8_t *I;
  int32_t step, nonstab;
  int32_t top, topb;
  uint8_t v[27];

  if (inhibit == NULL) 
  {
    inhibit = copyimage(image); 
    razimage(inhibit);
  }
  I = UCHARDATA(inhibit);

  if (nsteps == -1) nsteps = 1000000000;

  for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;

  init_topo3d();

  /* ================================================ */
  /*               DEBUT ALGO                         */
  /* ================================================ */

  step = 0;
  nonstab = 1;
  while (nonstab && (step < nsteps))
  {
    nonstab = 0;
    step++;
#ifdef VERBOSE
    printf("step %d\n", step);
#endif

    // MARQUE LES POINTS SIMPLES NON DANS I
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]) && !I[i] && simple26(S, i, rs, ps, N))
	SET_SIMPLE(S[i]);
    // MARQUE LES POINTS DE SURFACE (2)
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2s(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS DE SURFACE (3) ET LES POINTS INTERIEURS
    for (i = 0; i < N; i++)
    {
      if (IS_OBJECT(S[i]) && !IS_SIMPLE(S[i]))
      {    
	top26(S, i, rs, ps, N, &top, &topb);
	if (topb > 1) SET_SURF(S[i]);
	if (topb == 0) SET_SELECTED(S[i]);
      }
    }

    // DEMARQUE PTS ET REND "NON-SIMPLES" LES CANDIDATS
    for (i = 0; i < N; i++)
    { 
      UNSET_2M_CRUCIAL(S[i]);
      if (IS_OBJECT(S[i])) 
      {
#ifdef SAUVE
        for (k = 0; k < 12; k += 2)        /* parcourt les voisins en 6-connexite */
        {
          j = voisin6(i, k, rs, ps, N);
#endif
	for (k = 0; k < 26; k += 1)        /* parcourt les voisins en 26-connexite */
        {
	  j = voisin26(i, k, rs, ps, N);
          if ((j != -1) && IS_SELECTED(S[j])) break;
	}
#ifdef SAUVE
	if (k == 12) // le voxel est résiduel
#endif
	if (k == 26) // le voxel est résiduel
	{
	  for (k = 0; k < 26; k += 1)        /* parcourt les voisins en 26-connexite */
          {
	    j = voisin26(i, k, rs, ps, N);
	    if ((j != -1) && IS_SURF(S[j]))break;
	  }
	  if (k == 26) UNSET_SIMPLE(S[i]);
	}
      }
    }
    // MARQUE LES POINTS 2M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 1M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match1(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
    // MARQUE LES POINTS 0M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match0(v))
	  insert_vois(v, S, i, rs, ps, N);
      }

    memset(T, 0, N);
    for (i = 0; i < N; i++) // T := [S \ P] \cup M, où M représente les pts marqués
      if ((S[i] && !IS_SIMPLE(S[i])) || IS_2M_CRUCIAL(S[i]) || IS_1M_CRUCIAL(S[i]) || IS_0M_CRUCIAL(S[i]))
	T[i] = 1;
#ifdef DEBUG
writeimage(t,"_T");
#endif

    for (i = 0; i < N; i++)
      if (S[i] && !T[i]) 
      {
	S[i] = 0; 
	nonstab = 1; 
      }
    for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;
  }

#ifdef VERBOSE1
    printf("number of steps: %d\n", step);
#endif

  for (i = 0; i < N; i++) if (S[i]) S[i] = 255; // normalize values

  freeimage(t);
  termine_topo3d();
  return(1);
} /* lskelCK3() */

/* ==================================== */
int32_t lskelAK3(struct xvimage *image, 
	     int32_t nsteps,
	     struct xvimage *inhibit,
	     int32_t filter)
/* ==================================== */
/*
Amincissement symétrique avec inclusion de l'axe médian

A REVOIR : 
- pb des "échancrures"
- si on filtre l'axe médian alors l'homotopie n'est plus garantie ?
(cf. prop. sur les 0-cliques cruciales)

Algo AK3 données: S
K := \emptyset ; T := S
Répéter jusqu'à stabilité
  E := T \ominus \Gamma_6
  D := T \ [E \oplus \Gamma_6]
  T := E
  K := K \cup D
  P := voxels de S \ K simples pour S
  R := voxels de P qui s'apparient avec le masque C (bertrand_match3)
  S := [S  \  P]  \cup  R
*/
#undef F_NAME
#define F_NAME "lskelAK3"
{ 
  int32_t i, j, x, y, z;
  int32_t rs = rowsize(image);     /* taille ligne */
  int32_t cs = colsize(image);     /* taille colonne */
  int32_t ds = depth(image);       /* nb plans */
  int32_t ps = rs * cs;            /* taille plan */
  int32_t N = ps * ds;             /* taille image */
  uint8_t *S = UCHARDATA(image);      /* l'image de depart */
  struct xvimage *t = copyimage(image); 
  uint8_t *T = UCHARDATA(t);
  struct xvimage *r = copyimage(image); 
  uint8_t *R = UCHARDATA(r);
  struct xvimage *e = copyimage(image); 
  uint8_t *E = UCHARDATA(e);
  struct xvimage *d = copyimage(image); 
  uint8_t *D = UCHARDATA(d);
  struct xvimage *k = copyimage(image); 
  uint8_t *K = UCHARDATA(k);
  int32_t step, nonstab;
  uint8_t v[27];

  if (inhibit != NULL)
  {
    fprintf(stderr, "%s: inhibit image not implemented\n", F_NAME);
    return 0;
  }

  if (nsteps == -1) nsteps = 1000000000;

  for (i = 0; i < N; i++) if (S[i]) S[i] = 1; // normalize values

  init_topo3d();

  /* ================================================ */
  /*               DEBUT ALGO                         */
  /* ================================================ */

  step = 0;
  nonstab = 1;
  // K := \emptyset ; T := S
  memset(K, 0, N);
  memcpy(T, S, N);
  while (nonstab && (step < nsteps))
  {
    nonstab = 0;
    step++;
#ifdef VERBOSE
    printf("step %d\n", step);
#endif

    //  E := T \ominus \Gamma_6 
    memset(E, 0, N);
    for (i = 0; i < N; i++) 
      if (T[i] && (nbvoiso6(T, i, rs, ps, N) == 6)) E[i] = 1;

    //  D := E \oplus \Gamma_6
    memset(D, 0, N);
    for (i = 0; i < N; i++)
      if (E[i] || (nbvoiso6(E, i, rs, ps, N) >= 1)) D[i] = 1;

    //  D := T \ D
    for (i = 0; i < N; i++)
      if (T[i] && !D[i]) D[i] = 1; else D[i] = 0;

    //  T := E
    memcpy(T, E, N);

    //  K := K \cup D
    for (i = 0; i < N; i++)
      if (D[i]) K[i] = 1;

    // PREMIERE SOUS-ITERATION : MARQUE LES POINTS SIMPLES qui ne sont pas dans K
    if (step > filter)
    {
      for (i = 0; i < N; i++) 
	if (S[i] && !K[i] && simple26(S, i, rs, ps, N))
	  SET_SIMPLE(S[i]);
    }
    else
    {
      for (i = 0; i < N; i++) 
	if (S[i] && simple26(S, i, rs, ps, N))
	  SET_SIMPLE(S[i]);
    }

    // DEUXIEME SOUS-ITERATION : MARQUE LES CLIQUES CRUCIALES CORRESPONDANT AUX 2-FACES
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2(v))
	  insert_vois(v, S, i, rs, ps, N);
      }

    // TROISIEME SOUS-ITERATION : MARQUE LES CLIQUES CRUCIALES CORRESPONDANT AUX 1-FACES
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match1(v))
	  //if (match1s(v))  // VARIANTE POUR EVITER LES "ECHANCRURES" (à voir)
	  insert_vois(v, S, i, rs, ps, N);
      }

    // D := [S \ P] \cup  R, où R représente les pts marqués
    memset(D, 0, N);
    for (i = 0; i < N; i++) 
      if ((S[i] && !IS_SIMPLE(S[i])) || IS_2M_CRUCIAL(S[i]) || IS_1M_CRUCIAL(S[i]))
	D[i] = 1;

    for (i = 0; i < N; i++) // pour  tester la stabilité
      if (S[i] && !D[i]) 
      {
	S[i] = 0; 
	nonstab = 1; 
      }
    for (i = 0; i < N; i++) if (S[i]) S[i] = 1;
  }

#ifdef VERBOSE1
    printf("number of steps: %d\n", step);
#endif

  for (i = 0; i < N; i++) if (S[i]) S[i] = 255; // normalize values

  freeimage(t);
  freeimage(r);
  freeimage(e);
  freeimage(d);
  freeimage(k);
  termine_topo3d();
  return(1);
} /* lskelAK3() */

/* ==================================== */
int32_t lskelMK3(struct xvimage *image, 
	     int32_t nsteps,
	     struct xvimage *inhibit)
/* ==================================== */
/*
Squelette symétrique ultime avec ensemble de contrainte
Version révisée d'après le papier IWCIA 2006
Algo MK3 données: S
Répéter jusqu'à stabilité
  P := voxels simples pour S et non dans I
  C2 := voxels 2M-cruciaux (match2)
  C1 := voxels 1M-cruciaux (match1)
  C0 := voxels 0M-cruciaux (match0)
  P := P  \  [C2 \cup C1 \cup C0]
  S := S \ P
*/
#undef F_NAME
#define F_NAME "lskelMK3"
{ 
  int32_t i, j, x, y, z;
  int32_t rs = rowsize(image);     /* taille ligne */
  int32_t cs = colsize(image);     /* taille colonne */
  int32_t ds = depth(image);       /* nb plans */
  int32_t ps = rs * cs;            /* taille plan */
  int32_t N = ps * ds;             /* taille image */
  uint8_t *S = UCHARDATA(image);      /* l'image de depart */
  struct xvimage *t = copyimage(image); 
  uint8_t *T = UCHARDATA(t);
  uint8_t *I = NULL;
  int32_t step, nonstab;
  uint8_t v[27];

  if (nsteps == -1) nsteps = 1000000000;

  for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;

  if (inhibit != NULL) I = UCHARDATA(inhibit);

  init_topo3d();

  /* ================================================ */
  /*               DEBUT ALGO                         */
  /* ================================================ */

  step = 0;
  nonstab = 1;
  while (nonstab && (step < nsteps))
  {
    nonstab = 0;
    step++;
#ifdef VERBOSE
    printf("MK3b step %d\n", step);
#endif

    // PREMIERE SOUS-ITERATION : MARQUE LES POINTS SIMPLES ET PAS DANS I
    for (i = 0; i < N; i++) 
      if (IS_OBJECT(S[i]) && simple26(S, i, rs, ps, N) && (!I || !I[i]))
	SET_SIMPLE(S[i]);
#ifdef DEBUG
writeimage(image,"_S");
#endif
    // DEUXIEME SOUS-ITERATION : MARQUE LES POINTS 2M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
	extract_vois(S, i, rs, ps, N, v);
	if (match2(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
#ifdef DEBUG
memset(R, 0, N);
for (i = 0; i < N; i++) if (IS_2M_CRUCIAL(S[i])) R[i] = 255;
writeimage(r,"_C2");
#endif

    // TROISIEME SOUS-ITERATION : MARQUE LES POINTS 1M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
#ifdef DEBUG
printf("%d %d %d\n", i % rs, (i % ps) / rs, i / ps);
#endif
	extract_vois(S, i, rs, ps, N, v);
	if (match1(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
#ifdef DEBUG
memset(R, 0, N);
for (i = 0; i < N; i++) if (IS_1M_CRUCIAL(S[i])) R[i] = 255;
writeimage(r,"_C1");
#endif

    // QUATRIEME SOUS-ITERATION : MARQUE LES POINTS 0M-CRUCIAUX
    for (i = 0; i < N; i++) 
      if (IS_SIMPLE(S[i]))
      { 
#ifdef DEBUG
printf("%d %d %d\n", i % rs, (i % ps) / rs, i / ps);
#endif
	extract_vois(S, i, rs, ps, N, v);
	if (match0(v))
	  insert_vois(v, S, i, rs, ps, N);
      }
#ifdef DEBUG
memset(R, 0, N);
for (i = 0; i < N; i++) if (IS_0M_CRUCIAL(S[i])) R[i] = 255;
writeimage(r,"_C0");
#endif

    memset(T, 0, N);
    for (i = 0; i < N; i++) // T := [S \ P] \cup  R, où R représente les pts marqués
      if ((S[i] && !IS_SIMPLE(S[i])) || IS_2M_CRUCIAL(S[i]) || IS_1M_CRUCIAL(S[i]) || IS_0M_CRUCIAL(S[i]))
	T[i] = 1;
#ifdef DEBUG
writeimage(t,"_T");
#endif

    for (i = 0; i < N; i++)
      if (S[i] && !T[i]) 
      {
	S[i] = 0; 
	nonstab = 1; 
      }
    for (i = 0; i < N; i++) if (S[i]) S[i] = S_OBJECT;
  }

#ifdef VERBOSE1
    printf("number of steps: %d\n", step);
#endif

  for (i = 0; i < N; i++) if (S[i]) S[i] = 255; // normalize values

  freeimage(t);
  termine_topo3d();
  return(1);
} /* lskelMK3() */
