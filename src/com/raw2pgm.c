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
/* $Id: raw2pgm.c,v 1.9 2006/12/04 14:33:57 michel Exp $ */
/*! \file raw2pgm.c

\brief converts from raw format into pgm format

<B>Usage:</B> in.raw rs cs ds headersize nbytespervox littleendian [xdim ydim zdim] out.pgm

<B>Description:</B> Converts from raw format into pgm format. 

Parameters:
\li \b in.pgm : source file in raw format
\li \b rs (int32_t): row size (number of voxels in a row)
\li \b cs (int32_t): column size (number of voxels in a column)
\li \b ds (int32_t): number of planes
\li \b headersize (int32_t): size of the header in bytes (information in the header will be ignored)
\li \b nbytespervox (int32_t): number of bytes per voxel (1, 2 or 4)
\li \b littleendian (int32_t) 1: littleendian, 0: bigendian. Usual choice is 0.
\li \b xdim (float, optional) : gap (in the real world) between two adjacent voxels in a row.  
\li \b ydim (float, optional) : gap (in the real world) between two adjacent voxels in a column.  
\li \b zdim (float, optional) : gap (in the real world) between two adjacent planes.  

<B>Types supported:</B> byte 3D, int16_t 3D, int32_t 3D

\warning Signed integers are not supported.

<B>Category:</B> convert
\ingroup convert

\author Michel Couprie
*/

/* Michel Couprie - janvier 2000 */

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <mcimage.h>
#include <mccodimage.h>
#include <mcutil.h>

#define VERBOSE

/* =============================================================== */
int32_t main(argc, argv) 
/* =============================================================== */
  int32_t argc; char **argv; 
{
  FILE *fd = NULL;
  int32_t rs, cs, ds, N, nbytesperpix, headersize, littleendian;
  struct xvimage * image;
  uint16_t tmp;
  uint8_t tmp1;
  uint32_t * I;
  int32_t i;

  if ((argc != 9) && (argc != 12))
  {
    fprintf(stderr, "usage: %s in.raw rs cs ds headersize nbytespervox littleendian [xdim ydim zdim] out.pgm \n", argv[0]);
    exit(1);
  }

#ifdef UNIXIO
  fd = fopen(argv[1],"r");
#endif
#ifdef DOSIO
  fd = fopen(argv[1],"rb");
#endif

  if (fd == NULL)
  {
    fprintf(stderr, "%s: fopen failed for %s\n", argv[0], argv[1]);
    exit(1);
  }

  rs = atoi(argv[2]);
  cs = atoi(argv[3]);
  ds = atoi(argv[4]);
  N = rs * cs * ds;
  headersize = atoi(argv[5]);
  nbytesperpix = atoi(argv[6]);
  littleendian = atoi(argv[7]);

  if ((nbytesperpix != 1) && (nbytesperpix != 2) && (nbytesperpix != 4))
  {
    fprintf(stderr, "%s: bad value for pix size: %d\n", argv[0], nbytesperpix);
    exit(1);
  }

  if (nbytesperpix == 1)
  {
    image = allocimage(NULL, rs, cs, ds, VFF_TYP_1_BYTE);
    if (image == NULL)
    {   fprintf(stderr,"%s : allocimage failed\n", argv[0]);
        exit(1);
    }
    fread(UCHARDATA(image), sizeof(char), headersize, fd);
    // sera ecrase par la suite - attention: plante si la taille image est
    // < headersize
    fread(UCHARDATA(image), sizeof(char), N, fd);
  }
  else if (nbytesperpix == 2)
  {
    image = allocimage(NULL, rs, cs, ds, VFF_TYP_4_BYTE);
    if (image == NULL)
    {   fprintf(stderr,"%s : allocimage failed\n", argv[0]);
        exit(1);
    }
    fread(UCHARDATA(image), sizeof(char), headersize, fd);
    // sera ecrase par la suite - attention: plante si la taille image est
    // < headersize
    I = ULONGDATA(image);
    for (i = 0; i < N; i++) 
    { 
      fread(&tmp, 2 * sizeof(char), 1, fd);
      if (littleendian)
      {
        tmp1 = tmp & 0x00ff;
        tmp = tmp >> 8;
        tmp = tmp | (((uint32_t)tmp1) << 8);      
      }
      I[i] = (uint32_t)tmp; 
    }
  }
  else if (nbytesperpix == 4)
  {
    image = allocimage(NULL, rs, cs, ds, VFF_TYP_4_BYTE);
    if (image == NULL)
    {   fprintf(stderr,"%s : allocimage failed\n", argv[0]);
        exit(1);
    }
    fread(UCHARDATA(image), sizeof(char), headersize, fd);
    // sera ecrase par la suite - attention: plante si la taille image est
    // < headersize
    fread(ULONGDATA(image), 4 * sizeof(char), N, fd);
  }

  if (argc == 12)
  {
    image->xdim = atof(argv[8]);
    image->ydim = atof(argv[9]);
    image->zdim = atof(argv[10]);
  }

  writeimage(image, argv[argc - 1]);
  freeimage(image);
  fclose(fd);

  return 0;
} /* main */


