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
/* $Id: pgm2raw.c,v 1.7 2006/12/01 14:41:45 michel Exp $ */
/*! \file pgm2raw.c

\brief suppress the header from a pgm file

<B>Usage:</B> pgm2raw in.pgm out.raw

<B>Description:</B> suppress the header from a pgm file

<B>Types supported:</B> byte 2d, byte 3d

<B>Category:</B> convert
\ingroup  convert

\author Michel Couprie
*/

/* Michel Couprie - janvier 2000 */

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <mcimage.h>
#include <mccodimage.h>

#define VERBOSE

/* =============================================================== */
int32_t main(argc, argv) 
/* =============================================================== */
  int32_t argc; char **argv; 
{
  FILE *fd = NULL;
  int32_t rs, cs, ds, N, ret;
  struct xvimage * image;

  if (argc != 3)
  {
    fprintf(stderr, "usage: %s in.pgm out.raw \n", argv[0]);
    exit(1);
  }

  image = readimage(argv[1]);
  if (image == NULL)
  {
    fprintf(stderr, "%s: readimage failed\n", argv[0]);
    exit(1);
  }

  if (datatype(image) != VFF_TYP_1_BYTE)
  {
    fprintf(stderr, "%s: only byte images supported\n", argv[0]);
    exit(1);
  }

  rs = rowsize(image);
  cs = colsize(image);
  ds = depth(image);
  N = rs * cs * ds;

#ifdef VERBOSE
  printf("rs = %d ; cs = %d ; ds = %d ; N = rs * cs * ds = %d\n", rs, cs, ds, N);
#endif  

#ifdef UNIXIO
  fd = fopen(argv[argc - 1],"w");
#endif
#ifdef DOSIO
  fd = fopen(argv[argc - 1],"wb");
#endif

  ret = fwrite(UCHARDATA(image), sizeof(char), N, fd);
  if (ret != N)
  {
    fprintf(stderr, "%s: only %d items written\n", argv[0], ret);
    exit(1);
  }

  fclose(fd);

  return 0;
} /* main */


