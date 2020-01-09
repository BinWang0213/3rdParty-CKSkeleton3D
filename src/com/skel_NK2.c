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
/* $Id:  $ */
/*! \file skel_NK2.c

\brief parallel 2D binary curvilinear skeleton

<B>Usage:</B> skel_NK2 in.pgm nsteps [inhibit] out.pgm

<B>Description:</B>
Parallel 2D binary thinning or curvilinear skeleton. The parameter \b nsteps gives,
if positive, the number of parallel thinning steps to be processed.
If the value given for \b nsteps equals -1, the thinning is continued
until stability.

If the parameter \b inhibit is given and is a binary image name,
then the points of this image will be left unchanged. 

Reference:
G. Bertrand and M. Couprie,
"Two-dimensional thinning algorithms based on critical kernels", 
<I> Journal of Mathematical Imaging and Vision</I>, submitted, 2006.
<A HREF="http://igm.univ-mlv.fr/LabInfo/rapportsInternes/2006/02.pdf">Preprint: IGM2006-02</A>.

<B>Types supported:</B> byte 2d

<B>Category:</B> topobin
\ingroup  topobin

\author Michel Couprie
*/
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <mccodimage.h>
#include <mcimage.h>
#include <lskelpar.h>

/* =============================================================== */
int32_t main(argc, argv) 
/* =============================================================== */
  int32_t argc; char **argv; 
{
  struct xvimage * image;
  struct xvimage * inhibit = NULL;
  int32_t ret, nsteps;
  int32_t filter = 0;

  if ((argc != 4) && (argc != 5))
  {
    fprintf(stderr, "usage: %s in.pgm nsteps [inhibit] out.pgm\n", argv[0]);
    exit(1);
  }

  image = readimage(argv[1]);
  if (image == NULL)
  {
    fprintf(stderr, "%s: readimage failed\n", argv[0]);
    exit(1);
  }

  nsteps = atoi(argv[2]);
  if (argc == 5)
  {
    inhibit = readimage(argv[3]);
    if (inhibit == NULL)
    {
      fprintf(stderr, "%s: readimage failed\n", argv[0]);
      exit(1);
    }
  }

  if (depth(image) == 1)
  {
    if (! lskelNK2(image, nsteps, inhibit))
    {
      fprintf(stderr, "%s: lskelNK2 failed\n", argv[0]);
      exit(1);
    } 
  }
  else
  {
    fprintf(stderr, "%s: image must be 2D\n", argv[0]);
    exit(1);
  }

  writeimage(image, argv[argc-1]);
  freeimage(image);

  return 0;
} /* main */
