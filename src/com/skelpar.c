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
/*! \file skelpar.c

\brief parallel binary skeleton

<B>Usage:</B> skelpar in.pgm algorithm nsteps [inhibit] out.pgm

<B>Description:</B>
Parallel binary thinning or skeleton. The parameter \b nsteps gives,
if positive, the number of parallel thinning steps to be processed.
If the value given for \b nsteps equals -1, the thinning is continued
until stability.

The parameter \b algorithm is a numerical code
indicating which method will be used for the thinning.
The possible choices are:
\li 0: Pavlidis 1981
\li 1: Eckhardt and Maderlechner 1993
\li 2: Couprie and Bertrand, ultimate (MK2)
\li 3: Couprie and Bertrand, curvilinear
\li 4: Couprie and Bertrand, curvilinear with reconstruction property (AK2)
\li 5: Rutovitz 1966 [DOES NOT PRESERVE TOPOLOGY]
\li 6: Zhang and Wang 1988 [DOES NOT PRESERVE TOPOLOGY]
\li 7: Han, La and Rhee
\li 8: Guo and Hall 1992 (1)
\li 9: Guo and Hall 1992 (2)
\li 10: Guo and Hall 1992 (3)
\li 11: Chin and Wan and Stover and Iverson 1987
\li 12: Jang and Chin 1992 [DOES NOT PRESERVE TOPOLOGY]
\li 13: Jang and Chin (correction mc)
\li 14: MNS preservation
\li 16: Choy, Choy and Siu 1995 [DOES NOT PRESERVE TOPOLOGY]
\li 17: Bernard and Manzanera 1999
\li 18: Holt et al. 1987
\li 19: Hall 1989
\li 20: Wu and Tsai 1992 [DOES NOT PRESERVE TOPOLOGY]
\li 21: Manzanera and Bernard (variant by GB)
\li 23: Bertrand curvilinear with reconstruction, symmetrical
\li 24: Rosenfeld directionnel
\li 25: Experimental

If the parameter \b inhibit is given and is a binary image name,
then the points of this image will be left unchanged. 

<B>Warning:</B> The object must not have any point on the frame of the image.

References:<BR> 
[COU05] M. Couprie, Note on fifteen 2D parallel thinning algorithms,
IGM2006-01}, Universit\'e de Marne-la-Vall\'ee, 2006,
http://igm.univ-mlv.fr/LabInfo/rapportsInternes/2006/01.pdf<BR>

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
int main(int argc, char **argv)
/* =============================================================== */
{
  struct xvimage * image;
  struct xvimage * inhibit = NULL;
  int32_t mode, nsteps;

  if ((argc != 5) && (argc != 6))
  {
    fprintf(stderr, "usage: %s in.pgm algorithm nsteps [inhibit] out.pgm\n", argv[0]);
    fprintf(stderr, "  0: Pavlidis 1981\n");
    fprintf(stderr, "  1: Eckhardt and Maderlechner 1993\n");
    fprintf(stderr, "  2: Couprie and Bertrand, ultimate (MK2)\n");
    fprintf(stderr, "  3: Couprie and Bertrand, curvilinear\n");
    fprintf(stderr, "  4: Couprie and Bertrand, curvilinear with reconstruction property (AK2)\n");
    fprintf(stderr, "  5: Rutovitz 1966 [DOES NOT PRESERVE TOPOLOGY]\n");
    fprintf(stderr, "  6: Zhang and Wang 1988 [DOES NOT PRESERVE TOPOLOGY]\n");
    fprintf(stderr, "  7: Han, La and Rhee\n");
    fprintf(stderr, "  8: Guo and Hall 1992 (1)\n");
    fprintf(stderr, "  9: Guo and Hall 1992 (2)\n");
    fprintf(stderr, "  10: Guo and Hall 1992 (3)\n");
    fprintf(stderr, "  11: Chin and Wan and Stover and Iverson 1987\n");
    fprintf(stderr, "  12: Jang and Chin 1992 [DOES NOT PRESERVE TOPOLOGY]\n");
    fprintf(stderr, "  13: Jang and Chin (correction mc)\n");
    fprintf(stderr, "  14: MNS preservation\n");
    fprintf(stderr, "  15: Jang and Chin 1993 (reconstruction)\n");
    fprintf(stderr, "  16: Choy, Choy and Siu 1995 [DOES NOT PRESERVE TOPOLOGY]\n");
    fprintf(stderr, "  17: Bernard and Manzanera 1999\n");
    fprintf(stderr, "  18: Holt et al. 1987\n");
    fprintf(stderr, "  19: Hall 1989\n");
    fprintf(stderr, "  20: Wu and Tsai 1992 [DOES NOT PRESERVE TOPOLOGY]\n");
    fprintf(stderr, "  21: Manzanera and Bernard (variant by GB)\n");
    fprintf(stderr, "  22: Bertrand curvilinear with reconstruction, asymmetrical\n");
    fprintf(stderr, "  23: Bertrand curvilinear with reconstruction, symmetrical\n");
    fprintf(stderr, "  24: Rosenfeld directionnel\n");
    exit(1);
  }

  image = readimage(argv[1]);
  if (image == NULL)
  {
    fprintf(stderr, "%s: readimage failed\n", argv[0]);
    exit(1);
  }

  mode = atoi(argv[2]);
  nsteps = atoi(argv[3]);

  if (argc == 6)
  {
    inhibit = readimage(argv[4]);
    if (inhibit == NULL)
    {
      fprintf(stderr, "%s: readimage failed\n", argv[0]);
      exit(1);
    }
  }

  if (depth(image) == 1)
  {
    switch (mode)
    {
    case 0:
      if (! lskelpavlidis(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelpavlidis failed\n", argv[0]);
	exit(1);
      } break;
    case 1:
      if (! lskeleckhardt(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskeleckhardt failed\n", argv[0]);
	exit(1);
      } break;
    case 2:
      //      if (! lskelmcultime(image, nsteps, inhibit))
      if (! lskelMK2(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelMK2 failed\n", argv[0]);
	exit(1);
      } break;
    case 3:
      if (! lskelmccurv(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelmccurv failed\n", argv[0]);
	exit(1);
      } break;
    case 4:
      //      if (! lskelmccurvrec(image, nsteps, inhibit))
      if (! lskelAK2(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelAK2 failed\n", argv[0]);
	exit(1);
      } break;
    case 5:
      if (! lskelrutovitz(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelrutovitz failed\n", argv[0]);
	exit(1);
      } break;
    case 6:
      if (! lskelzhangwang(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelzhangwang failed\n", argv[0]);
	exit(1);
      } break;
    case 7:
      if (! lskelhanlarhee(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelhanlarhee failed\n", argv[0]);
	exit(1);
      } break;
    case 8:
      if (! lskelguohall(image, nsteps, inhibit, 1))
      {
	fprintf(stderr, "%s: lskelguohall failed\n", argv[0]);
	exit(1);
      } break;
    case 9:
      if (! lskelguohall(image, nsteps, inhibit, 2))
      {
	fprintf(stderr, "%s: lskelguohall failed\n", argv[0]);
	exit(1);
      } break;
    case 10:
      if (! lskelguohall(image, nsteps, inhibit, 3))
      {
	fprintf(stderr, "%s: lskelguohall failed\n", argv[0]);
	exit(1);
      } break;
    case 11:
      if (! lskelchinwan(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelchinwan failed\n", argv[0]);
	exit(1);
      } break;
    case 12:
      if (! lskeljang(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskeljang failed\n", argv[0]);
	exit(1);
      } break;
    case 13:
      if (! lskeljangcor(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskeljangcor failed\n", argv[0]);
	exit(1);
      } break;
    case 14:
      if (! lskelmns(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelmns failed\n", argv[0]);
	exit(1);
      } break;
    case 16:
      if (! lskelchoy(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelchoy failed\n", argv[0]);
	exit(1);
      } break;
    case 17:
      if (! lskelmanz(image, nsteps, inhibit, 1))
      {
	fprintf(stderr, "%s: lskelmanz failed\n", argv[0]);
	exit(1);
      } break;
    case 18:
      if (! lskelhall(image, nsteps, inhibit, 1))
      {
	fprintf(stderr, "%s: lskelhall failed\n", argv[0]);
	exit(1);
      } break;
    case 19:
      if (! lskelhall(image, nsteps, inhibit, 2))
      {
	fprintf(stderr, "%s: lskelhall failed\n", argv[0]);
	exit(1);
      } break;
    case 20:
      if (! lskelwutsai(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelwutsai failed\n", argv[0]);
	exit(1);
      } break;
    case 21:
      if (! lskelmanz(image, nsteps, inhibit, 2))
      {
	fprintf(stderr, "%s: lskelmanz failed\n", argv[0]);
	exit(1);
      } break;
    case 23:
      if (! lskelbertrand_sym(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelbertrand_sym failed\n", argv[0]);
	exit(1);
      } break;
    case 24:
      if (! lskelrosenfeld(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelbertrand_asym_s failed\n", argv[0]);
	exit(1);
      } break;
    case 25:
      if (! lskelbertrand_asym_s(image, nsteps, inhibit))
      {
	fprintf(stderr, "%s: lskelbertrand_asym_s failed\n", argv[0]);
	exit(1);
      } break;
    default: 
      fprintf(stderr, "%s: mode %d not implemented\n", argv[0], mode);
      exit(1);
    }
  }
  else
  {
    fprintf(stderr, "%s: 3D not implemented\n", argv[0]);
    exit(1);
  }

  writeimage(image, argv[argc-1]);
  freeimage(image);

  return 0;
} /* main */
