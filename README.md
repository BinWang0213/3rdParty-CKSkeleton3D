# 3rdParty-CKSkeleton3D
This is a windows version code (compiled using MinGW_64 come with Code Blocks) of IMA3D @ https://perso.esiee.fr/~coupriem/ck/CK_programs.html

More skeleton code can be found @ https://perso.esiee.fr/~coupriem/ck/index.html

Reference paper
G. Bertrand and M. Couprie, "New 3D parallel thinning algorithms based on critical kernels", Discrete geometry for computer imagery, Lecture Notes in Computer Science, Vol. 4245, pp. 580-591, Springer, 2006.

We provide C source code for several algorithms described in  [BC06a,BC06b,BC06c], namely MK², AK², NK² (2D), MK³, EK³ and CK³ (3D).

* FILE FORMAT

We use a standard file format for 2D images, called pgm. The format description can be found here. For 3D images, we use a "natural" extension of the pgm format. In addition to the "width" and "height" fields in the file header, a "depth" field (in ASCII decimal) indicates the number of planes of the 3D volume. The data is stored as a succession of 2D planes. Both raw data and ASCII data are supported.

For conversion from and to raw file format, we provide the programs raw2pgm and pgm2raw.

* TEST

To test the programs, you may use the images provided in the directory "test". These images were also used as illutrations in the papers [BC06a,BC06b,BC06c].

* DOCUMENTATION

A quite minimal documentation is provided. The html version can be found in "doc/html/index.html", and the pdf version in "doc/latex/refman.pdf".

* EXAMPLE - Surface skeleton extraction

Usage: skel_MK3 in.pgm nsteps [inhibit] out.pgm

Description: Parallel 3D binary thinning or ultimate skeleton. The parameter nsteps gives, if positive, the number of parallel thinning steps to be processed. If the value given for nsteps equals -1, the thinning is continued until stability.

If the parameter inhibit is given and is a binary image name, then the points of this image will be left unchanged.

Types supported: byte (unsigned char) 

Category: topobin

Author:
Michel Couprie