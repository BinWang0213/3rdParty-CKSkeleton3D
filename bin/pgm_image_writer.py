import numpy as np

def readPMG(fname):
    """ Read 3D image into binary numpy with unsigned char

    # This is basically used for input data of image skeleton tools M. Couprie, https://perso.esiee.fr/~coupriem/ck/index.html
    
    *.pgm format
    ---------
    P5/P7 #UCHAR code
    #xdim 1
    #ydim 1
    #zdim 1
    61 53 12 #NX,NY,NZ
    255 #Max value of UCHAR
    [binary raw image data.....]

    Author:Bin Wang(binwang.0213@gmail.com)
    Date: Jan. 2020
    """
    img=None
    with open(fname, mode='rb') as file: # b is important -> binary
        for li,line in enumerate(file):
            #Head
            if(li==0):
                DataType=line.decode('ascii').split()
                print("Only unsigned char is supported!")
                assert DataType or DataType, "only unsigned char is supported!"
            if(li==4): 
                NXYZ=[int(t) for t in line.decode('ascii').split()]
                NumVoxels=NXYZ[0]*NXYZ[1]*NXYZ[2]
                print("Image Size=",NXYZ,NumVoxels)
            #Raw data
            if(li==6):
                img=np.frombuffer(line,dtype=np.uint8)
                assert len(img)==NXYZ[0]*NXYZ[1]*NXYZ[2], "Wrong Data size %d->%d" % (len(img),NumVoxels)
                print("Image value range=",[np.min(img),np.max(img)])
    
    if(img is not None):
        img_3D=np.array(np.reshape(img,NXYZ,order='F'),dtype=np.uint8)
        return img_3D


def Image2PMG(fname,img):
    """ Write 3D image into binary PMG file with unsigned char

    # This is basically used for input data of image skeleton tools M. Couprie, https://perso.esiee.fr/~coupriem/ck/index.html
    
    *.pgm format
    ---------
    P5 #UCHAR code
    #xdim 1
    #ydim 1
    #zdim 1
    61 53 12 #NX,NY,NZ
    255 #Max value of UCHAR
    [binary raw image data.....]

    Author:Bin Wang(binwang.0213@gmail.com)
    Date: Jan. 2020
    """
    header=b"P7\n"
    header+=b"#xdim 1\n"
    header+=b"#ydim 1\n"
    header+=b"#zdim 1\n"
    header+=b"DATASET STRUCTURED_POINTS\n"
    header+=b"%d %d %d\n" % (img.shape[0],img.shape[1],img.shape[2])
    header+=b"255\n"
    with open(fname,'wb') as f:
        f.write(header)
        f.write(img.flatten(order="F").tobytes())