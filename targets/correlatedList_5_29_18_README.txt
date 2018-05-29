Hi David,

here is the latest version of the QSO-galaxy correlation table.
They are sorted by (1) classification (2) S/N ratio.


There are still some minor issues with some sightlines. E.g. two
STIS G140M sightlines (PG1149-110 and PG1103-006 are marked as not
completely ready, but I know I did them.
There are also a couple of good sightlines that I have not yet fully
IDed, but for now we'll skip them.


I classified each coincidence (see below).


I also realized that another way to approach this is to look at
each galaxy and find the nearby sightlines, so that we can also
include non-detections.


For now, I think that plan needs to be
(a) finish the SALT paper
(b) turn the catalogue paper into a thesis chapter
(c) finish the thesis as much as possible
(d) work on the last bit, which is the galaxy-QSO general study

We can discuss next week what we want to be in that last chapter

Bart









An entry consists of the following set of data.


For example the first one

MRK876                    QSO       98.27  40.38   0.12900      PI: Arav,Green       11524,11686   LB4F05,LB4Q03           SN=N:59.7,C:41.8,O:33.8 v3.0        nightonly=TODO  IGMPARS=yes FPN=yes

1st line: general data: name, type, galactic lon/lat, redshift
program PI, program ID (11524,11686), spectrum ID (LB4F05),
S/N ratio at 1238 (N:59.7), at 1548 (C:41.8), at 1031 (O:33.8),
calcos version (v3.0), whether nightonly data was made,
IGMPARS flag (yes or no)
and finally FPN flag. If this is 'yes', the sightline passed
all consistency tests I have (if 'STIS' or 'FUSE' in capital letters
there is no COS data).
If FPN=no/stis/fuse I still need to deal with it.

In the case of SBS0957+599, SBS1122+594, FBQSJ0751+2919
the FPN= flag is not yes/no, because
I need to fix the ID, but the Lyman alpha lines below 10000 km/s
are fine, so you can use these ones.
I will make finalizing them a priority, and if something changes
for the cz<10000 systems I'll tell you. I don't expect that to
happen.



2nd: set of lines with the IGM lines at cz<10000, both Lyman alpha
and metal lines.

name                                     wavel   classific.   line      lambd0     cz     z   helio/lst                 EW        dEW          system#
MRK876                    IGMpars wakk  1255.796 IGMABS       Lya       1215.6    9895  0.03301 H  1255.516 1255.977    21.0 pm   2.3 pm   0.0 10 
MRK876                    IGMpars wakk  1244.076 IGMABS       Lya       1215.6    7005  0.02337 H  1243.759 1244.379    19.8 pm   2.7 pm   0.0  9 
MRK876                    IGMpars wakk  1240.150 IGMABS       Lya       1215.6    6037  0.02014 H  1239.953 1240.308    85.1 pm   1.9 pm   0.0  8 
MRK876                    IGMpars wakk  1236.093 IGMABS       Lya       1215.6    5036  0.01680 H  1235.924 1236.279    20.7 pm   2.1 pm   0.0  7 
MRK876                    IGMpars wakk  1233.949 IGMABS       Lya       1215.6    4508  0.01504 H  1233.671 1233.990    17.0 pm   2.0 pm   0.0  6 
MRK876                    IGMpars wakk  1229.774 IGMABS       Lya       1215.6    3478  0.01160 H  1229.531 1230.034   230.0 pm   2.0 pm   0.0  5 
MRK876                    IGMpars wakk  1037.626 IGMABS       Lyb       1025.7    3479  0.01160 H  1037.482 1037.770    43.4 pm   4.6 pm   0.0  5 
MRK876                    IGMpars wakk  1220.562 IGMABS       SiIII     1206.5    3494  0.01160 H  1220.468 1220.638     5.0 pm   1.7 pm   0.0  5 
MRK876                    IGMpars wakk  1219.478 IGMABS       Lya       1215.6     939  0.00313 H  1219.077 1219.743   345.1 pm   2.9 pm   0.0  4 
MRK876                    IGMpars wakk  1028.913 IGMABS       Lyb       1025.7     933  0.00313 H  1028.684 1029.124    53.1 pm   4.4 pm   0.0  4 
MRK876                    IGMpars wakk  1035.176 IGMABS       OVI       1031.9     944  0.00313 H  1035.007 1035.305    20.3 pm   3.6 pm   0.0  4 


3rd: line given reason for classification

Classification for QSO-Gal correlation =   cl=119 @ 3504


cl=119
1st number = number of galaxies with L>0.1 with velocity within 500 km/s of
   the galaxy checked (the one with v=3504).
   This is looked at separately for each galaxy with L>0.1
2nd number = number of galaxies with L>0 within 500 km/s
3rd number = ratio of galaxy with highest L to galaxy with 2nd highest L
   (within 500 km/s)



4rd detailed data on the galaxies.


target                     RA      Dec   type      z      galaxy                          RA      Dec   type      vgal inc  pa  dist method  diam Rvir      lum   imp likel sta S/N  det diff L(v)
MRK876                   243.488  65.719 Sey1     0.12900 UGC10294                      243.342  65.418 Sm?       3504  49   4  51.1 HL     30.23 181.2   0.000   274 0.102 SI  59.7 LYA  26 0.100 CL=119
MRK876                   243.488  65.719 Sey1     0.12900 IC1214                        244.048  65.969 S0/a      7013  59  20 100.5 HL     38.73 232.2   1.500   594 0.001 SI  59.7 LYA   8 0.001 CL=119
MRK876                   243.488  65.719 Sey1     0.12900 NGC6140                       245.242  65.391 SB(s)cd_   910  44  60   8.1 TF      6.72  40.3   0.200   113 0.000 SI  59.7 LYA  29 0.000 CL=119


So, the coincidences to look at are the ones with cl=1#r, with r>5
(i.e. there is just one galaxy with L>0.1, and any number of other
galaxies with L>0 within 500 km/s, but the ratio of highest L to
next higher L is >5.



