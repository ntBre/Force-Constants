memory,50,m
nocompress;
geomtyp=xyz
angstrom
geometry={
O	0.0	0.0	-0.015744156799999992
H	0.0	0.7574590974	0.5217905143
H	0.0	-0.7074590973999999	0.5217905143
}
basis=cc-pVTZ-F12
set,charge=0
set,spin=0
hf
{CCSD(T)-F12}
