memory,50,m
nocompress;
geomtyp=xyz
angstrom
geometry={
O	0.0	0.0	-0.0657441568
H	0.0	0.7574590974000001	0.5217905143000001
H	0.0	-0.7574590973999998	0.4717905143000001
}
basis=cc-pVTZ-F12
set,charge=0
set,spin=0
hf
{CCSD(T)-F12}
