"""Utility functions for force constants"""

function readxyz(filename, comment=true)
    """Read .xyz file with name filename and return two arrays, one of
            the atom labels and one of their coordinates. If comment is true,
            skip the number of atoms and comment lines. If it is false, only skip
            the number of atoms"""
    infile = open(filename, "r")
    lines = readlines(infile)
    # Skip the number of atoms and comment lines
    atom_labels = []
    atom_coords = []
    if comment
        start = 3
    else
        start = 2
    end
    for line in lines[start:end]
        splitline = split(line)
        push!(atom_labels, splitline[1])
        push!(atom_coords, map(x->parse(Float64, x), splitline[2:end]))
    end
    return atom_labels, atom_coords
end

function printxyz(atom_labels, atom_coords)
    if length(atom_labels) == length(atom_coords)
        for atom in 1:length(atom_labels)
            println(atom_labels[atom], "\t", join(map(string, atom_coords[atom]), "\t"))
        end
    end
end

function writecom(filename, atomlabels, atomcoords, memory=50, test=true)
        """Take a file name, arrays of atom labels and atom coords, 
            and optionally an amount of memory and write a .com file"""
    if test
        println("file: ", filename)
        println("memory,$(memory),m")
        println("nocompress;")
        println("geomtyp=xyz")
        println("angstrom")
        println("geometry={")
        printxyz(atomlabels, atomcoords)
        println("}")
        println("basis=cc-pVTZ-F12")
        println("set,charge=0")
        println("set,spin=0")
        println("hf")
        println("{CCSD(T)-F12}")
    else
        open(filename, "w") do infile
        write(infile, "memory,$(memory),m\n")
        write(infile, "nocompress;\n")
        write(infile, "geomtyp=xyz\n")
        write(infile, "angstrom\n")
        write(infile, "geometry={\n")
        for atom in 1:length(atomlabels)
            write(infile, atomlabels[atom], "\t", join(map(string, atomcoords[atom]), "\t"), "\n")
        end
        write(infile, "}\n")
        write(infile, "basis=cc-pVTZ-F12\n")
        write(infile, "set,charge=0\n")
        write(infile, "set,spin=0\n")
        write(infile, "hf\n")
        write(infile, "{CCSD(T)-F12}\n")
        end
    end
end

#writecom("test.com", labels, coords, 50, false)

function writepbs(filenum, test=true)
        """Take a file number and write a .pbs file"""
    if test
        println("file: mp$(filenum).pbs ")
        println("#!/bin/sh\n")
        println("#PBS -N job$(filenum)\n")
        println("#PBS -S /bin/bash\n")
        println("#PBS -j oe\n")
        println("#PBS -W umask=022\n")
        println("#PBS -l walltime=00:30:00\n")
        println("#PBS -l ncpus=1\n")
        println("#PBS -l mem=50mb\n\n")
        println("module load intel\n")
        println("module load mvapich2\n")
        println("module load pbspro\n")
        println("export PATH=/usr/local/apps/molpro/2015.1.35/bin:\$PATH\n\n")
        println("export WORKDIR=\$PBS_O_WORKDIR\n")
        println("export TMPDIR=/tmp/\$USER/\$PBS_JOBID\n")
        println("cd \$WORKDIR\n")
        println("mkdir -p \$TMPDIR\n\n")
        println("date\n")
        println("molpro -t 1 input$(filenum).out\n")
        println("date\n\n")
        println("rm -rf \$TMPDIR")
    else
        open("mp$(filenum).pbs", "w") do infile
            write(infile, "#!/bin/sh\n")
            write(infile, "#PBS -N job$(filenum)\n")
            write(infile, "#PBS -S /bin/bash\n")
            write(infile, "#PBS -j oe\n")
            write(infile, "#PBS -W umask=022\n")
            write(infile, "#PBS -l walltime=00:30:00\n")
            write(infile, "#PBS -l ncpus=1\n")
            write(infile, "#PBS -l mem=50mb\n\n")
            write(infile, "module load intel\n")
            write(infile, "module load mvapich2\n")
            write(infile, "module load pbspro\n")
	    write(infile, "export PATH=/usr/local/apps/molpro/2015.1.35/bin:\$PATH\n\n")
            write(infile, "export WORKDIR=\$PBS_O_WORKDIR\n")
            write(infile, "export TMPDIR=/tmp/\$USER/\$PBS_JOBID\n")
            write(infile, "cd \$WORKDIR\n")
            write(infile, "mkdir -p \$TMPDIR\n\n")
            write(infile, "date\n")
            write(infile, "molpro -t 1 input$(filenum).com\n")
            write(infile, "date\n\n")
            write(infile, "rm -rf \$TMPDIR")
        end
    end
end
