include("./Utilities.jl")

function firstderivative(labels, coords, delta, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    tempcoords = copy(coords)
    for atom in 1:length(coords)
        for coord in 1:length(coords[atom])
            run(`mkdir Input`)
            # generating filenames
            filenum = parse(Int, string(atom) * string(coord))
            comfilename = "input$(filenum).com"
            # increment the point
            tempcoords[atom][coord] += delta
            # write files
            writecom(comfilename, labels, tempcoords, 50, test)
            writepbs(filenum, test)
            run(`mv $(comfilename) ./Input/.`)
            run(`mv mp$(filenum).pbs ./Input/.`)
            # decrement by twice delta to get
            # the negative steps
            tempcoords[atom][coord] -= 2*delta
            # update filenames
            comfilename = "input-$(filenum).com"
            # write more files
            writecom(comfilename, labels, tempcoords, 50, test)
            writepbs(-filenum, test)
            run(`mv $(comfilename) ./Input/.`)
            run(`mv mp-$(filenum).pbs ./Input/.`)
            # return the coordinates to their
            # original state
            tempcoords[atom][coord] += delta
        end
    end
end

labels, coords = readxyz("geom.xyz")
firstderivative(labels, coords, 0.05, false)
