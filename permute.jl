include("./Utilities.jl")
using Combinatorics

MAXFILES = 250

function getenergies(labels, coordarray, delta, transforms, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    rm("./Input", recursive=true)
    mkdir("Input")
    refjob = makejob(0, labels, coordarray, 50, test)
    refout = submitjob(refjob, 0)
    toread = []
    push!(toread, [0, refout])
    energies = []
    for i in 1:length(transforms)
        filenum = makefilenum(transforms[i])
        coords = coordarray + delta*reshape(transforms[i], (length(labels), 3))
        outfile = submitjob(makejob(filenum, labels, coords, 50, test), filenum)
        println(outfile, "\n")
        push!(toread, [filenum, outfile])
        if !test
            cd("./Input/")
            noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
            nomoretowrite = i == length(transforms)
            while (noroomtowrite || nomoretowrite) && length(toread) > 0
                if isfile(toread[1][2])
                    fileinfo = popfirst!(toread)
                    filenum = fileinfo[1]
                    file = fileinfo[2]
                    push!(energies, [filenum, energyfromfile("input$filenum.out")])
                    run(`rm input$filenum.com mp$filenum.pbs input$filenum.out $file`)
                end
            end
        cd("..")
        end
    end
    return energies
end
        
function first()
    """First derivatives"""
    labels, coords = readxyz("geom.xyz")
    rows = length(coords)
    cols = length(coords[1])
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    sampler = zeros(length(coordarray))
    psampler = copy(sampler)
    psampler[1] = 1
    nsampler = copy(sampler)
    nsampler[1] = -1
    transforms = vcat(unique(permutations(psampler)), unique(permutations(nsampler)))
    delta = 0.05
    energies = getenergies(labels, coordarray, delta, transforms, false)
    e0 = popfirst!(energies)[2]
    derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    return derivatives
end

function second()
    """Second derivatives"""
    labels, coords = readxyz("geom.xyz")
    rows = length(coords)
    cols = length(coords[1])
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    sampler = zeros(length(coordarray))
    # adding twice in the same place
    p2sampler = copy(sampler)
    p2sampler[1] = 2
    # subtracting twice in the same place
    n2sampler = copy(sampler)
    n2sampler[1] = -2
    # adding twice in different places 
    ppsampler = copy(sampler)
    ppsampler[1] = 1
    ppsampler[2] = 1
    # subtracting twice in different places 
    nnsampler = copy(sampler)
    nnsampler[1] = -1
    nnsampler[2] = -1
    # adding once, subtracting once
    npsampler = copy(sampler)
    npsampler[1] = -1
    npsampler[2] = 1
    transforms = vcat(unique(permutations(p2sampler)),
                      unique(permutations(n2sampler)),
                      unique(permutations(ppsampler)),
                      unique(permutations(nnsampler)),
                      unique(permutations(npsampler)))
    delta = 0.05
    energies = getenergies(labels, coordarray, delta, transforms, false)
    return energies
    #e0 = popfirst!(energies)[2]
    #derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    #return derivatives
end

function third()
    """Third derivatives"""
    labels, coords = readxyz("geom.xyz")
    rows = length(coords)
    cols = length(coords[1])
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    sampler = zeros(length(coordarray))
    # add thrice in the same place
    p3sampler = copy(sampler)
    p3sampler[1] = 3
    # subtract thrice in the same place
    n3sampler = copy(sampler)
    n3sampler[1] = -3
    # add twice, add once
    p2psampler = copy(sampler)
    p2psampler[1] = 2
    p2psampler[2] = 1
    # add twice, subtract once
    p2nsampler = copy(sampler)
    p2nsampler[1] = 2
    p2nsampler[2] = -1
    # subtract twice, add once
    n2psampler = copy(sampler)
    n2psampler[1] = -2
    n2psampler[2] = 1
    # subtract twice, subtract once
    n2nsampler = copy(sampler)
    n2nsampler[1] = -2
    n2nsampler[2] = -1
    # add thrice in different places
    pppsampler = copy(sampler)
    pppsampler[1] = 1
    pppsampler[2] = 1
    pppsampler[3] = 1
    # add twice, subtrace once different places
    ppnsampler = copy(sampler)
    ppnsampler[1] = 1
    ppnsampler[2] = 1
    ppnsampler[3] = -1
    # add once, subtrace twice different places
    pnnsampler = copy(sampler)
    pnnsampler[1] = 1
    pnnsampler[2] = -1
    pnnsampler[3] = -1
    # subtract thrice in different places
    nnnsampler = copy(sampler)
    nnnsampler[1] = -1
    nnnsampler[2] = -1
    nnnsampler[3] = -1
    transforms = vcat(unique(permutations(p3sampler)),
                      unique(permutations(n3sampler)),
                      unique(permutations(p2psampler)),
                      unique(permutations(p2nsampler)),
                      unique(permutations(n2psampler)),
                      unique(permutations(n2nsampler)),
                      unique(permutations(pppsampler)),
                      unique(permutations(ppnsampler)),
                      unique(permutations(pnnsampler)),
                      unique(permutations(nnnsampler)))
    delta = 0.05
    energies = getenergies(labels, coordarray, delta, transforms, false)
    return energies
    #e0 = popfirst!(energies)[2]
    #derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    #return derivatives
end

println(second())
#println(third())
