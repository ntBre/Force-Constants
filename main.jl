include("./Utilities.jl")
using Combinatorics
using DelimitedFiles

MAXFILES = 10

function getenergies(labels, coordarray, delta, transforms, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    rm("./Input", recursive=true)
    mkdir("Input")
    #refjob = makejob(0, labels, coordarray, 50, test)
    #refout = submitjob(refjob, 0)
    toread = []
    #push!(toread, [0, refout])
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
                    noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
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
    return energies
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
end

function fourth()
    """Fourth derivatives"""
    labels, coords = readxyz("geom.xyz")
    rows = length(coords)
    cols = length(coords[1])
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    sampler = zeros(length(coordarray))
    # add 4 once
    p4sampler = copy(sampler)
    p4sampler[1] = 4
    # subtract 4 once
    n4sampler = copy(sampler)
    n4sampler[1] = -4
    # add thrice, add once
    p3psampler = copy(sampler)
    p3psampler[1] = 3
    p3psampler[2] = 1
    # add thrice, subtract once
    p3nsampler = copy(sampler)
    p3nsampler[1] = 3
    p3nsampler[2] = -1
    # subtract thrice, add once
    n3psampler = copy(sampler)
    n3psampler[1] = -3
    n3psampler[2] = 1
    # subtract thrice, subtract once
    n3nsampler = copy(sampler)
    n3nsampler[1] = -3
    n3nsampler[2] = -1
    # add twice, add twice
    p2p2sampler = copy(sampler)
    p2p2sampler[1] = 2
    p2p2sampler[2] = 2
    # add twice, subtract twice
    p2n2sampler = copy(sampler)
    p2n2sampler[1] = 2
    p2n2sampler[2] = -2
    # subtract twice, subtract twice
    n2n2sampler = copy(sampler)
    n2n2sampler[1] = -2
    n2n2sampler[2] = -2
    # add four times in different places
    ppppsampler = copy(sampler)
    ppppsampler[1] = 1
    ppppsampler[2] = 1
    ppppsampler[3] = 1
    ppppsampler[4] = 1
    # add three times, subtract once in different places
    pppnsampler = copy(sampler)
    pppnsampler[1] = 1
    pppnsampler[2] = 1
    pppnsampler[3] = 1
    pppnsampler[4] = -1
    # add two times, subtract twice in different places
    ppnnsampler = copy(sampler)
    ppnnsampler[1] = 1
    ppnnsampler[2] = 1
    ppnnsampler[3] = -1
    ppnnsampler[4] = -1
    # add once, subtract thrice in different places
    pnnnsampler = copy(sampler)
    pnnnsampler[1] = 1
    pnnnsampler[2] = -1
    pnnnsampler[3] = -1
    pnnnsampler[4] = -1
    # subtract four times in different places
    nnnnsampler = copy(sampler)
    nnnnsampler[1] = -1
    nnnnsampler[2] = -1
    nnnnsampler[3] = -1
    nnnnsampler[4] = -1
    transforms = vcat(unique(permutations(p4sampler)),
                      unique(permutations(n4sampler)),
                      unique(permutations(p3psampler)),
                      unique(permutations(p3nsampler)),
                      unique(permutations(n3psampler)),
                      unique(permutations(n3nsampler)),
                      unique(permutations(p2p2sampler)),
                      unique(permutations(p2n2sampler)),
                      unique(permutations(n2n2sampler)),
                      unique(permutations(ppppsampler)),
                      unique(permutations(pppnsampler)),
                      unique(permutations(ppnnsampler)),
                      unique(permutations(pnnnsampler)),
                      unique(permutations(nnnnsampler)))
    delta = 0.05
    energies = getenergies(labels, coordarray, delta, transforms, false)
    return energies
end

function firstderivative()
    energies = readdlm("1st.dat")
    for i in 1:length(energies[1:end, 1])
        for j in i:length(energies[1:end, 1])
            if energies[i] == -energies[j]
                println(energies[i,2] - energies[j,2])
            end
        end
    end
end

firstderivative()

#writedlm("1st.dat", first())
#writedlm("2nd.dat", second())
#writedlm("3rd.dat", third())
#writedlm("4th.dat", fourth())
#println(first())
#println(second())
#println(third())
#println(fourth())
