include("./Utilities.jl")
using Combinatorics

MAXFILES = 250

function getenergies(labels, coords, delta, transforms, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    rows = length(coords)
    cols = length(coords[1])
    # vcat the splatted coords array and reshape it
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    rm("./Input", recursive=true)
    mkdir("Input")
    refjob = makejob(0, labels, coordarray, 50, test)
    refout = submitjob(refjob, 0)
    toread = []
    push!(toread, [0, refout])
    energies = []
    for i in 1:length(transforms)
        filenum = makefilenum(transforms[i])
        coords = coordarray + delta*reshape(transforms[i], (rows, cols))
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
    sampler = zeros(length(coordarray))
    psampler = copy(sampler)
    psampler[1] = 1
    nsampler = copy(sampler)
    nsampler[1] = -1
    transforms = vcat(unique(permutations(psampler)), unique(permutations(nsampler)))

    labels, coords = readxyz("geom.xyz")
    delta = 0.05
    energies = getenergies(labels, coords, delta, transforms, false)

    e0 = popfirst!(energies)[2]
    derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    return derivatives
end

println(first())
