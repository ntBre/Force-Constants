include("./Utilities.jl")
using Combinatorics

MAXFILES = 250

function firstderivative(labels, coords, delta, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    rows = length(coords)
    cols = length(coords[1])
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    rm("./Input", recursive=true)
    mkdir("Input")
    refjob = makejob(0, labels, coordarray, 50, test)
    refout = submitjob(refjob, 0)
    toread = []
    push!(toread, [0, refout])
    energies = []
    sampler = zeros(length(coordarray))
    sampler[1] = 1
    transforms = unique(permutations(sampler))
    for i in 1:length(transforms)
        posfilenum = i
        negfilenum = -i
        pcoords = coordarray + delta*reshape(transforms[i], (rows, cols))
        ncoords = coordarray - delta*reshape(transforms[i], (rows, cols))
        posout = submitjob(makejob(posfilenum, labels, pcoords, 50, test), posfilenum)
        negout = submitjob(makejob(negfilenum, labels, ncoords, 50, test), negfilenum)
        println(posout, "\n", negout, "\n\n")
        push!(toread, [posfilenum, posout], [negfilenum, negout])
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
    e0 = popfirst!(energies)[2]
    derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    return derivatives
end
        
labels, coords = readxyz("geom.xyz")
println(firstderivative(labels, coords, 0.05, false))
