include("../Utilities.jl")
using Combinatorics

MAXFILES = 10

function first(xyzfile, delta)
    """Returns the first derivatives"""
    test = false
    try 
        rm("./Input", recursive=true)
    catch e
    end
    mkdir("Input")
    labels, coords = readxyz(xyzfile)
    rows = length(coords)
    cols = length(coords[1])
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    sampler = zeros(length(coordarray))
    psampler = copy(sampler)
    psampler[1] = 1
    nsampler = copy(sampler)
    nsampler[1] = -1
    transforms = zipper(unique(permutations(psampler)), unique(permutations(nsampler)))
    toread = []
    d1s = []
    for i in 1:2:length(transforms)
        filenum = makefilenum(transforms[i])
        coords = coordarray + delta*reshape(transforms[i], (length(labels), 3))
        outfile = submitjob(makejob(filenum, labels, coords, 50, test), filenum)
        println(outfile, "\n")
        filenum2 = makefilenum(transforms[i+1])
        coords = coordarray + delta*reshape(transforms[i+1], (length(labels), 3))
        outfile2 = submitjob(makejob(filenum2, labels, coords, 50, test), filenum2)
        println(outfile2, "\n")
        push!(toread, ([filenum, outfile], [filenum2, outfile2]))
        if !test
            cd("./Input/")
            noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
            nomoretowrite = i == length(transforms) - 1
            while (noroomtowrite || nomoretowrite) && length(toread) > 0
                if isfile(toread[1][1][2]) && isfile(toread[1][2][2])
                    files = popfirst!(toread)
                    pfilenum = files[1][1]
                    pfile = files[1][2]
                    nfilenum = files[2][1]
                    nfile = files[2][2]
                    fxpx = energyfromfile("input$pfilenum.out")
                    fxmx = energyfromfile("input$nfilenum.out")
                    push!(d1s, (fxpx - fxmx) / (2 * delta))
                    noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
                    run(`rm input$pfilenum.com mp$pfilenum.pbs input$pfilenum.out $pfile`)
                    run(`rm input$nfilenum.com mp$nfilenum.pbs input$nfilenum.out $nfile`)
                end
            end
        cd("..")
        end
    end
    return d1s
end

function second(xyzfile, delta)
    """Returns the second derivatives"""
    test = false
    try 
        rm("./Input", recursive=true)
    catch e
    end
    mkdir("Input")
    labels, coords = readxyz(xyzfile)
    rows = length(coords)
    cols = length(coords[1])
    coordarray = reshape(vcat(copy(coords)...), (rows, cols))
    sampler = zeros(length(coordarray))
    # change for second
    psampler = copy(sampler)
    psampler[1] = 1
    nsampler = copy(sampler)
    nsampler[1] = -1
    transforms = zipper(unique(permutations(psampler)), unique(permutations(nsampler)))

    toread = []
    refjob = makejob(0, labels, coordarray, 50, test)
    refout = submitjob(refjob, 0)
    println(refout, "\n")
    push!(toread, [0, refout])
    # flag for checking if e0 has been read yet
    unread = true
    d2s = []
    for i in 1:2:length(transforms)
        # change for second
        filenum = makefilenum(transforms[i])
        coords = coordarray + delta*reshape(transforms[i], (length(labels), 3))
        outfile = submitjob(makejob(filenum, labels, coords, 50, test), filenum)
        println(outfile, "\n")
        filenum2 = makefilenum(transforms[i+1])
        coords = coordarray + delta*reshape(transforms[i+1], (length(labels), 3))
        outfile2 = submitjob(makejob(filenum2, labels, coords, 50, test), filenum2)
        println(outfile2, "\n")
        push!(toread, ([filenum, outfile], [filenum2, outfile2]))

        if !test
            cd("./Input/")
            noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
            nomoretowrite = i == length(transforms) - 1
            while (noroomtowrite || nomoretowrite) && length(toread) > 0
                # change for second
                if unread && isfile(toread[1][2])
                    fileinfo = popfirst!(toread)
                    filenum = fileinfo[1]
                    file = fileinfo[2]
                    e0 = energyfromfile("input$filenum.out")
                    noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
                    run(`rm input$filenum.com mp$filenum.pbs input$filenum.out $file`)
                    unread = false
                elseif isfile(toread[1][1][2]) && isfile(toread[1][2][2])
                    files = popfirst!(toread)
                    pfilenum = files[1][1]
                    pfile = files[1][2]
                    nfilenum = files[2][1]
                    nfile = files[2][2]
                    fxpx = energyfromfile("input$pfilenum.out")
                    fxmx = energyfromfile("input$nfilenum.out")
                    #if +2delta and -2delta pair
                    #push!(d2s, (fxp2x - 2*f0 + fxm2x)/(2*delta)^2)
                    #else if mixed
                    #push!(d2s, (fpp - fpm - fmp + fmm)/(4*delta^2))
                    noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
                    run(`rm input$pfilenum.com mp$pfilenum.pbs input$pfilenum.out $pfile`)
                    run(`rm input$nfilenum.com mp$nfilenum.pbs input$nfilenum.out $nfile`)
                end
            end
        cd("..")
        end
    end
    return d2s
end

#println(first("../geom.xyz", 0.005))
println(second("../geom.xyz", 0.005))
    
