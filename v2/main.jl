include("../Utilities.jl")
using Combinatorics

MAXFILES = 10

function first(xyzfile, delta)
    """Returns the first derivatives"""
    cd("/ddn/home1/r2518/research/Force-Constants/v2")
    test = false
    try 
        rm("./Input", recursive=true)
    catch e
    end
    numfiles = 0
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
        numfiles += 2
        if !test
            cd("./Input/")
            noroomtowrite = numfiles > MAXFILES
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
                    noroomtowrite = numfiles > MAXFILES
                    run(`rm input$pfilenum.com mp$pfilenum.pbs input$pfilenum.out $pfile`)
                    run(`rm input$nfilenum.com mp$nfilenum.pbs input$nfilenum.out $nfile`)
                    numfiles -= 2
                end
            end
        cd("..")
        end
    end
    return d1s
end

function second(xyzfile, delta)
    """Returns the second derivatives"""
    cd("/ddn/home1/r2518/research/Force-Constants/v2")
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
    psampler[1] = 2
    nsampler = copy(sampler)
    nsampler[1] = -2
    #firsttransforms = zipper(unique(permutations(psampler)), unique(permutations(nsampler)))
    transforms = zipper(unique(permutations(psampler)), unique(permutations(nsampler)))
    # TODO remember pp and mm, not just the mixed signs and put them in the right places
    #      in transforms too
    #      Should be set up below, just need to zip them together properly
    mixsampler = copy(sampler)
    mixsampler[1] = 1
    mixsampler[2] = -1
    x = unique(permutations(mixsampler))
    i = 1
    j = 8
    k = 8
    mix = []
    while i <= length(x)
        push!(mix, (x[i:j]...))
        i += k*2
        j += k*2 - 1
        k -= 1
    end
    mixed = zipper(mix, -mix)
    # making the double positives and double negatives
    pp = copy(mix)
    nn = copy(mix)
    for i in 1:length(pp)
        pp[i] = map(abs, pp[i])
        nn[i] = -map(abs, nn[i])
    end
    #transforms = vcat(firsttransforms, mixed)
    numfiles = 0
    toread = []
    refjob = makejob(0, labels, coordarray, 50, test)
    refout = submitjob(refjob, 0)
    println(refout, "\n")
    push!(toread, [0, refout])
    numfiles += 1
    # flag for checking if e0 has been read yet
    unread = true
    d2s = []
    e0 = nothing
    for i in 1:2:length(transforms)
        # change for second
        if 2 in transforms[i]
            filenum = makefilenum(transforms[i])
            coords = coordarray + delta*reshape(transforms[i], (length(labels), 3))
            outfile = submitjob(makejob(filenum, labels, coords, 50, test), filenum)
            println(outfile, "\n")
            filenum2 = makefilenum(transforms[i+1])
            coords = coordarray + delta*reshape(transforms[i+1], (length(labels), 3))
            outfile2 = submitjob(makejob(filenum2, labels, coords, 50, test), filenum2)
            println(outfile2, "\n")
            push!(toread, ([filenum, outfile], [filenum2, outfile2]))
            numfiles += 2
        end
        # TODO add in the mixed derivative file submissions

        if !test
            cd("./Input/")
            noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
            nomoretowrite = i == length(transforms) - 1
            while (noroomtowrite || nomoretowrite) && length(toread) > 0
                # e0 calculation 
                if unread && isfile(toread[1][2])
                    fileinfo = popfirst!(toread)
                    filenum = fileinfo[1]
                    file = fileinfo[2]
                    e0 = energyfromfile("input$filenum.out")
                    noroomtowrite = numfiles > MAXFILES
                    run(`rm input$filenum.com mp$filenum.pbs input$filenum.out $file`)
                    unread = false
                    numfiles -= 1
                # unmixed derivatives
                elseif !unread && length(toread[1]) == 2
                    if isfile(toread[1][1][2]) && isfile(toread[1][2][2])
                        files = popfirst!(toread)
                        pfilenum = files[1][1]
                        pfile = files[1][2]
                        nfilenum = files[2][1]
                        nfile = files[2][2]
                        fxp2x = energyfromfile("input$pfilenum.out")
                        fxm2x = energyfromfile("input$nfilenum.out")
                        push!(d2s, (fxp2x - 2*e0 + fxm2x)/(2*delta)^2)
                        run(`rm input$pfilenum.com mp$pfilenum.pbs input$pfilenum.out $pfile`)
                        run(`rm input$nfilenum.com mp$nfilenum.pbs input$nfilenum.out $nfile`)
                        numfiles -= 2
                    end
                # mixed derivatives
                elseif !unread && length(toread[1]) == 4
                    if isfile(toread[1][1][2]) && isfile(toread[1][2][2]) && isfile(toread[1][2][3]) && isfile(toread[1][2][4])
                        files = popfirst!(toread)
                        ppfilenum = files[1][1]
                        ppfile = files[1][2]
                        pnfilenum = files[2][1]
                        pnfile = files[2][2]
                        npfilenum = files[3][1]
                        npfile = files[3][2]
                        nnfilenum = files[4][1]
                        nnfile = files[4][2]
                        fpp = energyfromfile("input$ppfilenum.out")
                        fpm = energyfromfile("input$pnfilenum.out")
                        fmp = energyfromfile("input$npfilenum.out")
                        fmm = energyfromfile("input$nnfilenum.out")
                        numfiles -=4
                    end
                    #push!(d2s, (fpp - fpm - fmp + fmm)/(4*delta^2))
                end
            end
        cd("..")
        end
    end
    return d2s
end

#println(first("../geom.xyz", 0.005))
println(second("../geom.xyz", 0.005))
