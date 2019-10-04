include("../Utilities.jl")
using Combinatorics

MAXFILES = 10

macro unperm(arr)
    return :( unique(permutations($(esc(arr)))) )
end

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
            j = i
            files = []
            while j < i+2
                filenum = makefilenum(transforms[j])
                coords = coordarray + delta*reshape(transforms[j], (length(labels), 3))
                outfile = submitjob(makejob(filenum, labels, coords, 312, test), filenum)
                println(outfile, "\n")
                j += 1
                push!(files, [filenum, outfile])
            end
            push!(toread, files)
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
    firsttransforms = zipper(unique(permutations(psampler)), unique(permutations(nsampler)))
    mixsampler = copy(sampler)
    mixsampler[1] = 1
    mixsampler[2] = -1
    x = unique(permutations(mixsampler))
    i = 1
    j = length(x[1]) - 1
    k = length(x[1]) - 1
    pn = []
    while i <= length(x)
        push!(pn, (x[i:j]...))
        i += k*2
        j += k*2 - 1
        k -= 1
    end
    np = -copy(pn)
    pp = copy(pn)
    nn = copy(pn)
    for i in 1:length(pp)
        pp[i] = map(abs, pp[i])
        nn[i] = -map(abs, nn[i])
    end
    zip1 = zipper(pp, np)
    zip2 = zipper(pn, nn)
    mixed = zipper(zip1, zip2)
    transforms = vcat(firsttransforms, mixed)
    numfiles = 0
    toread = []
    refjob = makejob(0, labels, coordarray, 312, test)
    refout = submitjob(refjob, 0)
    println(refout, "\n")
    push!(toread, [0, refout])
    numfiles += 1
    # flag for checking if e0 has been read yet
    unread = true
    d2s = []
    e0 = nothing
    i = 1
    while i <= length(transforms)
        # change for second
        j = i
        if 2 in transforms[i]
            files = []
            while j < i+2
                filenum = makefilenum(transforms[j])
                coords = coordarray + delta*reshape(transforms[j], (length(labels), 3))
                outfile = submitjob(makejob(filenum, labels, coords, 312, test), filenum)
                println(outfile, "\n")
                j += 1
                push!(files, [filenum, outfile])
            end
            push!(toread, files)
            numfiles += 2
            i += 2
        else
            files = []
            while j < i+4
                filenum = makefilenum(transforms[j])
                coords = coordarray + delta*reshape(transforms[j], (length(labels), 3))
                outfile = submitjob(makejob(filenum, labels, coords, 312, test), filenum)
                println(outfile, "\n")
                j += 1
                push!(files, [filenum, outfile])
            end
            push!(toread, files)
            numfiles += 4
            i += 4
        end

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
                    if isfile(toread[1][1][2]) && isfile(toread[1][2][2]) &&
                        isfile(toread[1][3][2]) && isfile(toread[1][4][2])

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
                        run(`rm input$ppfilenum.com mp$ppfilenum.pbs input$ppfilenum.out $ppfile`)
                        run(`rm input$pnfilenum.com mp$pnfilenum.pbs input$pnfilenum.out $pnfile`)
                        run(`rm input$npfilenum.com mp$npfilenum.pbs input$npfilenum.out $npfile`)
                        run(`rm input$nnfilenum.com mp$nnfilenum.pbs input$nnfilenum.out $nnfile`)
                        numfiles -=4
                        push!(d2s, (fpp - fpm - fmp + fmm)/(4*delta^2))
                    end
                end
            end
        cd("..")
        end
    end
    return d2s
end

function third(xyzfile, delta)
    """Returns the third derivatives"""
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
    # unmixed third derivatives
    # (E(3delta) - 3E(delta) + 3E(-delta) - E(-3delta)) / (2delta)^3
    psampler = copy(sampler)
    psampler[1] = 1
    nsampler = copy(sampler)
    nsampler[1] = -1
    p3sampler = 3*copy(psampler)
    n3sampler = 3*copy(nsampler)
    unmixed = zipper(@unperm(p3sampler), @unperm(psampler), @unperm(nsampler),
                     @unperm(n3sampler))
    # 2/1 derivatives
    # E(2dx;dy) - 2*E(0;dy) + E(-2dx;dy) - E(2dx;-dy) + 2*E(0;-dy) - E(-2dx;-dy)) / (2*delta)^3
    # add twice, add once
    p2psampler = copy(sampler)
    p2psampler[1] = 2
    p2psampler[2] = 1
    p2pperms = @unperm(p2psampler)
    zpperms = copy(p2pperms)
    # get (0;dy) from (2dx;dy) by setting the 2 to 0
    for i in 1:length(zpperms)
        zpperms[i] = map(x -> x%2, zpperms[i])
    end
    # subtract twice, add once
    n2psampler = copy(sampler)
    n2psampler[1] = -2
    n2psampler[2] = 1
    n2pperms = @unperm(n2psampler)
    # add twice, subtract once
    p2nsampler = copy(sampler)
    p2nsampler[1] = 2
    p2nsampler[2] = -1
    p2nperms = @unperm(p2nsampler)
    znperms = copy(p2nperms)
    # get (0;-dy) from (2dx;-dy) by setting the 2 to 0
    for i in 1:length(znperms)
        znperms[i] = map(x -> x%2, znperms[i])
    end
    # subtract twice, subtract once
    n2nsampler = copy(sampler)
    n2nsampler[1] = -2
    n2nsampler[2] = -1
    n2nperms = @unperm(n2nsampler)

    twomixed = zipper(p2pperms, zpperms, n2pperms, p2nperms, znperms, n2nperms)

    # 1;1;1 derivatives
    # (E(ppp) - E(pmp) - E(mpp) + E(mmp) - E(ppm) + E(pmm) + E(mpm) - E(mmm)) / (2*delta)^3
    # add thrice in different places
    pppsampler = copy(sampler)
    pppsampler[1] = 1
    pppsampler[2] = 1
    pppsampler[3] = 1
    pppperms = @unperm(pppsampler)

    nppperms = negate(pppperms, 1, 1)
    pnpperms = negate(pppperms, 1, 2)
    ppnperms = negate(pppperms, 1, 3)

    nnnsampler = -copy(pppsampler)
    nnnperms = @unperm(nnnsampler)

    pnnperms = negate(nnnperms, -1, 1)
    npnperms = negate(nnnperms, -1, 2)
    nnpperms = negate(nnnperms, -1, 3)

    threemixed = zipper(pppperms, pnpperms, nppperms, nnpperms, ppnperms, pnnperms,
                        npnperms, nnnperms)

    # cat them all together
    transforms = vcat(unmixed, twomixed, threemixed)

    numfiles = 0
    toread = []
    d3s = []
    i = 1
    while i <= length(transforms)
        if 3 in transforms[i]
            # (E(3delta) - 3E(delta) + 3E(-delta)-E(-3delta)) / (2delta)^3
            j = i
            files = []
            while j < i+4
                filenum = makefilenum(transforms[j])
                coords = coordarray + delta*reshape(transforms[j], (length(labels), 3))
                outfile = submitjob(makejob(filenum, labels, coords, 312, test), filenum)
                println(outfile, "\n")
                j += 1
                push!(files, [filenum, outfile])
            end
            push!(toread, files)
            numfiles += 4
            i += 4
        elseif 2 in transforms[i]
            # TODO 2/1 derivatives
            # E(2dx;dy) - 2*E(0;dy) + E(-2dx;dy) - E(2dx;-dy) + 2*E(0;-dy) - E(-2dx;-dy)) / (2*delta)^3
            ppnum = makefilenum(transforms[i])
            coords = coordarray + delta*reshape(transforms[i], (length(labels), 3))
            ppout = submitjob(makejob(ppnum, labels, coords, 312, test), ppnum)
            println(ppout, "\n")

            pnnum = makefilenum(transforms[i+1])
            coords = coordarray + delta*reshape(transforms[i+1], (length(labels), 3))
            pnout = submitjob(makejob(pnnum, labels, coords, 312, test), pnnum)
            println(pnout, "\n")

            npnum = makefilenum(transforms[i+2])
            coords = coordarray + delta*reshape(transforms[i+2], (length(labels), 3))
            npout = submitjob(makejob(npnum, labels, coords, 312, test), npnum)
            println(npout, "\n")

            nnnum = makefilenum(transforms[i+3])
            coords = coordarray + delta*reshape(transforms[i+3], (length(labels), 3))
            nnout = submitjob(makejob(nnnum, labels, coords, 312, test), nnnum)
            println(nnout, "\n")

            push!(toread, ([ppnum, ppout], [pnnum, pnout], [npnum, npout],
                           [nnnum, nnout]))
            numfiles += 4
            i += 4

        else
            # TODO 1;1;1 derivatives
            # (E(ppp) - E(pmp) - E(mpp) + E(mmp) - E(ppm) + E(pmm) + E(mpm) - E(mmm)) / (2*delta)^3
        end

        if !test
            cd("./Input/")
            noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
            nomoretowrite = i == length(transforms) - 1
            while (noroomtowrite || nomoretowrite) && length(toread) > 0
                if length(toread[1]) == 4
                    # TODO unmixed third derivatives
                    # (E(3delta) - 3E(delta) + 3E(-delta)-E(-3delta)) / (2delta)^3
                    # push!(d3s, (ep3d - 3*epd + 3*emd - em3d)/(2*delta)^3)
                    # use first derivatives and multiply by 3 to get the first and last ones
                    if isfile(toread[1][1][2]) && isfile(toread[1][2][2]) &&
                        isfile(toread[1][3][2]) && isfile(toread[1][4][2])

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
                        push!(d2s, (fpp - fpm - fmp + fmm)/(4*delta^2))
                    end
                elseif length(toread[1]) == 6
                    # 2/1 mixed derivatives
                    # E(2dx;dy) - 2*E(0;dy) + E(-2dx;dy) - E(2dx;-dy) + 2*E(0;-dy) - E(-2dx;-dy)) / (2*delta)^3
                    # push!(d3s, (ep2p - 2*e0p + em2p - ep2m + 2*e0m - em2m) / (2*delta)^3)
                elseif length(toread[1]) == 8
                    # 1;1;1 derivatives
                    # (E(ppp) - E(pmp) - E(mpp) + E(mmp) - E(ppm) + E(pmm) + E(mpm) - E(mmm)) / (2*delta)^3
                    # push!(d3s, (eppp - epmp - empp + emmp - eppm + epmm + empm - emmm) / (2*delta)^3)
                end
            end
        cd("..")
        end
    end
    return d3s
end

#println(first("../geom.xyz", 0.005))
#println(second("../geom.xyz", 0.005))
println(third("../geom.xyz", 0.005))
