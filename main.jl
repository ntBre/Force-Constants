include("./Utilities.jl")

MAXFILES = 250

function firstderivative(labels, coords, delta, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    tempcoords = copy(coords)
    rm("./Input", recursive=true)
    mkdir("Input")
    refjob = makejob(0, labels, tempcoords, 50, test)
    refout = submitjob(refjob, 0)
    toread = []
    push!(toread, [0, refout])
    energies = []
    for atom in 1:length(coords)
        for coord in 1:length(coords[atom])
            posfilenum = parse(Int, string(atom) * string(coord))
            negfilenum = -parse(Int, string(atom) * string(coord))
            tempcoords[atom][coord] += delta
            posout = submitjob(makejob(posfilenum, labels, tempcoords, 50, test), posfilenum)
            tempcoords[atom][coord] -= 2*delta
            negout = submitjob(makejob(negfilenum, labels, tempcoords, 50, test), negfilenum)
            println(posout, "\n", negout, "\n\n")
            push!(toread, [posfilenum, posout], [negfilenum, negout])
            if !test
                cd("./Input/")
                noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
                nomoretowrite = atom == length(coords) & coord == length(coords[atom])
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
            tempcoords[atom][coord] += delta
        end
    end
    e0 = popfirst!(energies)[2]
    derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    return derivatives
end

# TODO calculate the right derivatives
function secondderivative(labels, coords, delta, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the second derivative energy values"""
    tempcoords = copy(coords)
    rm("./Input", recursive=true)
    mkdir("./Input")
    refjob = makejob(0, labels, tempcoords, 50, test)
    refout = submitjob(refjob, 0)
    toread = []
    push!(toread, [0, refout])
    energies = []
    for atom in 1:length(coords)
        for coord in 1:length(coords[atom])
            for repeatcoord in 1:length(coords[atom])
                ppfilenum = string(atom) * string(coord) * string(atom) * string(repeatcoord)
                nnfilenum = "-" * string(atom) * string(coord) * "-" * string(atom) * string(repeatcoord)
                tempcoords[atom][coord] += delta
                tempcoords[atom][repeatcoord] += delta
                ppout = submitjob(makejob(ppfilenum, labels, tempcoords, 50, test), ppfilenum)
                tempcoords[atom][coord] -= 2*delta
                tempcoords[atom][repeatcoord] -= 2*delta
                nnout = submitjob(makejob(nnfilenum, labels, tempcoords, 50, test), nnfilenum)
                println(ppout, "\n", nnout, "\n\n")
                push!(toread, [ppfilenum, ppout], [nnfilenum, nnout])
                if coord != repeatcoord
                    pnfilenum = string(atom) * string(coord) * "-" * string(atom) * string(repeatcoord)
                    npfilenum = "-" * string(atom) * string(coord) * string(atom) * string(repeatcoord)
                    tempcoords[atom][coord] += 2*delta
                    pnout = submitjob(makejob(pnfilenum, labels, tempcoords, 50, test), pnfilenum)
                    tempcoords[atom][coord] -= 2*delta
                    tempcoords[atom][repeatcoord] += 2*delta
                    npout = submitjob(makejob(npfilenum, labels, tempcoords, 50, test), npfilenum)
                    println(pnout, "\n", npout, "\n\n")
                    push!(toread, [pnfilenum, pnout], [npfilenum, npout])
                end
                if !test
                    cd("./Input/")
                    noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
                    nomoretowrite = atom == length(coords) & coord == length(coords[atom])
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
            tempcoords = copy(coords)
        end
    end
    e0 = popfirst!(energies)[2]
    derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    return derivatives
end

function thirdderivative(labels, coords, delta, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the third derivative energy values"""
    tempcoords = copy(coords)
    rm("./Input", recursive=true)
    mkdir("./Input")
    refjob = makejob(0, labels, tempcoords, 50, test)
    refout = submitjob(refjob, 0)
    toread = []
    push!(toread, [0, refout])
    energies = []
    for atom in 1:length(coords)
        for coord in 1:length(coords[atom])
            for repeatcoord in 1:length(coords[atom])
                pppfilenum = string(atom) * string(coord) * string(atom) * string(repeatcoord)
                nnnfilenum = "-" * string(atom) * string(coord) * "-" * string(atom) * string(repeatcoord)
                tempcoords[atom][coord] += delta
                tempcoords[atom][repeatcoord] += delta
                ppout = submitjob(makejob(ppfilenum, labels, tempcoords, 50, test), ppfilenum)
                tempcoords[atom][coord] -= 2*delta
                tempcoords[atom][repeatcoord] -= 2*delta
                nnout = submitjob(makejob(nnfilenum, labels, tempcoords, 50, test), nnfilenum)
                println(ppout, "\n", nnout, "\n\n")
                push!(toread, [ppfilenum, ppout], [nnfilenum, nnout])
                if coord != repeatcoord
                    pnfilenum = string(atom) * string(coord) * "-" * string(atom) * string(repeatcoord)
                    npfilenum = "-" * string(atom) * string(coord) * string(atom) * string(repeatcoord)
                    tempcoords[atom][coord] += 2*delta
                    pnout = submitjob(makejob(pnfilenum, labels, tempcoords, 50, test), pnfilenum)
                    tempcoords[atom][coord] -= 2*delta
                    tempcoords[atom][repeatcoord] += 2*delta
                    npout = submitjob(makejob(npfilenum, labels, tempcoords, 50, test), npfilenum)
                    println(pnout, "\n", npout, "\n\n")
                    push!(toread, [pnfilenum, pnout], [npfilenum, npout])
                end
                if !test
                    cd("./Input/")
                    noroomtowrite = parse(Int64, read(pipeline(`ls`, `grep pbs`, `wc -l`), String)) > MAXFILES
                    nomoretowrite = atom == length(coords) & coord == length(coords[atom])
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
            tempcoords = copy(coords)
        end
    end
    e0 = popfirst!(energies)[2]
    derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    return derivatives
end

labels, coords = readxyz("geom.xyz")
println(firstderivative(labels, coords, 0.05, false))
#println(secondderivative(labels, coords, 0.05, false))
