include("./Utilities.jl")

# DONE check if there is space to write more files
# DONE check if output file is finished
# DONE capture energy from output
# DONE perform calculation on energy and save
# TODO refactor to submit one job at a time and remove pos and neg prefixes
# 	This should make it easier to extend to higher dimensions
#       Mostly done

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
            end
            cd("..")
            tempcoords[atom][coord] += delta
        end
    end
    e0 = popfirst!(energies)[2]
    derivatives = map(x -> [x[1], x[2] - e0 / delta], energies)
    return derivatives
end

labels, coords = readxyz("geom.xyz")
println(firstderivative(labels, coords, 0.05, false))

