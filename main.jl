include("./Utilities.jl")

# DONE check if there is space to write more files
# DONE check if output file is finished
# DONE capture energy from output
# TODO perform calculation on energy and save
# TODO refactor to submit one job at a time and remove pos and neg prefixes
# 	This should make it easier to extend to higher dimensions

MAXFILES = 250

function firstderivative(labels, coords, delta, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    tempcoords = copy(coords)
    rm("./Input", recursive=true)
    mkdir("Input")
    energies = []
    toread = []
    for atom in 1:length(coords)
        for coord in 1:length(coords[atom])
            # makenames
            posfilenum = parse(Int, string(atom) * string(coord))
            negfilenum = -parse(Int, string(atom) * string(coord))
            poscomname = "input$posfilenum.com"
            negcomname = "input$negfilenum.com"
            pospbsname = "mp$posfilenum.pbs"
            negpbsname = "mp$negfilenum.pbs"
            # writefiles
            tempcoords[atom][coord] += delta
            writecom(poscomname, labels, tempcoords, 50, test)
            writepbs(posfilenum, test)
            tempcoords[atom][coord] -= 2*delta
            writecom(negcomname, labels, tempcoords, 50, test)
            writepbs(negfilenum, test)
            # submit jobs
            run(`mv $poscomname ./Input/.`)
            run(`mv $pospbsname ./Input/.`)
            run(`mv $negcomname ./Input/.`)
            run(`mv $negpbsname ./Input/.`)
            cd("./Input")
            posjobnum = read(`qsub $pospbsname`, String)
            negjobnum = read(`qsub $negpbsname`, String)
            println(posjobnum, negjobnum)
            posout = "job$posfilenum.o$(posjobnum[1:5])"
            negout = "job$negfilenum.o$(negjobnum[1:5])"
            push!(toread, [posfilenum, posout], [negfilenum, negout])
            if !test
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
    return energies
end

labels, coords = readxyz("geom.xyz")
println(firstderivative(labels, coords, 0.05, false))

