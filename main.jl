include("./Utilities.jl")

# TODO check if there is space to write more files
# TODO check if output file is finished
# TODO capture energy from output
# TODO perform calculation on energy and save

function firstderivative(labels, coords, delta, test=true)
    """Takes cartesian atomic coordinates and a step size,
        writes a .pbs and .com file for each possible step, 
        runs the necessary calculations in Molpro, and returns
        an array of the first derivative energy values"""
    tempcoords = copy(coords)
    rm("./Input", recursive=true)
    mkdir("Input")
    posenergies = []
    negenergies = []
    for atom in 1:length(coords)
        for coord in 1:length(coords[atom])
            filenum = parse(Int, string(atom) * string(coord))
            comfilename = "input$(filenum).com"
            tempcoords[atom][coord] += delta
            writecom(comfilename, labels, tempcoords, 50, test)
            writepbs(filenum, test)
            run(`mv $(comfilename) ./Input/.`)
            run(`mv mp$(filenum).pbs ./Input/.`)
            tempcoords[atom][coord] -= 2*delta
            comfilename = "input-$(filenum).com"
            writecom(comfilename, labels, tempcoords, 50, test)
            writepbs(-filenum, test)
            run(`mv $(comfilename) ./Input/.`)
            run(`mv mp-$(filenum).pbs ./Input/.`)
            if !test
		cd("./Input")
		posjobnum = read(`qsub mp$(filenum).pbs`, String)
		negjobnum = read(`qsub mp-$(filenum).pbs`, String)
		println(posjobnum, negjobnum)
		posout = "job$filenum.o$(posjobnum[1:5])"
		negout = "job-$filenum.o$(negjobnum[1:5])"
		while !(isfile(posout) & isfile(negout))
		    sleep(1)
		end
		push!(posenergies, energyfromfile("input$filenum.out"))
		push!(negenergies, energyfromfile("input-$filenum.out"))
		run(`rm input$filenum.com input-$filenum.com mp$filenum.pbs mp-$filenum.pbs`)
		run(`rm $posout $negout input$filenum.out input-$filenum.out`)
		cd("..")
            end
            tempcoords[atom][coord] += delta
        end
    end
    return posenergies, negenergies
end

labels, coords = readxyz("geom.xyz")
println(firstderivative(labels, coords, 0.05, false))

