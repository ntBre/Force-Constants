labels, coords = readxyz("../geom.xyz")
rows = length(coords)
cols = length(coords[1])
coordarray = reshape(vcat(copy(coords)...), (rows, cols))
sampler = zeros(length(coordarray))
# change for second
psampler = copy(sampler)
psampler[1] = 2
nsampler = copy(sampler)
nsampler[1] = -2
transforms = zipper(unique(permutations(psampler)), unique(permutations(nsampler)))
mixsampler = copy(sampler)
mixsampler[1] = 1
mixsampler[2] = -1
# need to do some kind of processing to get the right ones next to each other
x = unique(permutations(mixsampler))
i = 1
j = 8
k = 8
mix = []
while i <= length(x)
    push!(mix, (x[i:j]...))
    global i += k*2
    global j += k*2 - 1
    global k -= 1
end
mixed = zipper(mix, -mix)
for line in 1:length(mixed)
    println(mixed[line])
end
