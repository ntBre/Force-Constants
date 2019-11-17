using DelimitedFiles
using Printf

data = readdlm("practice.dat")
data = transpose(data)
#println(display(readdlm("practice.dat")))
for i in 1:3:length(data)-2
    s = map(x -> @sprintf("%20.10f", x), data[i:i+2])
    println(join(s))
end

