using Weave
fileset = readdir(@__DIR__)
for file in fileset
    if occursin(".jmd", file) 
        Weave.notebook(file)
        println("Weaved $file")
    end
end
