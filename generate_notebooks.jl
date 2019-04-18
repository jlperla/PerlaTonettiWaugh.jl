using Weave
fileset = readdir(@__DIR__)
for file in fileset
    if occursin(".jmd", file) && occursin("compare_steady_states", file) == false # disable this pending recalibration
        Weave.notebook(file)
        println("Weaved $file")
    end
end
