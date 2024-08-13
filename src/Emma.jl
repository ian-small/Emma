module Emma

using Serialization
using Artifacts
using BioSequences
using FASTX
using Logging
using Unicode
using GenomicAnnotations
using UUIDs

export main, emma, emmaone, writeGFF, writeGB, tempfilename, TempFile, drawgenome, rotate

#const emmamodels = "/data/Emma/emma_vertebrate_models"
const emmamodels = joinpath(artifact"Emma_vertebrate_models", "emma-models-1.0.0")

include("tempfile.jl")
include("circularity.jl")
include("feature.jl")
include("orfs.jl")
include("tRNAs.jl")
include("rRNAs.jl")
include("gff.jl")
include("gb.jl")
include("visuals.jl")
include("process.jl")
include("cmd.jl")

end
