
struct TempFile
    directory::String
    uuid::UUID
    ext::Vector{String}
    function TempFile(directory::String=".")
        new(directory, uuid1(), [])
    end
end
function tempfilename(tempfile::TempFile, ext::String)
    push!(tempfile.ext, ext)
    joinpath(tempfile.directory, "$(tempfile.uuid).$(ext)")
end

function cleanfiles(tempfile::TempFile)
    while length(tempfile.ext) > 0
        ext = pop!(tempfile.ext)
        # don't add to tempfile.ext with tempfilename!!
        path = joinpath(tempfile.directory, "$(tempfile.uuid).$(ext)")
        rm(path, force=true)
    end

end

function which(cmd::String)::String
    # place holder... do we have cmd in our PATH or not...
    # this should be like python which.
    cmd
end
