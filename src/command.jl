
module EmmaCommand

using ArgMacros, Emma

args = @dictarguments begin
    @helpusage "Emma/src/command.jl [options] <FASTA_file>"
    @helpdescription """
        Note: Use consistant inputs/outputs. If you wish
        to annotate a directory of fasta files, ensure that
        the output parameters are also directories.
        """
    @argumentdefault Int16 2 translation_table "--transl_table" 
    @arghelp "NCBI translation table; 2 for vertebrates (the default), 5 for invertebrates"
    @argumentoptional String rotate_to "--rotate" 
    @arghelp "rotate genome and annotations to start with this feature"
    @argumentoptional String FA_out "--fa"
    @arghelp "file/dir for fasta output"
    @argumentoptional String GFF_out "--gff" 
    @arghelp "file/dir for gff output"
    @argumentoptional String GB_out "--tbl"
    @arghelp "file/dir for .tbl output (for GenBank submissions)"
    @argumentoptional String SVG_out "--svg"
    @arghelp "file/dir for svg output"
    @positionalrequired String FASTA_file
    @arghelp "file/dir for fasta input"
end
println(ARGS)
file_args = [args[:FA_out], args[:GFF_out], args[:GB_out], args[:SVG_out], args[:FASTA_file]]
filtered_args = filter(x -> ~isnothing(x), file_args)
all_dirs = all(isdir, filtered_args)
all_files = !any(isdir, filtered_args)
if all_dirs
    fafiles = filter!(x->endswith(x,".fa") || endswith(x,".fasta"), readdir(args[:FASTA_file], join=true))
    for fasta in fafiles
        accession = first(split(basename(fasta),"."))
        outfile_gff = haskey(filtered_args, :GFF_out) ? joinpath(args[:GFF_out], accession * ".gff") : nothing
        outfile_gb = haskey(filtered_args, :GB_out) ? joinpath(args[:GB_out], accession * ".tbl") : nothing
        outfile_fa = haskey(filtered_args, :FA_out) ? joinpath(args[:FA_out], accession * ".fa") : nothing
        outfile_svg = haskey(filtered_args, :SVG_out) ? joinpath(args[:SVG_out], accession * ".svg") : nothing
        emma(fasta; translation_table = args[:translation_table], rotate_to = args[:rotate_to], outfile_gff=outfile_gff, outfile_gb=outfile_gb, outfile_fa=outfile_fa, outfile_svg=outfile_svg)
    end
elseif all_files
    emma(args[:FASTA_file]; translation_table = args[:translation_table], rotate_to = args[:rotate_to], outfile_gff = args[:GFF_out], outfile_gb = args[:GB_out], outfile_fa = args[:FA_out], outfile_svg = args[:SVG_out])
else
    throw("Inputs must be consistant; all directories or all files")
end

end