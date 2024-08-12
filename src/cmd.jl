import ArgParse: ArgParseSettings, @add_arg_table!, parse_args
import Logging

MayBeString = Union{Nothing,String}

function mainone(infile::String; gff::MayBeString=nothing, gb::MayBeString=nothing,
    fa::MayBeString=nothing, svg::MayBeString=nothing,
    tempdir::MayBeString=nothing, rotate::Bool=false)
    if isnothing(tempdir)
        tempdir = "."
    end
    tempfile = TempFile(tempdir)
    id, gffs, genome = doone(infile, tempfile)
    if rotate
        gffs, genome = trnF_start(gffs, genome)
    end

    if fa !== nothing
        open(FASTA.Writer, fa) do w
            write(w, FASTA.Record(id, genome))
        end
    end
    if gff !== nothing
        writeGFF(id, gffs, gff)
    else
        writeGFF(id, gffs, stdout)
    end
    if gb !== nothing
        writeGB(id, gffs, gb)
    end
    if svg !== nothing
        glength = length(genome)
        mRNAless = filter(x -> x.ftype != "mRNA" && x.ftype != "CDS", gffs)
        drawgenome(svg, id, glength, mRNAless)
    end
end

function emma_args(args::Vector{String}=ARGS)
    emma_args = ArgParseSettings(prog="emma", autofix_names=true)

    @add_arg_table! emma_args begin
        "--level", "-l"
        arg_type = String
        default = "info"
        help = "log level (info,warn,error,debug)"
        "--rotate"
        action = :store_true
        help = "begin with tRNA-Phe"
        "--gff"
        arg_type = String
        help = "output gff file"
        "--fa"
        arg_type = String
        help = "output FASTA file (only makes sense with --rotate)"
        "--svg"
        arg_type = String
        help = "output SVG image"
        "--gb"
        arg_type = String
        help = "output GB file"
        "--tempdir", "-t"
        arg_type = String
        help = "directory to use for temporary files [default is current directory]"
        "fastafiles"
        arg_type = String
        nargs = '+'
        action = :store_arg
        help = "fasta files to process"
    end

    emma_args.epilog = """
    If there are multiple FASTA files then (--gff|--fa|--svg|--gb)
    can refer to a directory in which case the respective files will be
    placed there under the original filename (with new extension). Otherwise
    the values will be used as a suffix.
    """
    parse_args(args, emma_args; as_symbols=true)
end


const LOGLEVELS = Dict("info" => Logging.Info, "debug" => Logging.Debug, "warn" => Logging.Warn,
    "error" => Logging.Error)

function main(args::Vector{String}=ARGS)

    function getout(accession, out, ext)
        if out === nothing
            return nothing
        end
        if isdir(out)
            return joinpath(out, basename(accession) * ext)
        end
        return accession * out
    end
    eargs = emma_args(args)
    llevel = get(LOGLEVELS, lowercase(eargs[:level]), Logging.Warn)
    global_logger(ConsoleLogger(stderr, llevel, meta_formatter=Logging.default_metafmt))

    fastafiles = eargs[:fastafiles]
    if length(fastafiles) == 1
        # all_files = all(maybefile, values(filtered_args))
        mainone(fastafiles[1], gff=eargs[:gff], fa=eargs[:fa], gb=eargs[:gb], svg=eargs[:svg], tempdir=eargs[:tempdir],
            rotate=eargs[:rotate])
    else
        for fasta in fastafiles
            accession = first(splitext(fasta))

            gff = getout(accession, eargs[:gff], ".gff")
            fa = getout(accession, eargs[:fa], ".fa")
            svg = getout(accession, eargs[:svg], ".svg")
            gb = getout(accession, eargs[:gb], ".gb")

            mainone(fasta, gff=gff, fa=fa, gb=gb, svg=svg, tempdir=eargs[:tempdir],
                rotate=eargs[:rotate])
        end
    end
end
