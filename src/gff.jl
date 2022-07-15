struct GFF
    source::String
    ftype::String
    fstart::String
    fend::String
    score::String
    strand::Char
    phase::String
    attributes::String
end

function HMMmatch2GFF(cds::HMMmatch, genome_length::Integer)
    bits = split(cds.orf, "_")
    ends = split(bits[3], "-")
    strand = bits[2][1]
    start = parse(Int, ends[1])
    finish = parse(Int, ends[2])
    hmmstart = start + 3*(cds.ali_from-1)
    hmmfinish = start + 3*(cds.ali_to - cds.ali_from + 1) - 1
    if strand == '-'
        tmp = start
        start = reverse_complement(finish, genome_length)
        finish = reverse_complement(tmp, genome_length)
        tmp = hmmstart
        hmmstart = reverse_complement(hmmfinish, genome_length)
        hmmfinish = reverse_complement(tmp, genome_length)
    end
    return [GFF("Emma", "ORF", string(start), string(finish), string(cds.Evalue), strand, "0", "Name=" * cds.query * "_ORF"),
        GFF("Emma", "HMM alignment", string(hmmstart), string(hmmfinish), string(cds.Evalue), strand, "0", "Name=" * cds.query * "_HMM")]
end

function CMAlignment2GFF(trn::CMAlignment_trn, glength::Integer)
    startstring = trn.tstrand =='+' ? string(trn.tfrom) : string(reverse_complement(trn.tto, glength))
    finishstring = trn.tstrand =='+' ? string(trn.tto) : string(reverse_complement(trn.tfrom, glength))
    attributes = "Name=" * trn.query * "-" * trn.anticodon
    if trn.polyA; attributes *= ";Note=Sequence completed by polyadenylation";end
    return GFF("Emma", "tRNA match", startstring, finishstring, string(trn.Evalue), trn.tstrand, "0", attributes)
end

function rRNA2GFF(rrn::rRNA, glength::Integer)
    #sort starts and ends by Evalue
    sort!(rrn.start, by = x -> x.Evalue)
    sort!(rrn.stop, by = x -> x.Evalue)
    rrnstart = rrn.start[1]
    rrnstop = rrn.stop[1]
    startstring = rrnstop.tstrand =='+' ? string(rrnstart.ali_from) : string(reverse_complement(rrnstop.tto, glength))
    finishstring = rrnstop.tstrand =='+' ? string(rrnstop.tto) : string(rrnstart.ali_from)
    attributes = "Name=" * rrnstop.query
    return GFF("Emma", "rRNA match", startstring, finishstring, string(min(rrnstart.Evalue, rrnstop.Evalue)), rrnstop.tstrand, "0", attributes)
end

function writeGFF(outfile::String, id::String, genome_length::Integer, cds_matches::Vector{HMMmatch},
             trn_matches::Vector{CMAlignment_trn}, rRNAs::NamedTuple{(:rrnL, :rrnS), Tuple{rRNA, rRNA}})

    function writeone(out::IO, gff::GFF)
        write(out, join([id, gff.source, gff.ftype,gff.fstart,gff.fend,gff.score,gff.strand,gff.phase,gff.attributes], "\t"), "\n")
    end

    open(outfile, "w") do out
        for cds in cds_matches
            gffs = HMMmatch2GFF(cds, genome_length)
            for gff in gffs
                writeone(out, gff)
            end
        end
        for trn in trn_matches
            gff = CMAlignment2GFF(trn, genome_length)
            writeone(out, gff)
        end
        gff = rRNA2GFF(rRNAs.rrnL, genome_length)
        writeone(out, gff)
        gff = rRNA2GFF(rRNAs.rrnS, genome_length)
        writeone(out, gff)
    end
end
