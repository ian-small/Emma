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
    cdsstart = cds.ali_from
    cdsstop = cds.ali_to + 2
    if cds.strand == '-'
        tmp = cdsstart
        cdsstart = reverse_complement(cdsstop, genome_length)
        cdsstop = reverse_complement(tmp, genome_length)
    end
    if cdsstart > genome_length
        cdsstart = mod1(cdsstart, genome_length)
        cdsstop = mod1(cdsstop, genome_length)
    end
    return GFF("Emma", "CDS", string(cdsstart), string(cdsstop), string(cds.Evalue), cds.strand, "0", "Name=" * cds.query * "_HMM")
end

function CMAlignment2GFF(trn::CMAlignment_trn, glength::Integer)
    startstring = trn.tstrand =='+' ? string(trn.tfrom) : string(reverse_complement(trn.tto, glength))
    finishstring = trn.tstrand =='+' ? string(trn.tto) : string(reverse_complement(trn.tfrom, glength))
    attributes = "Name=" * trn.query * "-" * trn.anticodon
    if trn.polyA; attributes *= ";Note=Sequence completed by polyadenylation";end
    return GFF("Emma", "tRNA", startstring, finishstring, string(trn.Evalue), trn.tstrand, "0", attributes)
end

function rRNA2GFF(rrn::rRNA, glength::Integer)
    start, stop = gettermini(rrn, glength)
    attributes = "Name=" * rrn.stop[1].query
    evalue = getevalue(rrn)
    return GFF("Emma", "rRNA", string(start), string(stop), string(evalue), rrn.stop[1].tstrand, "0", attributes)
end

function writeGFF(outfile::String, id::String, genome_length::Integer, cds_matches::Vector{HMMmatch},
             trn_matches::Vector{CMAlignment_trn}, rRNAs::NamedTuple{(:rrnL, :rrnS), Tuple{rRNA, rRNA}})

    function writeone(out::IO, gff::GFF)
        write(out, join([id, gff.source, gff.ftype,gff.fstart,gff.fend,gff.score,gff.strand,gff.phase,gff.attributes], "\t"), "\n")
    end

    open(outfile, "w") do out
        for cds in cds_matches
            gff = HMMmatch2GFF(cds, genome_length)
            writeone(out, gff)
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