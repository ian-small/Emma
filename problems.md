
*  gb.jl::34   vector not defined `second_cds = vector[findfirst(x -> x.fstart == string(cdsend + 2), vector)]`
* orfs.jl::63    no ncbi_trans_table        `translation = BioSequences.translate(genome.sequence[start:(nextstop-1)], code=ncbi_trans_table[translation_table])`

* Luxor.jl (SVG drawing) might be problematic in a multi-threaded system (Like EmmaServer.jl)
  (It has global dictionaries indexed by threadid *but!* julia (1.7) can now have executing tasks
  move to another thread [see here](https://docs.julialang.org/en/v1/manual/multi-threading/#man-task-migration))
