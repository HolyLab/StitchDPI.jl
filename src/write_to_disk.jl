#(When threading works, not yet):
#       To get the benefit of threading start julia with multiple threads:
#       JULIA_NUM_THREADS=x julia (where you choose x)

#A multiprocess write method for saving large image sequences in NRRD format
#This function is useful when retrieving a stack from img takes a lot of time due to lazy computation
function multiproc_write(fname::AbstractString, img::AbstractArray{T}, out_type=T) where {T}
    mmapa = prep_nrrd_write(fname, img, out_type)
    #pp = out_type == T ? x->x : x->out_type.(x)
    multiproc_write!(mmapa, img)
end

function prep_nrrd_write(fname::AbstractString, img::AbstractArray{T,N}, out_type=T) where {T,N}
    bname, ext = splitext(fname)
    if !isempty(ext) && !in(ext, [".nhdr", ".nrrd"])
        error("Must write in NRRD format (.nhdr or .nrrd extension)")
    end
    mmapa = Mmap.mmap(bname*".raw", Array{out_type, N}, size(img); shared=true)
    header = NRRD.headerinfo(out_type, axes(img))
    header["datafile"] = bname*".raw"
    open(fname, "w") do io
        write(io, magic(format"NRRD"))
        NRRD.write_header(io, "0004", header)
    end
    return mmapa
end

#Multithreading would be nice, but currently this doesn't seem to work well, perhaps because IO operations and external C library calls are restricted to one thread (see julia parallel docs as off march 2018)
function multiproc_write!(mmapA::AbstractArray{T,4}, img::AbstractArray{T2,4}) where {T, T2}
    @assert size(mmapA) == size(img)
    np = nworkers()
    if np == 0
        error("No workers available.  Add them with addprocs() (also make sure code and the input image are loaded on the worker processes")
    end
    print("Defining image and output array on other processes...\n")
    if np > 1
        #@eval @everywhere mmapA = $(mmapA)
        @eval @everywhere img = $(img)
    end
    #nthreads = Base.Threads.nthreads()
    #print("Writing image with $(nthreads) threads\n")
    print("Writing image with $(np) procs\n")
    cur_i = 1
    while cur_i<=size(img,4)
        i_start = cur_i
        cur_i = min(size(img,4), i_start + np-1)
        idxs = [i_start:cur_i...]
        rrefs = []
        ws = workers()
#    Threads.@threads for i=1:size(img,4)
        for i=1:length(idxs)
            push!(rrefs, remotecall(getindex, ws[i], img, :,:,:,idxs[i]))
        end
        for i=1:length(idxs)
            mmapA[:,:,:,idxs[i]] = fetch(rrefs[i])
        end
        print("Wrote $cur_i stacks out of $(size(img,4)) so far.\n")
        cur_i+=1
    end
    print("Finished\n")
end
