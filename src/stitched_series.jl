mutable struct StitchedSeries{TO,TI, N} <: AbstractArray{TO,N}
    img_top::AbstractArray{TI,N}
    img_bottom::AbstractArray{TI,N}
    tfms::Vector{CoordinateTransformations.Transformation} #This is a vector to allow for a separate transform for each z slice
    #tfm::CoordinateTransformations.Transformation
    y_sz::Int
    cached::AbstractArray{TO,2}
    cache_idxs::Tuple
    correct_bias::Bool
    sqrt_tfm::Bool
    flipy_bottom::Bool
end

function StitchedSeries(img_top::Array34{T}, img_bottom::Array34{T}, tfm::CoordinateTransformations.Transformation, out_type=Float64; correct_bias=true, sqrt_tfm=false, flipy_bottom=true) where {T}
    tfms = [tfm]
    for i= 2:size(img_bottom,3) #set all ztransforms to be equal
        push!(tfms, tfm) 
    end
    return StitchedSeries(img_top, img_bottom, tfms, out_type; correct_bias=correct_bias, sqrt_tfm=sqrt_tfm, flipy_bottom=flipy_bottom)
end

function StitchedSeries(img_top::Array34{T}, img_bottom::Array34{T}, tfms, out_type=Float64; correct_bias=true, sqrt_tfm=false, flipy_bottom=true) where {T}
    szt, szb = size(img_top), size(img_bottom)
    for i =1:4
        if szt[i] != szb[i]
            error("Input image size mismatch in dimension $i")
        end
    end
    if size(img_top,3) !== length(tfms)
        error("Input image size in the Z-slice dimension (3) does not match the number of transforms provided")
    end
    if sqrt_tfm && !correct_bias
        warn("Square root transformation will yield incorrect results if bias was not subtracted first")
    end
    #do a trial transform to get the y size (nevermind, just setting it to 2x input image y size)
    #temp = stitch(view(img_top,:, :, 1, 2), view(img_bottom,:,:,1,2), tfm, out_type; flip_y_img2=flipy_bottom)
    ss = StitchedSeries{out_type, T, ndims(img_top)}(img_top, img_bottom, tfms, 2*size(img_top,2), zeros(out_type, size(img_top,1), 2*size(img_top,2)), (ones(ndims(img_top)-2)...), correct_bias, sqrt_tfm, flipy_bottom)
    update_cache!(ss, (ones(Int, ndims(img_top)-2)...))
    return ss
end

size(A::StitchedSeries{T}) where {T} = (size(A.img_top,1), A.y_sz, size(A.img_top,3), size(A.img_top,4))
eltype(A::StitchedSeries{T}) where {T} = T
show(io::IO, A::StitchedSeries{T}) where {T} = print(io, "Stitched image series of size $(size(A)) mapped to element type $T\n")
show(io::IO, ::MIME"text/plain", A::StitchedSeries{T}) where {T} = show(io, A)

Base.IndexStyle(::Type{<:StitchedSeries}) = IndexCartesian()

function update_cache!(A::StitchedSeries{TO, TI, N}, inds::NTuple{N2, Int}) where {TO, TI, N, N2}
    slt = view(A.img_top, :, :, inds...)
    slb = view(A.img_bottom, :, :, inds...)
    tfm = A.tfms[inds[1]] #get transform for this slice
    if A.correct_bias && A.sqrt_tfm
        slt = sqrt.(Float64.(correctbias.(slt)))
	slb = sqrt.(Float64.(correctbias.(slb)))
    elseif A.correct_bias
        slt = correctbias.(slt)
	slb = correctbias.(slb)
    elseif A.sqrt_tfm
        slt = sqrt.(Float64.(slt))
        slb = sqrt.(Float64.(slb))
    end
    A.cached = stitch(slt, slb, tfm, TO; flip_y_img2=A.flipy_bottom)
    A.cache_idxs = inds
end

function slow_getindex(A::StitchedSeries, I)
    sl = A[:,:,I[3:end]...]
    return sl[I[1],I[2]]
end

getindex(A::StitchedSeries, dim1::Int, otherdims...) = slow_getindex(A, (dim1, otherdims...))
getindex(A::StitchedSeries, dim1, dim2::Int, otherdims...) = slow_getindex(A, (dim1, dim2, otherdims...))
getindex(A::StitchedSeries, dim1::Int, dim2::Int, otherdims...) = slow_getindex(A, (dim1, dim2, otherdims...))

getindex(A::StitchedSeries, dim1::Union{UnitRange, Colon}, dim2::Union{UnitRange,Colon}, dim3::Int, dim4::Int) = _getslice(A, dim1, dim2, (dim3,dim4))
getindex(A::StitchedSeries, dim1::Union{UnitRange, Colon}, dim2::Union{UnitRange,Colon}, dim3::Int) = _getslice(A, dim1, dim2, (dim3,))

#NOTE: duplicated in ZWarpedArray
function _getslice(A::StitchedSeries, dim1::Union{UnitRange, Colon}, dim2::Union{UnitRange,Colon}, queryidxs::NTuple{N,Int}) where {N}
#    if A.cache_zidx != dim3 || A.cache_tidx != dim4
#
#        if A.correct_bias
#            A.cached = stitch(correctbias(view(A.img_top, dim1, dim2, dim3, dim4)), correctbias(view(A.img_bottom,dim1,dim2,dim3,dim4)), A.tfm, T; flip_y_img2=A.flipy_bottom)
#        else
#            A.cached = stitch(view(A.img_top, dim1, dim2, dim3, dim4), view(A.img_bottom,dim1,dim2,dim3,dim4), A.tfm, T; flip_y_img2=A.flipy_bottom)
#        end
#        if A.sqrt_tfm
#            A.cached = sqrt.(A.cached)
#        end
#        A.cache_zidx = dim3
#        A.cache_tidx = dim4
#    end
#    return A.cached[dim1,dim2]
    if A.cache_idxs != queryidxs
        update_cache!(A, queryidxs)
    end
    return A.cached[dim1,dim2]
end

#returns multiple slices
#NOTE: duplicated in ZWarpedArray
function getindex(A::StitchedSeries, dim1::Union{UnitRange,Colon}, dim2::Union{UnitRange,Colon}, queryidxs)
    allidxs = [dim1;dim2;queryidxs...]
    idxs = []
    for i = 1:ndims(A)
        push!(idxs, _idxs(A, i, allidxs[i]))
    end
    output = zeros(T, map(length, idxs)...)
    _getindex!(output, A, idxs1, idxs2, queryidxs)
end

#NOTE: duplicated in ZWarpedArray
function _getindex!(prealloc, A::StitchedSeries{T}, idxs1, idxs2, queryidxs::NTuple{N}) where {T,N}
    itr = enumerate(Base.product(queryidxs))
    itrsz = size(itr)
    for (ipre, ia) in itr
        prealloc[:, :, ind2sub(itersz, ipre)...] = A[idxs1,idxs2,ia...]
    end
    return prealloc
end

StitchedSeries(img_top::ImageMeta, img_bottom::ImageMeta, tfm, out_type=Float64; kwargs...) = ImageMeta(StitchedSeries(data(img_top), data(img_bottom), tfm, out_type; kwargs...), properties(img_top))
StitchedSeries(img_top::AxisArray{T,4}, img_bottom::AxisArray{T,4}, tfm, out_type=Float64; kwargs...) where {T} = match_axisspacing(StitchedSeries(data(img_top),data(img_bottom),tfm,out_type; kwargs...), img_top)

