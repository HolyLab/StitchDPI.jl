#using CoordinateTransformations, Images, AxisArrays, Interpolations

#import Base: size, eltype, show, getindex
#include("preprocess.jl")
#include("util.jl")

mutable struct ZWarpedArray{TO,TI,N} <: AbstractArray{TO,N}
    parent::AbstractArray{TI,N}
    tfms::Vector{CoordinateTransformations.Transformation}
    cached::AbstractArray{TO,2}
    cache_idxs::Tuple
    correct_bias::Bool
    sqrt_tfm::Bool
end

function ZWarpedArray(img::Array34{T}, tfms, out_type=Float64; correct_bias=true, sqrt_tfm=false) where {T}
    if size(img,3) !== length(tfms)
        error("Input image size in the Z-slice dimension (3) does not match the number of transforms provided")
    end
    if sqrt_tfm && !correct_bias
        warn("Square root transformation will yield incorrect results if bias was not subtracted first")
    end
    za = ZWarpedArray{out_type, T, ndims(img)}(img, tfms, zeros(out_type, size(img)[1:2]...), (ones(ndims(img)-2)...), correct_bias, sqrt_tfm)
    update_cache!(za, (ones(Int, ndims(img)-2)...))
    return za
end

size(A::ZWarpedArray) = size(A.parent)
eltype(A::ZWarpedArray{TO}) where {TO} = TO
show(io::IO, A::ZWarpedArray{TO}) where {TO} = print(io, "ZWarpedArray of size $(size(A)) mapped to element type $TO\n")
show(io::IO, ::MIME"text/plain", A::ZWarpedArray{TO}) where {TO} = show(io, A)

Base.IndexStyle(::Type{<:ZWarpedArray}) = IndexCartesian()

function update_cache!(A::ZWarpedArray{TO, TI, N}, inds::NTuple{N2, Int}) where {TO, TI, N, N2}
    pslice = view(A.parent, :, :, inds...)
    tfm = A.tfms[inds[1]] #get transform for this slice
    if A.correct_bias && A.sqrt_tfm
        pslice = sqrt.(Float64.(correctbias.(pslice)))
    elseif A.correct_bias
        pslice = Float64.(correctbias.(pslice))
    elseif A.sqrt_tfm
        pslice = sqrt.(Float64.(pslice))
    end
    warped = warp(pslice, tfm)
    itp = interpolate(warped, BSpline(Linear()), OnGrid())
    etp = extrapolate(itp, NaN)
    A.cached = etp[indices(pslice)...] 
    A.cache_idxs = inds
end

function slow_getindex(A::ZWarpedArray, I)
    sl = A[:,:,I[3:end]...]
    return sl[I[1],I[2]]
end

getindex(A::ZWarpedArray, dim1::Int, otherdims...) = slow_getindex(A, (dim1, otherdims...))
getindex(A::ZWarpedArray, dim1, dim2::Int, otherdims...) = slow_getindex(A, (dim1, dim2, otherdims...))
getindex(A::ZWarpedArray, dim1::Int, dim2::Int, otherdims...) = slow_getindex(A, (dim1, dim2, otherdims...))

getindex(A::ZWarpedArray, dim1::Union{UnitRange, Colon}, dim2::Union{UnitRange,Colon}, dim3::Int, dim4::Int) = _getslice(A, dim1, dim2, (dim3,dim4))
getindex(A::ZWarpedArray, dim1::Union{UnitRange, Colon}, dim2::Union{UnitRange,Colon}, dim3::Int) = _getslice(A, dim1, dim2, (dim3,))

function _getslice(A::ZWarpedArray, dim1::Union{UnitRange, Colon}, dim2::Union{UnitRange,Colon}, queryidxs::NTuple{N,Int}) where {N}
    if A.cache_idxs != queryidxs
        update_cache!(A, queryidxs)
    end
    return A.cached[dim1,dim2]
end

#returns multiple slices
function getindex(A::ZWarpedArray, dim1::Union{UnitRange,Colon}, dim2::Union{UnitRange,Colon}, queryidxs)
    allidxs = [dim1;dim2;queryidxs...]
    idxs = []
    for i = 1:ndims(A)
        push!(idxs, _idxs(A, i, allidxs[i]))
    end
    output = zeros(T, map(length, idxs)...)
    _getindex!(output, A, idxs1, idxs2, queryidxs)
end

function _getindex!(prealloc, A::ZWarpedArray{T}, idxs1, idxs2, queryidxs::NTuple{N}) where {T,N}
    itr = enumerate(Base.product(queryidxs))
    itrsz = size(itr)
    for (ipre, ia) in itr
        prealloc[:, :, ind2sub(itersz, ipre)...] = A[idxs1,idxs2,ia...]
    end
    return prealloc
end

ZWarpedArray(img::ImageMeta, tfms, out_type=Float64; kwargs...) = ImageMeta(ZWarpedArray(data(img), tfms, out_type; kwargs...), properties(img))
ZWarpedArray(img::AxisArray, tfms, out_type=Float64; kwargs...) = match_axisspacing(ZWarpedArray(data(img),tfms,out_type; kwargs...), img)
