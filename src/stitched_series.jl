mutable struct StitchedSeries{TO,TI} <: AbstractArray{TO,4}
    img_top::AbstractArray{TI,4}
    img_bottom::AbstractArray{TI,4}
    tfm::CoordinateTransformations.Transformation
    y_sz::Int
    cached::AbstractArray{TO,2}
    cache_zidx::Int
    cache_tidx::Int
    correct_bias::Bool
    sqrt_tfm::Bool
    flipy_bottom::Bool
end

correctbias(pix::Normed{UInt16, 16}, bias=100/2^16) = Normed{UInt16,16}(max(0.0, pix-bias))
correctbias(img::AbstractArray{Normed{UInt16, 16}, N}, bias=100/2^16) where {N} = mappedarray(x->correctbias(x), img)

sqrtimg(img::AbstractArray{T,N}) where {T,N} = mappedarray(x->sqrt(x), img)

function StitchedSeries(img_top::AbstractArray{T,4}, img_bottom::AbstractArray{T,4}, tfm, out_type=Float64; correct_bias=true, sqrt_tfm=false, flipy_bottom=true) where {T}
    szt, szb = size(img_top), size(img_bottom)
    for i =1:4
        if szt[i] != szb[i]
            error("Input image size mismatch in dimension $i")
        end
    end
    #do a trial transform to get the y size
    temp = stitch(view(img_top,:, :, 1, 2), view(img_bottom,:,:,1,2), tfm, out_type; flip_y_img2=flipy_bottom)
    ss = StitchedSeries{out_type, T}(img_top, img_bottom, tfm, size(temp,2), temp, 1, 2, correct_bias, sqrt_tfm, flipy_bottom)
    temp = ss[:,:,1,1]; #to update cached slice
    #ss.cached = ss[:,:,1,1]; #Not sure why the assignment is necessary
    #ss.cache_tidx = ss.cache_zidx = 1
    return ss
end

size(A::StitchedSeries{T}) where {T} = (size(A.img_top,1), A.y_sz, size(A.img_top,3), size(A.img_top,4))
eltype(A::StitchedSeries{T}) where {T} = T
show(io::IO, A::StitchedSeries{T}) where {T} = print(io, "Stitched image series of size $(size(A)) mapped to element type $T\n")
show(io::IO, ::MIME"text/plain", A::StitchedSeries{T}) where {T} = show(io, A)

Base.IndexStyle(::Type{<:StitchedSeries}) = IndexCartesian()

#idx_warn() = warn("Iterating through a StitchedSeries this way is very inefficient.  Query an entire frame (Colons in dims 1 and 2) instead")

function slow_getindex(A::StitchedSeries, I)
    sl = A[:,:,I[3], I[4]]
    return sl[I[1],I[2]]
end

getindex(A::StitchedSeries, dim1::Int, otherdims...) = slow_getindex(A, (dim1, otherdims...))
getindex(A::StitchedSeries, dim1, dim2::Int, otherdims...) = slow_getindex(A, (dim1, dim2, otherdims...))
getindex(A::StitchedSeries, dim1::Int, dim2::Int, otherdims...) = slow_getindex(A, (dim1, dim2, otherdims...))

function getindex(A::StitchedSeries{T}, dim1::Union{UnitRange, Colon}, dim2::Union{UnitRange,Colon}, dim3::Int, dim4::Int) where {T}
    if A.cache_zidx != dim3 || A.cache_tidx != dim4
        if A.correct_bias
            A.cached = stitch(correctbias(view(A.img_top, dim1, dim2, dim3, dim4)), correctbias(view(A.img_bottom,dim1,dim2,dim3,dim4)), A.tfm, T; flip_y_img2=A.flipy_bottom)
        else
            A.cached = stitch(view(A.img_top, dim1, dim2, dim3, dim4), view(A.img_bottom,dim1,dim2,dim3,dim4), A.tfm, T; flip_y_img2=A.flipy_bottom)
        end
        if A.sqrt_tfm
            A.cached = sqrt.(A.cached)
        end
        A.cache_zidx = dim3
        A.cache_tidx = dim4
    end
    return A.cached[dim1,dim2]
end

_idxs(A, dim, idxs) = idxs
_idxs(A, dim, idxs::Colon) = indices(A,dim)

function getindex(A::StitchedSeries{T}, dim1::Union{UnitRange,Colon}, dim2::Union{UnitRange,Colon}, dim3, dim4) where {T}
    idxs1 = _idxs(A, 1, dim1)
    idxs2 = _idxs(A, 2, dim2)
    idxs3 = _idxs(A, 3, dim3)
    idxs4 = _idxs(A, 4, dim4)
    output = zeros(T, length(idxs1), length(idxs2), length(idxs3), length(idxs4))
    _getindex!(output, A, idxs1, idxs2, idxs3, idxs4)
end

function _getindex!(prealloc, A::StitchedSeries{T}, idxs1, idxs2, idxs3, idxs4) where {T}
    for i3 = 1:length(idxs3)
        for i4 = 1:length(idxs4)
            sl = A[idxs1,idxs2,idxs3[i3],idxs4[i4]]
            prealloc[:, :, i3, i4] = sl
        end
    end
    return prealloc
end

axisspacing(A::AxisArray) = map(step, axisvalues(A))
match_axisspacing(B::AbstractArray{T,N}, A::IM) where {T,N, IM<:ImageMeta} = ImageMeta(match_axisspacing(B, data(A)), properties(A))
function match_axisspacing(B::AbstractArray{T,N}, A::AxisArray{T2,N}) where {T,T2,N}
    sp = axisspacing(A)
    nms = axisnames(A)
    newaxs = []
    for ax = 1:length(sp)
        u = unit(sp[ax])
        push!(newaxs, Axis{nms[ax]}(linspace(0.0*u, (size(B,ax)-1)*sp[ax], size(B,ax))))
    end
    return AxisArray(B, newaxs...)
end

StitchedSeries(img_top::ImageMeta, img_bottom::ImageMeta, tfm, out_type=Float64; kwargs...) = ImageMeta(StitchedSeries(data(img_top), data(img_bottom), tfm, out_type; kwargs...), properties(img_top))
StitchedSeries(img_top::AxisArray{T,4}, img_bottom::AxisArray{T,4}, tfm, out_type=Float64; kwargs...) where {T} = match_axisspacing(StitchedSeries(data(img_top),data(img_bottom),tfm,out_type; kwargs...), img_top)

