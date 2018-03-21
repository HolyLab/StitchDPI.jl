const Tform = CoordinateTransformations.Transformation

mutable struct StitchedSeries{TO,TI,N} <: CachedSeries2D{TO,TI,N}
    img_top::AbstractArray{TI,N}
    img_bottom::AbstractArray{TI,N}
    tfms::Vector{Tform} #This is a vector to allow for a separate transform for each z slice
    y_sz::Int
    cached::AbstractArray{TO,2}
    cache_idxs::Tuple
    correct_bias::Bool
    sqrt_tfm::Bool
    flipy_bottom::Bool
end

cache(A::StitchedSeries) = A.cached
cache_idxs(A::StitchedSeries) = A.cache_idxs

function StitchedSeries(img_top::Array34{T}, img_bottom::Array34{Normed{UInt16,16}}, tfm, out_type=Float64; kwargs...) where {T<:Union{Float32, Float64}}
    StitchedSeries(img_top, mappedarray(T, img_bottom), tfm, out_type; kwargs...)
end

function StitchedSeries(img_top::Array34{T}, img_bottom::Array34{T}, tfm::Tform, out_type=Float64; correct_bias=true, sqrt_tfm=false, flip_y_bottom=true) where {T}
    tfms = [tfm]
    for i= 2:size(img_bottom,3) #set all ztransforms to be equal
        push!(tfms, tfm) 
    end
    return StitchedSeries(img_top, img_bottom, tfms, out_type; correct_bias=correct_bias, sqrt_tfm=sqrt_tfm, flip_y_bottom=flip_y_bottom)
end

#function StitchedSeries(img_top::Array34{T}, img_bottom::Array34{T}, tfms::Vector{T2}, out_type=Float64; correct_bias=true, sqrt_tfm=false, flip_y_bottom=true) where {T, T2<:Tform}
function StitchedSeries(img_top::Array34{T}, img_bottom::Array34{T}, tfms, out_type=Float64; correct_bias=true, sqrt_tfm=false, flip_y_bottom=true) where {T}
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
    #Cache y size is always 2x input image y size (sometimes this will result in harmless empty space at the edge of the field of view)
    ss = StitchedSeries{out_type, T, ndims(img_top)}(img_top, img_bottom, tfms, 2*size(img_top,2), zeros(out_type, size(img_top,1), 2*size(img_top,2)), (ones(ndims(img_top)-2)...), correct_bias, sqrt_tfm, flip_y_bottom)
    update_cache!(ss, (ones(Int, ndims(img_top)-2)...))
    return ss
end

function update_cache!(A::StitchedSeries{TO, TI, N}, inds::NTuple{N2, Int}) where {TO, TI, N, N2}
    slt = view(A.img_top, :, :, inds...)
    slb = view(A.img_bottom, :, :, inds...)
    tfm = A.tfms[inds[1]] #get transform for this slice
    if A.correct_bias
        slt = mappedarray(correctbias, slt)
        slb = mappedarray(correctbias, slb)
    end
    if A.sqrt_tfm
        A.cached = sqrt.(stitch(slt[:,:], slb[:,:], tfm, TO; flip_y_bottom=A.flipy_bottom))
    else
        A.cached = stitch(slt[:,:], slb[:,:], tfm, TO; flip_y_bottom=A.flipy_bottom)
    end
    A.cache_idxs = inds
end

size(A::StitchedSeries{T}) where {T} = (size(A.img_top,1), A.y_sz, size(A.img_top,3), size(A.img_top,4))
show(io::IO, A::StitchedSeries{T}) where {T} = print(io, "Stitched image series of size $(size(A)) mapped to element type $T\n")
show(io::IO, ::MIME"text/plain", A::StitchedSeries{T}) where {T} = show(io, A)

#Yuck, had to write two versions to resolve ambiguity error (not sure why because AxisArray{T,4} should be more specific than AbstractArray{T,4}?)
function StitchedSeries(img_top::ImageMeta, img_bottom::ImageMeta, tfm::Tform, out_type=Float64; kwargs...)
    ImageMeta(StitchedSeries(data(img_top), data(img_bottom), tfm, out_type; kwargs...), properties(img_top))
end
function StitchedSeries(img_top::ImageMeta, img_bottom::ImageMeta, tfm::Vector{Tform}, out_type=Float64; kwargs...)
    ImageMeta(StitchedSeries(data(img_top), data(img_bottom), tfm, out_type; kwargs...), properties(img_top))
end
function StitchedSeries(img_top::AxisArray{T,4}, img_bottom::AxisArray{T,4}, tfm::Tform, out_type=Float64; kwargs...) where {T}
    match_axisspacing(StitchedSeries(data(img_top),data(img_bottom),tfm,out_type; kwargs...), img_top)
end
function StitchedSeries(img_top::AxisArray{T,4}, img_bottom::AxisArray{T,4}, tfm::Vector{Tform}, out_type=Float64; kwargs...) where {T}
    match_axisspacing(StitchedSeries(data(img_top),data(img_bottom),tfm,out_type; kwargs...), img_top)
end
