trunc_above(v::T,thr) where {T} = ifelse(v<thr || isnan(v), v, T(thr))
nanz(x) = ifelse(isnan(x), zero(x), x)
#Assume that float values are also normed
correctbias(pix::T, bias=100/2^16) where T<:Union{Normed{UInt16, 16}, Float32, Float64}  =
    T(max(0.0, pix-bias))
correctbias(img::AbstractArray{T, N}, bias=100/2^16) where {T<:Union{Normed{UInt16, 16},Float32,Float64},N} =
    mappedarray(x->correctbias(x), img)
sqrtimg(img::AbstractArray{T,N}) where {T,N} = mappedarray(x->sqrt(x), img)

function ypad(img, padded_size::Int; flip_y = true, fillval = eltype(img)(NaN), pad_side=:both)
    @assert ndims(img) == 2
    @assert padded_size >= size(img,2)
    npix = div(padded_size - size(img,2), 2)
    if flip_y
        img = flipy(img)
    end
    first_idx = (1,1)
    if pad_side == :both
        first_idx = (1,npix+1)
    elseif pad_side == :top
        first_idx = (1,2*npix+1)
    elseif pad_side != :bottom
        error("Unrecognized padding argument $pad_side")
    end
    return PaddedView(fillval, img, (size(img,1),padded_size), first_idx)
end

flipy(img::AbstractArray{T,2}) where {T} = view(img, :, reverse(last(axes(img))))

#Crops the image to the minimal rectangle without val along any edge
function crop_vals(img::AbstractArray{T,2}, val, dims::Tuple=(1,2)) where {T}  
    if dims == (1,)
        return crop_vals_x(img, val)
    elseif dims == (2,)
        return crop_vals_y(img, val)
    elseif dims == (1,2)
        return crop_vals_xy(img, val)
    else
        error("Unrecognized dims Tuple")
    end
end

function crop_vals_xy(img::AbstractArray{T,2}, val) where {T}
    xmax = ymax = 1
    xmin = size(img,1)
    ymin = size(img,2)
    for y = 1:size(img,2)
        for x = 1:size(img,1)
            if !isequal(img[x,y], val)
                xmin = min(x,xmin)
                xmax = max(x,xmax)
                ymin = min(y,ymin)
                ymax = max(y,ymax)
            end
        end
    end
    if xmax < xmin || ymax < ymin
        error("All elements of image seem to be $val")
    end
    return view(img, xmin:xmax, ymin:ymax)
end

function crop_vals_x(img::AbstractArray{T,2}, val) where {T}
    xmax = 1
    xmin = size(img,1)
    min_found = false
    max_found = false
    for x = 1:size(img,1)
        for y = 1:size(img,2)
            if !isequal(img[x,y], val)
                xmin = x
                min_found = true
                break
            end
        end
        if min_found
            break
        end
    end
    for x = size(img,1):-1:1
        for y = 1:size(img,2)
            if !isequal(img[x,y], val)
                xmax = x
                max_found = true
                break
            end
        end
        if max_found
            break
        end
    end
    if xmax < xmin
        error("All elements of image seem to be $val")
    end
    return view(img, xmin:xmax, :)
end

function crop_vals_y(img::AbstractArray{T,2}, val) where {T}
    ymax = 1
    ymin = size(img,2)
    min_found = false
    max_found = false
    for y = 1:size(img,2)
        for x = 1:size(img,1)
            if !isequal(img[x,y], val)
                ymin = y
                min_found = true
                break
            end
        end
        if min_found
            break
        end
    end
    for y = size(img,2):-1:1
        for x = 1:size(img,1)
            if !isequal(img[x,y], val)
                ymax = y
                max_found = true
                break
            end
        end
        if min_found
            break
        end
    end
    if ymax < ymin
        error("All elements of image seem to be $val")
    end
    return view(img, :, ymin:ymax)
end

function nanplus!(prealloc, img1, img2)
    for i in eachindex(img1)
        val1 = img1[i]
        val2 = img2[i]
        isn1 = isnan(val1)
        if isn1 || isnan(val2)
            prealloc[i] = ifelse(isn1, val2 , val1)
        else
            prealloc[i] = val1 + val2
        end
    end
    return prealloc
end

nanplus(img1::AbstractArray{T}, img2::AbstractArray{T}) where {T<:AbstractFloat} = nanplus!(similar(img1), img1, img2) 
nanplus(img1, img2) = img1.+img2

##Not used
##Crops the image to the specified y size (symmetrically crops top and bottom)
#function crop_y(img::AbstractArray{T,2}, ysz::Int) where {T}
#    @assert iseven(ysz)
#    cursz = size(img,2)
#    @assert iseven(cursz) && cursz >= ysz
#    halfcrop = div(cursz-ysz, 2)
#    return view(img, :, halfcrop+1:cursz-halfcrop)
#end
