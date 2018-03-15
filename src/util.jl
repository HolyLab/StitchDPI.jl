const Array34{T} = Union{AbstractArray{T,3}, AbstractArray{T,4}}

function ypad(img, padded_size::Int; flip_y = true, fillval = eltype(img)(NaN), pad_side=:both)
    @assert ndims(img) == 2
    @assert padded_size >= size(img,2)
    @assert iseven(padded_size - size(img,2))
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
    return PaddedView(fillval, img, (size(img,1),size(img,2)+2*npix), first_idx) # (2, 2) is the location of the first element from a
end

flipy(img::AbstractArray{T,2}) where {T} = view(img, :, reverse(last(indices(img))))

function crop_vals(img::AbstractArray{T,2}, val) where {T}
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

function nanplus!(prealloc, img1, img2)
    for i in eachindex(img1)
        val1 = img1[i]
        val2 = img2[i]
        if isnan(val1) || isnan(val2)
            prealloc[i] = isnan(val1) ? val2 : val1
        else
            prealloc[i] = val1 + val2
        end
    end
    return prealloc
end

nanplus(img1, img2) = nanplus!(similar(img1), img1, img2)

#fake split images with requested split
#returns the top full/split images first, then bottom
#ysz_roi is the size of one camera's ROI (so total combined image will have dimension 2x ysz_roi)
function fake_split(ysz_chip, ysz_roi; frac_overlap=0.0, xsz=20)
    @assert iseven(ysz_chip)
    @assert iseven(ysz_roi)
    @assert 0.0 <= frac_overlap <= 1.0
    npix_extra = div(ysz_chip - ysz_roi, 2) #unused vertical roi  (symmetric above and below)
    npix_overlap = ceil(Int, frac_overlap * ysz_roi) #number of pixels overlapping between shared images (zero would mean no redundant pixels)
    chip_halfsz = div(ysz_chip,2)
    ysz_full = 2*ysz_roi - npix_overlap + 2*npix_extra
    full_img = rand(xsz, ysz_full)
    full_top = full_img[:,1:ysz_chip]
    full_bottom = full_img[:,(ysz_full-ysz_chip+1):ysz_full]
    img_top = full_top[:,(npix_extra+1):(npix_extra+ysz_roi)]
    img_bottom = full_bottom[:,(npix_extra+1):(npix_extra+ysz_roi)]
    return full_top, img_top, full_bottom, img_bottom
end


_idxs(A, dim, idxs) = idxs
_idxs(A, dim, idxs::Colon) = indices(A,dim)


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
