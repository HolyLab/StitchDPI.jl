correctbias(pix::Normed{UInt16, 16}, bias=100/2^16) = Normed{UInt16,16}(max(0.0, pix-bias))
correctbias(img::AbstractArray{Normed{UInt16, 16}, N}, bias=100/2^16) where {N} = mappedarray(x->correctbias(x), img)

sqrtimg(img::AbstractArray{T,N}) where {T,N} = mappedarray(x->sqrt(x), img)
