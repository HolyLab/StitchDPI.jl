#Assume that float values are also normed
correctbias(pix::T, bias=100/2^16) where T<:Union{Normed{UInt16, 16}, Float32, Float64}  = T(max(0.0, pix-bias))
correctbias(img::AbstractArray{T, N}, bias=100/2^16) where {T<:Union{Normed{UInt16, 16}, Float32, Float64}, N} = mappedarray(x->correctbias(x), img)

sqrtimg(img::AbstractArray{T,N}) where {T,N} = mappedarray(x->sqrt(x), img)
