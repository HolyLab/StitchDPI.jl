module TeamedImaging

using CoordinateTransformations, Interpolations, PaddedViews, Unitful, AxisArrays, Images, MappedArrays
using BlockRegistration, RegisterOptimize

import Base: size, getindex, show, eltype

export stitch_tfm, stitch, StitchedSeries

include("align.jl")
include("stitched_series.jl")

end
