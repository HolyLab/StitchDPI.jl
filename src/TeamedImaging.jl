module TeamedImaging

using CoordinateTransformations, Interpolations, PaddedViews, Unitful, AxisArrays, Images, MappedArrays
using BlockRegistration, RegisterOptimize

import Base: size, getindex, show, eltype

export stitch_tfm, stitch, StitchedSeries, ZWarpedArray

include("util.jl")
include("preprocess.jl")
include("align.jl")
include("stitched_series.jl")
include("zwarpedarray.jl")

end
