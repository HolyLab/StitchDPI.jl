module TeamedImaging

using CoordinateTransformations, Interpolations, PaddedViews, Unitful, AxisArrays, Images, MappedArrays
using BlockRegistration, RegisterOptimize
using NRRD, FileIO


import Base: size, getindex, setindex!, show, eltype

export stitch_tfm, stitch, warp_and_resample, register_padmatched, full2full, 
        StitchedSeries, ZWarpedArray, InterlacedStackSeries,
        lazy_load, multiproc_write

include("util.jl")
include("preprocess.jl")
include("align.jl")
include("cached_series_2d.jl")
include("stitched_series.jl")
include("zwarpedarray.jl")
include("slicewise.jl")
include("interlaced.jl")
include("load.jl")
include("write_to_disk.jl")

end
