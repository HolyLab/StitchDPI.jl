module TeamedImaging

using CoordinateTransformations, Interpolations, PaddedViews, Unitful, AxisArrays, Images, MappedArrays
using BlockRegistration, RegisterOptimize
using NRRD, FileIO
using CachedSeries, ZWarpedArrays, InterlacedStacks

import Base: size, getindex, setindex!, show
import CachedSeries: update_cache!, cache, cache_idxs
import ZWarpedArrays: Array34

export stitch_tfm, stitch, mg_overlay, stitched_mg_overlay,# register_padmatched, full2full, 
        StitchedSeries,
        lazy_load, multiproc_write

include("util.jl")
include("preprocess.jl")
include("align.jl")
include("stitched_series.jl")
include("slicewise.jl")
include("load.jl")
include("write_to_disk.jl")

end
