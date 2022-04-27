module StitchDPI

using CoordinateTransformations, Interpolations, PaddedViews, Unitful, AxisArrays, ImageMetadata, ImageBase, MappedArrays
const axes = Base.axes #for name conflict with AxisArrays

const Tform = CoordinateTransformations.Transformation

using RegisterQD
using NRRD, FileIO
using CachedArrays, ZWarpedArrays, InterleavedImages

using Mmap, Distributed #write_to_disk.jl

import Base: size, axes, getindex, setindex!, show
import CachedArrays: AbstractCachedArray,
                        update_cache!,
                        set_I!,
                        parent,
                        #cached_axes,
                        #noncached_axes,
                        axisspacing,
                        match_axisspacing

export stitch_tfm, stitch, mg_overlay, stitched_mg_overlay,# register_padmatched, full2full,
        StitchedSeries,
        lazy_load, multiproc_write

include("preprocess.jl")
include("align.jl")
include("stitched_series.jl")
include("slicewise.jl")
include("load.jl")
include("write_to_disk.jl")

end
