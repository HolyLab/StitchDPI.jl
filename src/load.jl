#Lazy loading of stitched high-speed experiment (applies both a stitching tfm and a slicewise calibration)
function lazy_load(cam1::AbstractArray, cam2::AbstractArray, octr, stitch_tfm::T, fwd_tfms, bck_tfms;
                        sqrt_tfm=false, correct_bias=true) where{T<:CoordinateTransformations.Transformation}
    #the match_axisspacing call works around a breakage in AxisArrays in this case
    cam1fwd = match_axisspacing(view(cam1.data.data,:,:,:,1:2:size(cam1,4)), cam1.data)
    cam2fwd = match_axisspacing(view(cam2.data.data,:,:,:,1:2:size(cam2,4)), cam2.data)
    cam1bck = match_axisspacing(view(cam1.data.data,:,:,:,2:2:size(cam1,4)), cam1.data)
    cam2bck = match_axisspacing(view(cam2.data.data,:,:,:,2:2:size(cam2,4)), cam2.data)

    #(This workaround isn't ideal because it doesn't create a StepRangeLen, so the notion of uniform axis spacing is lost
    #cam1fwd = view(cam1,:,:,:,[1:2:size(cam1,4)...])
    #cam2fwd = view(cam2,:,:,:,[1:2:size(cam2,4)...])
    #cam1bck = view(cam1,:,:,:,[2:2:size(cam1,4)...])
    #cam2bck = view(cam2,:,:,:,[2:2:size(cam2,4)...])
i
    ssfwd = TeamedImaging.stitch_slicewise(cam1fwd, cam2fwd, stitch_tfm, fwd_tfms, Float64;
        octr=octr, correct_bias=correct_bias, sqrt_tfm=sqrt_tfm, flip_y_bottom=true)
    ssbck = TeamedImaging.stitch_slicewise(cam1bck, cam2bck, stitch_tfm, bck_tfms, Float64;
        octr=octr, correct_bias=correct_bias, sqrt_tfm=sqrt_tfm, flip_y_bottom=true)
    ssi = InterlacedStackSeries(ssfwd, ssbck)
    return ssi
end
