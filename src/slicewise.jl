#This produces a StitchedSeries that incorporates both the stitch transform and slicewise transformations that compensate for
#wiggle in the field of view resulting from fast dynamic piezo operation
#slice_tfms is the set of slicewise transformations that must be applied to img_top (and also to img_bottom after it is aligned to the space of img_top)
#octr is the point around which rotations encoded by slice_tfms are centered.  This is the (x,y) center of the image sequence used for calibration
#This function assumes that sequence was acquired on camera1.  If the (x,y) size of that sequence differes from the (x,y) sizes of the top and bottom
#images then the default octr will not be optimal.  Override it with the true center
function stitch_slicewise(img_top, img_bottom, st_tfm, slice_tfms, out_type=Float64; octr=RegisterQD.center(view(img_top,:,:,1,1)), correct_bias=true, sqrt_tfm=false, flip_y_bottom=true)
    tfms_top, tfms_bottom = recenter_tfms(slice_tfms, st_tfm, img_top, octr)
    za_top = ZWarpedArray(img_top, tfms_top, out_type; correct_bias=false, sqrt_tfm=false)
    return StitchedSeries(za_top, img_bottom, tfms_bottom, Float64; correct_bias=correct_bias, sqrt_tfm=sqrt_tfm, flip_y_bottom=flip_y_bottom)
end

#assume that slice_tfms are already centered on the calibration images, so need to be recentered by the difference between that point and the
#center of camera1's image.
function recenter_tfms(slice_tfms, st_tfm, cam1, octr)
    tfmstop = []
    tfmsbottom = []
    cam1ctr = RegisterQD.center(view(cam1, :,:,1,1)).-octr
    #The below calculation of cam2ctr isn't perfect but is close enough that it's not noticeable by eye
    cam2ctr = cam1ctr + [0.0; -(size(cam1,2)+1)/2] #The center must move negatively in y (assumes top and bottom strips have equal y sz and perfect alignment)
    #Future: maybe for cam2ctr we can just apply the stitching transform to the center of the image? watch for sign errors
    for i =1:length(slice_tfms)
        push!(tfmstop, recenter(slice_tfms[i], cam2ctr))
        push!(tfmsbottom, st_tfm ∘ recenter(slice_tfms[i], cam2ctr))
    end
    return tfmstop, tfmsbottom
end
