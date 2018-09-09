#shared image alignment procedure on the microscope
#1. Install KEM.  Make sure the KEM is positioned so that the image stays crisp on camera2 (may want to compare with a dichroic to make sure)
#2. Set camera ROIs as desired (same size for both cameras)
#3. Move camera1 up so that the KEM line corresponds with the bottom side of the ROI (larger ROIs will require moving the camera farther)
#4. Move camera2 left so that the KEM line corresponds with the bottom side of the ROI
#5. Note the ROI settings and change the ROI to full chip
#6. Install 50/50 dichroic instead of KEM
#7. Take a snapshot of beads with camera2
#8. Reset ROIs as noted in #5, replace KEM mirror, and take a snapshot with each camera.  Make sure the beads don't move during this sequence!

#shared image alignment computational procedure
#1.  Compute transform to align camera2's flipped padded half-image (step 8) with its flipped full image (step 7) (should be almost aligned already)
#2.  Compute transform to align camera2's flipped full image (step 7) with camera1's padded half-image (step 8)
#3.  Compose transform from 2 with transform from 1 and call recenter() on the composed transform.
#4.  For each camera2 half-image:
    # -flip y
    # -pad vertically to full chip size
    # -apply transform from #3 
    # -nansum with camera1's corresponding padded half-image
    # -find minimal bounding box that contains all non-nan pixels and return image cropped to that region

function stitch(img_top_fixed::AbstractArray{T,2}, img_bottom_moving::AbstractArray{T,2}, tfm, out_type::DataType; kwargs...) where {T}
    if T!=out_type
        img_top_fixed = out_type.(img_top_fixed)
        img_bottom_moving = out_type.(img_bottom_moving)
    end
    return stitch(img_top_fixed, img_bottom_moving, tfm; kwargs...)
end

#Stitches img_top_fixed and img_bottom_moving using a transformation tfm that, when applied to img_bottom_moving, aligns it with img_top_fixed.
#Since img_top_fixed and img_bottom_moving are expected to have little-to-no overlap they are padded before applying the transformation.
#The transformed img_bottom_moving is added to img_top_fixed and the interesting (non-NaN) region of the summed image is returned.
#For most accurate stitching you must provide ysz_full, the size in the y dimension of full images used to compute tfm with the stitch_tfm function
function stitch(img_top_fixed::AbstractArray{T,2}, img_bottom_moving::AbstractArray{T,2}, tfm; flip_y_bottom = true, ysz_full=2048) where {T}
    yroi = size(img_top_fixed,2)
    @assert size(img_bottom_moving,2) == yroi
    yextra = ceil(Int, abs(tfm.v[2]))
    #fillvals should be NaN, but currently interpolations doesn't handle it well (imputes NaNs at gridpoints)
    img_bottom_moving_pre = ypad(img_bottom_moving, ysz_full; pad_side=:both, flip_y=flip_y_bottom, fillval=0.0)
    img_top_fixed = ypad(img_top_fixed, ysz_full; pad_side=:both, flip_y=false, fillval=0.0)
    img_bottom_moving = warp_and_resample(img_bottom_moving_pre, tfm)
    imgsummed = nanplus(img_top_fixed,img_bottom_moving)
    ystart = div(ysz_full,2) - div(yroi,2) + 1
    return view(imgsummed, :, ystart:ystart+2*yroi-1)
end

#returns a transform that stitches a pair of corresponding shared images into a single image when passed to the stitch function
#Assumes that moving_roi contains the bottom half of the image and that it needs to be flipped vertically to align with fixed_roi.
function stitch_tfm(moving_full, fixed_roi, moving_roi; trunc_frac=0.1, kwargs...)
    @assert size(fixed_roi) == size(moving_roi)

    #bright bead clumps can bias calculation, so truncate them at this value:
    trunc_thr = trunc_frac * maxfinite(fixed_roi)
    moving_roi = trunc_above.(moving_roi, trunc_thr)
    fixed_roi = trunc_above.(fixed_roi, trunc_thr)
    moving_full = trunc_above.(moving_full, trunc_thr)
    
    #flip and pad
    moving_full = flipy(moving_full)
    fixed_roi_pad = ypad(fixed_roi, size(moving_full,2); flip_y = false, fillval=NaN)
    moving_roi_pad = ypad(moving_roi, size(moving_full,2); flip_y = true, fillval=NaN)

    #################### Moving trip to moving full
    print("Registering the moving strip to the moving full image\n")
    mxshift_moving = [50; size(moving_roi,2)+50]
    m_roi_pwr = sum(abs2.(moving_roi_pad[.!(isnan.(moving_roi_pad))]))
    moving_pre_rigid, pre_rigid_mm = qd_rigid(moving_full, moving_roi_pad, mxshift_moving, [pi/60], [pi/10000]; thresh=0.85*m_roi_pwr, kwargs...)
    @show moving_pre_rigid

    #refine it further, allowing full affine
    moving_pre_tfm, pre_mm = qd_affine(moving_full, moving_roi_pad, [50;50]; thresh=0.5*m_roi_pwr, initial_tfm=moving_pre_rigid, kwargs...)
    @show moving_pre_tfm

    #################### Moving trip to moving full
    print("Registering the moving full image to the fixed strip\n")
    mxshift_xcam = [50; round(Int, size(moving_roi, 2)+50)]
    f_roi_pwr = sum(abs2.(fixed_roi_pad[.!(isnan.(fixed_roi_pad))]))
    moving_post_rigid, post_rigid_mm = qd_rigid(fixed_roi_pad, moving_full, mxshift_xcam, [pi/60], [pi/10000]; thresh=0.5*f_roi_pwr, kwargs...)
    @show moving_post_rigid

    #refine it further, allowing full affine
    print("Found rigid transform, now optimizing allowing general affine transforms\n")
    moving_post_tfm, post_mm = qd_affine(fixed_roi_pad, moving_full, [50;50]; thresh=0.5*f_roi_pwr, initial_tfm=moving_post_rigid, kwargs...)
    @show moving_post_tfm

    return recenter(moving_post_tfm ∘ moving_pre_tfm, center(moving_full)), post_mm
end

trunc_above(v::T,thr) where {T} = ifelse(v<thr || isnan(v), v, T(thr))

nanz(x) = ifelse(isnan(x), zero(x), x)
#green and magenta overlay
mg_overlay(img1, img2) = colorview(RGB, nanz.(img1), nanz.(img2), nanz.(img1))
#green and magenta overlay, one color per half-image
function stitched_mg_overlay(strip1::AbstractArray{T,2}, strip2::AbstractArray{T,2}, tfm, ysz_full; flip_y2=true) where {T}
    strip0 = zeros(T, size(strip1))
    simg1 = stitch(strip1, strip0, tfm; flip_y_bottom = false, ysz_full=ysz_full)
    simg2 = stitch(strip0, strip2, tfm; flip_y_bottom = flip_y2, ysz_full=ysz_full)
    return mg_overlay(simg1,simg2)
end

#By convention, we align images to camera1, which has the top portion of the split image on OCPI2
#This returns an initial guess translation mapping camera 2's image to camera1's image after camera2's image has been flipped (flip required on OCPI2)
#ysz_roi is the pixel size of the image in the vertical direction (the dimension that gets split between cameras)
function initial_share_tfm(ysz_chip::Int, ysz_roi::Int)
    @assert iseven(ysz_chip)
    @assert iseven(ysz_roi)
    if ysz_roi > ysz_chip
        error("ROI size cannot be larger than chip size")
    end
    return Translation(0, -ysz_roi)
end
#
#Makes an initial guess based on stripsz
#function full2full(fixed, moving, stripsz::Int, mxshift, mxrot, rotstep; thresh = 0.1*sum(abs2.(moving[.!(isnan.(moving))])), tfm0 = IdentityTransformation(), kwargs...)
#    @assert size(fixed) == size(moving)
#    @show tfm_guess = initial_share_tfm(size(fixed,2), stripsz) ∘ tfm0
#    #tfm, mm = register_padmatched(fixed, moving, mxshift, mxrot, rotstep; thresh=thresh, tfm0=tfm_guess, kwargs...)
#    tfm, mm = register_padmatched(fixed, moving, mxshift, mxrot, rotstep; thresh=thresh, kwargs...)
#end
#
#function register_padmatched(fixed::AbstractMatrix, moving::AbstractMatrix, mxshift, mxrot, rotstep; thresh = 0.1*sum(abs2.(moving[.!(isnan.(moving))])), tfm0 = IdentityTransformation(), kwargs...)
#    if size(fixed,1) != size(moving, 1)
#        error("Images must be equally sized in the first dimension")
#    end
#    fixed, moving = match_ypadding(fixed, moving)
#    tfm, mm = qd_rigid(copy(fixed), copy(moving), mxshift, [mxrot...], [rotstep...]; thresh=thresh, tfm0=tfm0, kwargs...)
#end
#
#function match_ypadding(img1, img2; padval=NaN)
#    if size(img1,2) > size(img2, 2)
#        img2 = ypad(img2, size(img1,2); flip_y = false, fillval=padval)
#    end
#    if size(img1,2) < size(img2, 2)
#        img1 = ypad(img1, size(img2,2); flip_y = false, fillval=padval)
#    end
#    return img1, img2
#end
