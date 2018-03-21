#shared image alignment procedure on the microscope
#1. Install KEM.  Make sure the KEM is positioned so that the image stays crisp on camera2 (may want to compare with a dichroic to make sure)
#2. Set camera ROIs as desired (same size for both cameras)
#3. Move camera1 up so that the KEM line corresponds with the bottom side of the ROI (larger ROIs will require moving the camera farther)
#4. Move camera2 left so that the KEM line corresponds with the bottom side of the ROI
#5. Note the ROI settings and change the ROI to full chip temporarily
#6. Install 50/50 dichroic instead of KEM
#7. Take a snapshot of beads with both cameras
#8. Reset ROI settings, replace KEM mirror, and take another pair of snapshots.  Make sure the beads don't move during this sequence!

#shared image alignment computational procedure
#1.  Compute transform to align camera2's flipped padded half-image (step 8) with its flipped full image (step 7) (should be almost perfect already)
#2.  Compute transform to align camera2's flipped full image (step 7) with camera1's full image (step 7) (will require a downward shift roughly equal to ROI ysz)
#3.  Compute transform to align camera1's full image (step 7) with camera1's padded half-image image (step 8)
#4.  For each camera2 half-image:
    # -flip y
    # -pad vertically to full chip size
    # -apply transform from #1 
    # -apply transform from #2
    # -apply transform from #3 (Actually all 3 transforms should be composed into one and applied just once)
    # -nansum with camera1's corresponding padded half-image
    # -find minimal bounding box that contains all non-nan pixels and return image cropped to that region

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

#Stitches img_top_fixed and img_bottom_moving using a tranformation tfm that, when applied to img_bottom_moving, aligns it with img_top_fixed.  Since img_top_fixed and img_bottom_moving are expected to have little-to-no overlap
#they are padded before applying the transformation.  The transformed img_bottom_moving is added to img_top_fixed and the interesting (non-NaN) region of the summed image is returned.
#function stitch(img_top_fixed, img_bottom_moving, tfm, y_sz; flip_y_bottom = true)
function stitch(img_top_fixed::AbstractArray{T,2}, img_bottom_moving::AbstractArray{T,2}, tfm, out_type::DataType; flip_y_bottom = true) where {T}
    if T!=out_type
        img_top_fixed = out_type.(img_top_fixed)
        img_bottom_moving = out_type.(img_bottom_moving)
    end
    return stitch(img_top_fixed, img_bottom_moving, tfm; flip_y_bottom=flip_y_bottom)
end

function stitch(img_top_fixed::AbstractArray{T,2}, img_bottom_moving::AbstractArray{T,2}, tfm; flip_y_bottom = true) where {T}
    yroi = size(img_top_fixed,2)
    @assert size(img_bottom_moving,2) == yroi
    yextra = ceil(Int, abs(tfm.v[2]))
    #Note: should probably recenter tfm because it was centered on a full-chip image when found with stitch_tfm()
    img_bottom_moving_pre = ypad(img_bottom_moving, 3*yroi; pad_side=:both, flip_y=flip_y_bottom, fillval=0.0) #fillvals should be NaN, but currently interpolations doesn't handle it well (imputes NaNs at gridpoints)
    img_top_fixed = ypad(img_top_fixed, 3*yroi; pad_side=:both, flip_y=false, fillval=0.0)
    img_bottom_moving = warp_and_resample(img_bottom_moving_pre, tfm)
    imgsummed = nanplus(img_top_fixed,img_bottom_moving)
    return view(imgsummed, :, (yroi+1):3*yroi)
end

#returns a transform that stitches a pair of corresponding shared images into a single image when passed to the stitch function
#Assumes that moving_roi contains the bottom half of the image and that it needs to be flipped vertically to align with fixed_full and fixed_roi.
#We can't assume that fixed_full is perfectly aligned with the half-image (fixed_roi) from its camera. That's because the removal of the dichroic results in a 
#lateral shift of the image, so we must compensate for that.
function stitch_tfm(fixed_full, moving_full, fixed_roi, moving_roi; kwargs...)
    @assert size(fixed_full) == size(moving_full)
    moving_full = flipy(moving_full)
    fixed_roi_pad = ypad(fixed_roi, size(fixed_full,2); flip_y = false, fillval=NaN)
    moving_roi_pad = ypad(moving_roi, size(moving_full,2); flip_y = true, fillval=NaN)
    print("Registering the moving strip to the moving full image\n")
    mxshift_moving = [10;10]
    moving_pre_tfm, pre_mm = register_padmatched(moving_full, moving_roi_pad, mxshift_moving, [pi/2000;], [pi/1000]; thresh=0.4*sum(abs2.(moving_roi_pad[.!(isnan.(moving_roi_pad))])), kwargs...)
    mxshift_xcam = [50; ceil(Int, size(moving_roi, 2)+100)]
    print("Registering the moving full image to the fixed full image\n")
    moving_post_tfm, post_mm = full2full(fixed_full, moving_roi_pad, size(moving_roi, 2), mxshift_xcam, [pi/10], [pi/10000]; tfm0=moving_pre_tfm, thresh=0.4*sum(abs2.(moving_roi_pad[.!(isnan.(moving_roi_pad))])), kwargs...)
    mxshift_fixed = [100;100]
    print("Registering the fixed full image to the fixed strip\n")
    fixed_self_tfm, mm_self1 = register_padmatched(fixed_roi_pad, fixed_full, mxshift_fixed, [pi/2000], [pi/1000]; thresh=0.4*sum(abs2.(fixed_roi_pad[.!(isnan.(fixed_roi_pad))])), kwargs...)
    return fixed_self_tfm ∘ moving_post_tfm
end

#Makes an initial guess based on stripsz
function full2full(fixed, moving, stripsz::Int, mxshift, mxrot, rotstep; thresh = 0.1*sum(abs2.(moving[.!(isnan.(moving))])), tfm0 = IdentityTransformation(), kwargs...)
    @assert size(fixed) == size(moving)
    tfm_guess = tfm0 ∘ initial_share_tfm(size(fixed,2), stripsz)
    tfm, mm = register_padmatched(fixed, moving, mxshift, mxrot, rotstep; thresh=thresh, tfm0=tfm_guess, kwargs...)
end

function register_padmatched(fixed::AbstractMatrix, moving::AbstractMatrix, mxshift, mxrot, rotstep; thresh = 0.1*sum(abs2.(moving[.!(isnan.(moving))])), tfm0 = IdentityTransformation(), kwargs...)
    if size(fixed,1) != size(moving, 1)
        error("Images must be equally sized in the first dimension")
    end
    fixed, moving = match_ypadding(fixed, moving)
    tfm, mm = qd_rigid(fixed, moving, mxshift, [mxrot...], [rotstep...]; thresh=thresh, tfm0=tfm0, kwargs...)
end

function match_ypadding(img1, img2; padval=NaN)
    if size(img1,2) > size(img2, 2)
        img2 = ypad(img2, size(img1,2); flip_y = false, fillval=padval)
    end
    if size(img1,2) < size(img2, 2)
        img1 = ypad(img1, size(img2,2); flip_y = false, fillval=padval)
    end
    return img1, img2
end
