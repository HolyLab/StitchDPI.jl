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

#Stitches img1 and img2 using a tranformation tfm that, when applied to img2, aligns it with img1.  Since img1 and img2 are expected to have little-to-no overlap
#they are padded before applying the transformation.  The transformed img2 is added to img1 and the interesting (non-NaN) region of the summed image is returned.
#function stitch(img1, img2, tfm, y_sz; flip_y_img2 = true)
function stitch(img1::AbstractArray{T,2}, img2::AbstractArray{T,2}, tfm, out_type::DataType; flip_y_img2 = true) where {T}
    if T!=out_type
        img1 = out_type.(img1)
        img2 = out_type.(img2)
    end
    return stitch(img1, img2, tfm; flip_y_img2=flip_y_img2)
end

function stitch(img1::AbstractArray{T,2}, img2::AbstractArray{T,2}, tfm; flip_y_img2 = true) where {T}
    yroi = size(img1,2)
    @assert size(img2,2) == yroi
    yextra = ceil(Int, abs(tfm.v[2]))
    img2_pre = ypad(img2, yroi+2*yextra; pad_side=:both, flip_y=flip_y_img2, fillval=0.0) #fillvals should be NaN, but currently interpolations doesn't handle it well (imputes NaNs at gridpoints)
    img1 = ypad(img1, yroi+2*yextra; pad_side=:both, flip_y=false, fillval=0.0)
    img2 = warp(img2_pre, tfm)
    itp = interpolate(img2, BSpline(Linear()), OnGrid())
    etp = extrapolate(itp, 0.0)
    img2 = etp[indices(img2_pre)...]
    imgsummed = nanplus(img1,img2)
    return crop_vals(imgsummed, 0.0)
end

#returns a function f(img1,img2) that stitches a pair of corresponding shared images into a single image.
#Assumes that img2_roi contains the bottom half of the image and that it needs to be flipped vertically to align with img1_full
#We also assume that img1_full is perfectly aligned with the half-image from its camera.  This is valid on OCPI2 because camera1's portion of the image does not hit the KEM.  SCRATCH THAT!  We can't assume this, the removal of the dichroic results in a 
#lateral shift of the image, we must compensate for that.
function stitch_tfm(img1_full, img2_full, img1_roi, img2_roi; kwargs...)
    @assert size(img1_full) == size(img2_full)
    img2_full = flipy(img2_full)
    img2_roi_pad = ypad(img2_roi, size(img2_full,2); flip_y = true, fillval=NaN)
    img1_roi_pad = ypad(img1_roi, size(img1_full,2); flip_y = false, fillval=NaN)
    mxshift_img2 = [10;10]
    print("Registering the img2 strip to the img2 full image\n")
    img2_pre_tfm, pre_mm = qd_rigid(img2_full, img2_roi_pad, mxshift_img2, [pi/2000], [pi/1000]; thresh=0.4*sum(abs2.(img2_roi_pad[.!(isnan.(img2_roi_pad))])), kwargs...)
    post_tfm_guess = initial_share_tfm(size(img2_full,2), size(img2_roi,2)) ∘ img2_pre_tfm 
    mxshift_xcam = [50; ceil(Int, size(img2_roi, 2)+100)]
    print("Registering the img2 full image to the img1 full image\n")
    img2_post_tfm, post_mm = qd_rigid(img1_full, img2_roi_pad, mxshift_xcam, [pi/10], [pi/10000]; tfm0= post_tfm_guess, thresh=0.4*sum(abs2.(img2_roi_pad[.!(isnan.(img2_roi_pad))])), kwargs...)
    mxshift_img1 = [100;100]
    print("Registering the img1 full image to the img1 strip\n")
    img1_self_tfm, mm_self1 = qd_rigid(img1_roi_pad, img1_full, mxshift_img1, [pi/2000], [pi/1000]; thresh=0.4*sum(abs2.(img1_roi_pad[.!(isnan.(img1_roi_pad))])), kwargs...)
    return img1_self_tfm ∘ img2_post_tfm
end

#fake split images with requested split
#returns the top full/split images first, then bottom
#ysz_roi is the size of one camera's ROI (so total combined image will have dimension 2x ysz_roi)
function fake_split(ysz_chip, ysz_roi; frac_overlap=0.0, xsz=20)
    @assert iseven(ysz_chip)
    @assert iseven(ysz_roi)
    @assert 0.0 <= frac_overlap <= 1.0
    npix_extra = div(ysz_chip - ysz_roi, 2) #unused vertical roi  (symmetric above and below)
    npix_overlap = ceil(Int, frac_overlap * ysz_roi) #number of pixels overlapping between shared images (zero would mean no redundant pixels)
    chip_halfsz = div(ysz_chip,2)
    ysz_full = 2*ysz_roi - npix_overlap + 2*npix_extra
    full_img = rand(xsz, ysz_full)
    full_top = full_img[:,1:ysz_chip]
    full_bottom = full_img[:,(ysz_full-ysz_chip+1):ysz_full]
    img_top = full_top[:,(npix_extra+1):(npix_extra+ysz_roi)]
    img_bottom = full_bottom[:,(npix_extra+1):(npix_extra+ysz_roi)]
    return full_top, img_top, full_bottom, img_bottom
end
