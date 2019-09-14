using StitchDPI, Test

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
