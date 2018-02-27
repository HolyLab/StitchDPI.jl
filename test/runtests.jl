using TeamedImaging, Base.Test

#nanplus
a = ones(10,10)
b = deepcopy(a)
a[:,1] = NaN
c = TeamedImaging.nanplus(a,b)
@test all(c[:,1] .== 1.0)
@test all(c[:,2:end] .== 2.0)

#crop_vals
cr = TeamedImaging.crop_vals(a, NaN)
@test size(cr,2) == size(a,2) - 1
@test size(cr,1) == size(a,1)
@test all(cr.==a[:,2:end])
a[8:end,:] = NaN
cr = TeamedImaging.crop_vals(a, NaN)
@test size(cr,2) == size(a,2) - 1
@test size(cr,1) == size(a,1) - 3
@test all(cr.==a[1:7,2:end])

#ypad
r = rand(10,10)
rp = TeamedImaging.ypad(r, 12; flip_y = false, fillval=NaN)
@test all(isnan.(rp[:,1]))
@test all(isnan.(rp[:,12]))
@test all(rp[:,2:11].==r)

#flipy
rpf = TeamedImaging.ypad(r, 12; flip_y = true, fillval=NaN)
@test all(reverse(rpf[1,2:11]).==rp[1,2:11])
@test all(reverse(rpf[10,2:11]).==rp[10,2:11])

#initial_share_tfm
ysz_chip = 1000
ysz_roi = 500
tfm0 = TeamedImaging.initial_share_tfm(ysz_chip, ysz_roi)

#fake split
ft, rt, fb, rb = TeamedImaging.fake_split(8, 4; frac_overlap = 0.0, xsz = 4)
@test all(ft[:,3:6].==rt)
@test all(fb[:,3:6].==rb)
@test all(ft[:,7:8].==rb[:,1:2])
@test all(fb[:,1:2].==rt[:,3:4])
#one pixel overlap
ft, rt, fb, rb = TeamedImaging.fake_split(8, 4; frac_overlap = 0.25, xsz = 4)
@test all(ft[:,3:6].==rt)
@test all(fb[:,3:6].==rb)
@test all(rt[:,4].==rb[:,1])

#stitching
ft, rt, fb, rb = TeamedImaging.fake_split(80, 40; frac_overlap = 0.0, xsz = 40)
tfm = TeamedImaging.initial_share_tfm(size(ft,2), size(rt,2))
#imgst = stitch(rt, TeamedImaging.flipy(rb), tfm, size(ft,2); flip_y_img2 = true)
imgst = stitch(rt, TeamedImaging.flipy(rb), tfm; flip_y_img2 = true)
@test all(imgst[:,1:40].==rt)
@test all(imgst[:,41:80].==rb)

#stitch function generation
tfm = stitch_tfm(ft, TeamedImaging.flipy(fb), rt, TeamedImaging.flipy(rb); maxevals=500)
#simg = stitch(rt, TeamedImaging.flipy(rb), tfm, size(ft,2); flip_y_img2 = true)
simg = stitch(rt, TeamedImaging.flipy(rb), tfm; flip_y_img2 = true)
@test all(simg.==imgst)

#with ten pixels of overlap
ft, rt, fb, rb = TeamedImaging.fake_split(80, 40; frac_overlap = 0.25, xsz = 40)
tfm = stitch_tfm(ft, TeamedImaging.flipy(fb), rt, TeamedImaging.flipy(rb); maxevals=500)
#simg = stitch(rt, TeamedImaging.flipy(rb), tfm, size(ft,2); flip_y_img2 = true)
simg = stitch(rt, TeamedImaging.flipy(rb), tfm; flip_y_img2 = true)
@test all(simg[:,31:40].==(rt[:,31:40].*2))


##bug reports
#using Images, CoordinateTransformations, Interpolations
#img = rand(4,8)
#img[:,1:2] = NaN
#img[:,7:8] = NaN
#imgw = warp(img, IdentityTransformation(), indices(img), Linear(), Flat())
#@show img
##Notice that column 6 gets set to NaN
#@show parent(imgw)
#
#img[:,1:2] = NaN
#img[:,7:8] = NaN
#using Interpolations
#img = rand(5,5)
#img[3,3] = NaN
#itp = interpolate(img, BSpline(Linear()), OnGrid())
#img2 = itp[indices(img)...]
#@show img
#@show img2

#Notice that column 6 gets set to NaN
#itp = interpolate(img, BSpline(Linear()), OnCell())
