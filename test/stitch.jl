using StitchDPI, LinearAlgebra, Test

#initial_share_tfm
ysz_chip = 1000
ysz_roi = 500
tfm0 = StitchDPI.initial_share_tfm(ysz_chip, ysz_roi)
tfm_tol = 1e-6

#fake split
ft, rt, fb, rb = StitchDPI.fake_split(8, 4; frac_overlap = 0.0, xsz = 4)
@test all(ft[:,3:6].==rt)
@test all(fb[:,3:6].==rb)
@test all(ft[:,7:8].==rb[:,1:2])
@test all(fb[:,1:2].==rt[:,3:4])
#one pixel overlap
ft, rt, fb, rb = StitchDPI.fake_split(8, 4; frac_overlap = 0.25, xsz = 4)
@test all(ft[:,3:6].==rt)
@test all(fb[:,3:6].==rb)
@test all(rt[:,4].==rb[:,1])

#stitching
ft, rt, fb, rb = StitchDPI.fake_split(80, 40; frac_overlap = 0.0, xsz = 40)
tfm = StitchDPI.initial_share_tfm(size(ft,2), size(rt,2))
#imgst = stitch(rt, StitchDPI.flipy(rb), tfm, size(ft,2); flip_y_bottom = true)
imgst = stitch(rt, StitchDPI.flipy(rb), tfm; flip_y_bottom = true)
@test all(imgst[:,1:40].==rt)
@test all(imgst[:,41:80].==rb)

#stitch transform generation
tfmfound, mm = stitch_tfm(StitchDPI.flipy(fb), rt, StitchDPI.flipy(rb); maxevals=500)
@test sum(abs.(tfmfound.linear - LinearAlgebra.I(2))) < tfm_tol
@test sum(abs.(tfmfound.translation - tfm.translation)) < tfm_tol
#simg = stitch(rt, StitchDPI.flipy(rb), tfmfound, size(ft,2); flip_y_bottom = true)
#simg = stitch(rt, StitchDPI.flipy(rb), tfmfound; flip_y_bottom = true)
#@test all(simg.==imgst)

#with ten pixels of overlap
ft, rt, fb, rb = StitchDPI.fake_split(80, 40; frac_overlap = 0.25, xsz = 40)
tfmfound, mm = stitch_tfm(StitchDPI.flipy(fb), rt, StitchDPI.flipy(rb); maxevals=500)
@test sum(abs.(tfmfound.linear - LinearAlgebra.I(2))) < tfm_tol
@test sum(abs.(tfmfound.translation - (tfm.translation.+[0.0; 10.0]))) < tfm_tol
#simg = stitch(rt, StitchDPI.flipy(rb), tfmfound, size(ft,2); flip_y_bottom = true)
#simg = stitch(rt, StitchDPI.flipy(rb), tfmfound; flip_y_bottom = true)
#@test all(simg[:,31:40].==(rt[:,31:40].*2))
