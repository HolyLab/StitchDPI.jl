using TeamedImaging, Base.Test

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

#stitch transform generation
tfmfound = stitch_tfm(ft, TeamedImaging.flipy(fb), rt, TeamedImaging.flipy(rb); maxevals=500)
@test sum(abs.(tfmfound.m - eye(2))) < 1e-7
@test sum(abs.(tfmfound.v - tfm.v)) < 1e-7
#simg = stitch(rt, TeamedImaging.flipy(rb), tfmfound, size(ft,2); flip_y_img2 = true)
#simg = stitch(rt, TeamedImaging.flipy(rb), tfmfound; flip_y_img2 = true)
#@test all(simg.==imgst)

#with ten pixels of overlap
ft, rt, fb, rb = TeamedImaging.fake_split(80, 40; frac_overlap = 0.25, xsz = 40)
tfmfound = stitch_tfm(ft, TeamedImaging.flipy(fb), rt, TeamedImaging.flipy(rb); maxevals=500)
@test sum(abs.(tfmfound.m - eye(2))) < 1e-7
@test sum(abs.(tfmfound.v - (tfm.v.+[0.0; 10.0]))) < 1e-7
#simg = stitch(rt, TeamedImaging.flipy(rb), tfmfound, size(ft,2); flip_y_img2 = true)
#simg = stitch(rt, TeamedImaging.flipy(rb), tfmfound; flip_y_img2 = true)
#@test all(simg[:,31:40].==(rt[:,31:40].*2))
