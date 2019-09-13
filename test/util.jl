using StitchDPI, Test

#nanplus
a = ones(10,10)
b = deepcopy(a)
a[:,1] .= NaN
c = StitchDPI.nanplus(a,b)
@test all(c[:,1] .== 1.0)
@test all(c[:,2:end] .== 2.0)

#crop_vals
cr = StitchDPI.crop_vals(a, NaN)
@test size(cr,2) == size(a,2) - 1
@test size(cr,1) == size(a,1)
@test all(cr.==a[:,2:end])
a[8:end,:] .= NaN
cr = StitchDPI.crop_vals(a, NaN)
@test size(cr,2) == size(a,2) - 1
@test size(cr,1) == size(a,1) - 3
@test all(cr.==a[1:7,2:end])

#ypad
r = rand(10,10)
rp = StitchDPI.ypad(r, 12; flip_y = false, fillval=NaN)
@test all(isnan.(rp[:,1]))
@test all(isnan.(rp[:,12]))
@test all(rp[:,2:11].==r)

#flipy
rpf = StitchDPI.ypad(r, 12; flip_y = true, fillval=NaN)
@test all(reverse(rpf[1,2:11]).==rp[1,2:11])
@test all(reverse(rpf[10,2:11]).==rp[10,2:11])
