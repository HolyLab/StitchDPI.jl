using TeamedImaging, Images, Base.Test, CoordinateTransformations

img0 = zeros(Normed{UInt16, 16}, 5,5,3,3)
imgpost = zeros(5,5,3,3)
normed1 = Normed{UInt16,16}(1.0)
imgpost[3,3,:,:] = normed1 #We hope the warped array looks like this
img0[3,3,1,:] = normed1 #centered
img0[4,3,2,:] = normed1 #shifted +1 in x
img0[4,4,3,:] = normed1 #shifted +1 in x and y
tfms = [IdentityTransformation(); Translation(1,0); Translation(1,1)]
za  = ZWarpedArray(img0, tfms, Float64; correct_bias=false, sqrt_tfm=false);

for t = 1:size(imgpost,4)
    for z = 1:size(imgpost,3)
        slw = za[:,:,z,t]
        slp = imgpost[:,:,z,t]
        @test !isnan(slw[3,3])
        for i in eachindex(slw)
            if !isnan(slw[i])
                @test slw[i] == slp[i]
            end
        end
    end
end

#test 3D (single stack)
img0 = img0[:,:,:,1]
imgpost = imgpost[:,:,:,1]
za  = ZWarpedArray(img0, tfms, Float64; correct_bias=false, sqrt_tfm=false);

for z = 1:size(imgpost,3)
    slw = za[:,:,z]
    slp = imgpost[:,:,z]
    @test !isnan(slw[3,3])
    for i in eachindex(slw)
        if !isnan(slw[i])
            @test slw[i] == slp[i]
        end
    end
end

#@test all(TeamedImaging.correctbias(rb) .== map(x->eltype(rb)(max(x-100/(2^16), 0.0)), rb))
#
##With bias correction
#ssb = StitchedSeries(imgtop, imgbottom, tfm; correct_bias=true, sqrt_tfm=false, flipy_bottom=true);
#for t = 1:size(imgtop,4)
#    for z = 1:size(imgtop,3)
#        @test all(ssb[:,1:size(imgtop,2),z,t] .== TeamedImaging.correctbias(rt))
#        @test all(ssb[:,size(imgtop,2)+1:2*size(imgtop,2),z,t] .== TeamedImaging.correctbias(rb))
#    end
#end
#
##With bias correction and square root
#ssq = StitchedSeries(imgtop, imgbottom, tfm; correct_bias=true, sqrt_tfm=true, flipy_bottom=true);
#for t = 1:size(imgtop,4)
#    for z = 1:size(imgtop,3)
#        @test all(ssq[:,1:size(imgtop,2),z,t] .== TeamedImaging.sqrtimg(Float64.(TeamedImaging.correctbias(rt))))
#        @test all(ssq[:,size(imgtop,2)+1:2*size(imgtop,2),z,t] .== TeamedImaging.sqrtimg(Float64.(TeamedImaging.correctbias(rb))))
#    end
#end
