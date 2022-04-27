using StitchDPI, ImageBase, Test

ft, rt, fb, rb = fake_split(8, 4; frac_overlap = 0.0, xsz = 4)
tfm = StitchDPI.initial_guess_tfm(size(ft,2), size(rt,2))
rt = Normed{UInt16,16}.(rt)
rb = Normed{UInt16,16}.(rb)

imgtop = zeros(Normed{UInt16,16}, 4,4,3,3)
imgbottom = zeros(Normed{UInt16, 16}, 4,4,3,3)
for t = 1:size(imgtop,4)
    for z = 1:size(imgtop,3)
        imgtop[:,:,z,t] = rt
        imgbottom[:,:,z,t] = StitchDPI.flipy(rb)
    end
end

ss = StitchedSeries(imgtop, imgbottom, tfm; correct_bias=false, sqrt_tfm=false, flip_y_bottom=true);
for t = 1:size(imgtop,4)
    for z = 1:size(imgtop,3)
        @test all(ss[:,1:size(imgtop,2),z,t] .== rt)
        @test all(ss[:,size(imgtop,2)+1:2*size(imgtop,2),z,t].== rb)
    end
end

@test all(StitchDPI.correctbias(rb) .== map(x->eltype(rb)(max(x-100/(2^16), 0.0)), rb))

#With bias correction
ssb = StitchedSeries(imgtop, imgbottom, tfm; correct_bias=true, sqrt_tfm=false, flip_y_bottom=true);
for t = 1:size(imgtop,4)
    for z = 1:size(imgtop,3)
        @test all(ssb[:,1:size(imgtop,2),z,t] .== StitchDPI.correctbias(rt))
        @test all(ssb[:,size(imgtop,2)+1:2*size(imgtop,2),z,t] .== StitchDPI.correctbias(rb))
    end
end

#With bias correction and square root
ssq = StitchedSeries(imgtop, imgbottom, tfm; correct_bias=true, sqrt_tfm=true, flip_y_bottom=true);
for t = 1:size(imgtop,4)
    for z = 1:size(imgtop,3)
        @test all(ssq[:,1:size(imgtop,2),z,t] .== StitchDPI.sqrtimg(Float64.(StitchDPI.correctbias(rt))))
        @test all(ssq[:,size(imgtop,2)+1:2*size(imgtop,2),z,t] .== StitchDPI.sqrtimg(Float64.(StitchDPI.correctbias(rb))))
    end
end
