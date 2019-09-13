using StitchDPI, Images, Test

img = rand(3,3,3,4)
multiproc_write("test.nhdr", img, Float32)

img2 = load("test.nhdr")

@test all(img2.==Float32.(img))
