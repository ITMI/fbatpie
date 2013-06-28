import fbat
f = fbat.FBAT()
f.load("large_test_set/camp.tfam", "large_test_set/camp.tped")
f.setOffset(0.5)
f.fbat()



