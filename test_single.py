from src import fbat
f = fbat.FBAT()
f.load("data/largeset/camp.fam", "data/largeset/camp.tped")
f.setOffset(0)
f.single()

#i = 4
#s = fbat.SingleTest(f.markers[i], f.tped[i], map(lambda x: f.applyOffset(x), f.phenotypes), f.famidx, f.childidx, f.paridx)
#s.computeS()
#s.computeEofXandV()
#s.computeUandV()
#s.printTest()



