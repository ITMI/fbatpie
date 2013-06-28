from src import fbat
f = fbat.FBAT()
f.load("data/largeset/camp.tfam", "data/largeset/camp.tped")
f.setOffset(0.0)
f.rare("data/largeset/camp_regions.txt",
       "data/largeset/camp_freqs.txt",
       "data")
