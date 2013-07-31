from src import fbat
f = fbat.FBAT()
f.load("data/muc17test/muc17.tfam", "data/muc17test/muc17.tped")
f.setOffset(0.5)
f.rare("data/muc17test/rarerunlist.txt",  # region
       "none", # index
       "none", # freq
       "madsen") # weighting
