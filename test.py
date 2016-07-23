# *-* coding: utf-8 *-*

import cotlib


c = cotlib.species("H2O")

print c.myReactions.__len__()

c.myReactions.add("C")
c.myReactions.add("D")

print c.myReactions
print c.myReactions.__len__()


