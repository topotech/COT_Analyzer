# *-* coding: utf-8 *-*

"""
COTlib
======

Provides three classes needed to build a reaction network.

"""



class species:
    myReactions = set()
    def __init__(self, name):
        if isinstance(name,str):
            self.name = name
        else:
            raise TypeError


class renet:
    X = set()
    R_x = set()

class reaction:
    MyOwners = set()