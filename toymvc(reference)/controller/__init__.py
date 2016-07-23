from model import *
from view import *


class ChangerWidget(tk.Toplevel):
    def __init__(self, master):
        tk.Toplevel.__init__(self, master)
        self.addButton = tk.Button(self, text='Add', width=8)
        self.addButton.pack(side='left')
        self.removeButton = tk.Button(self, text='Remove', width=8)
        self.removeButton.pack(side='left')



class Controller:
    def __init__(self, root):
        self.model = Model()
        self.model.myMoney.addCallback(self.MoneyChanged)
        self.view1 = View(root)
        self.view2 = ChangerWidget(self.view1)
        self.view2.addButton.config(command=self.AddMoney)
        self.view2.removeButton.config(command=self.RemoveMoney)
        self.MoneyChanged(self.model.myMoney.get())

    def AddMoney(self):
        self.model.addMoney(10)

    def RemoveMoney(self):
        self.model.removeMoney(10)

    def MoneyChanged(self, money):
        self.view1.SetMoney(money)