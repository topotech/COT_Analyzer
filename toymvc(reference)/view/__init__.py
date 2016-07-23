import Tkinter as tk

class View(tk.Toplevel):
    def __init__(self, master):
        tk.Toplevel.__init__(self, master)
        self.protocol('WM_DELETE_WINDOW', self.master.destroy)
        tk.Label(self, text='My Money').pack(side='left')
        self.moneyCtrl = tk.Entry(self, width=8)
        self.moneyCtrl.pack(side='left')

    def SetMoney(self, money):
        self.moneyCtrl.delete(0,'end')
        self.moneyCtrl.insert('end', str(money))

