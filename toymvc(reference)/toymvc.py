'''
Original version by AnthonyMuss. Modified as a reference for the current project.

'''


#import Tkinter as tk ---> Imported by the view module
from controller import *


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    app = Controller(root)
    root.mainloop()