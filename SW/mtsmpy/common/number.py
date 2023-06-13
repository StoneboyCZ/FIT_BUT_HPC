import numpy as np
from mpmath import mp
from enum import Enum

class Type(Enum):
    NORMAL = 0
    VPA = 1

class Number:

    def __init__(self, v,d) -> None:
        if d > 15:
            self.val = mp.mpf(v)
            self.mode = Type.VPA 
        else:
            self.val = v
            self.mode = Type.NORMAL

    def __add__(self, o):
        return self.val + o.val 

    def __sub__(self, o):
        return self.val - o.val

    def __mult__(self, o):
        return self.val * o.val

    def __div__(self, o):
        if self.mode == Type.NORMAL:
            return self.val/o.val
        else:
            return mp.fdiv(self.val,o.val)
    
    def __str__(self) -> str:
        return str(self.val)