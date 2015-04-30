cdef extern from "numerics.h" namespace "tbem":
    double from_11_to_01(double)

def call(x):
    return from_11_to_01(x)
print("Hello World")
