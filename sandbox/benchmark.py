import sys; sys.path.append("..")

import random
import time

import dedaLES


def timerfunc(func):
    """
    A timer decorator
    """
    def function_timer(*args, **kwargs):
        """
        A nested function for timing other functions
        """
        start = time.time()
        value = func(*args, **kwargs)
        end = time.time()
        runtime = end - start
        msg = "The runtime for {func} took {time} seconds to complete"
        print(msg.format(func=func.__name__,
                         time=runtime))
        return value
    return function_timer

model_dns = dedaLES.build_rayleigh_benard_benchmark(64, 64, 16, closure=None)
model_constsmag = dedaLES.build_rayleigh_benard_benchmark(64, 64, 16, closure=dedaLES.ConstantSmagorinsky())

@timerfunc
def dns_benchmark():
    dedaLES.run_rayleigh_benard_benchmark(model_dns)

@timerfunc
def const_smag_benchmark():
    dedaLES.run_rayleigh_benard_benchmark(model_constsmag)
    
if __name__ == '__main__':
    dns_benchmark()
    const_smag_benchmark()
