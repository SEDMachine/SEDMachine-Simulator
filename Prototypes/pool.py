
from multiprocessing import Pool

def f(x):
    return x*x*x


po = Pool()

po.map_async(f,range(100))

po.close()
po.join()