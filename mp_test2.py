from multiprocessing import Process, Array
import numpy as np

def main():
    im_arr = np.array([[1,2,3],[4,6,7]])
    print('Array in main before process:',im_arr)

    shape = im_arr.shape
    size = im_arr.size
    im_arr.shape = size
    arr = Array('B', im_arr)   
    p1 = Process(target=fun, args=(arr,shape, 0))
    p2 = Process(target=fun, args=(arr,shape, 1))
    p1.start()
    p2.start()
    p1.join()
    p2.join()

    arr = np.frombuffer(arr.get_obj(), dtype=np.uint8)
    arr.shape = shape
    print('Array in main after process:',arr)

def fun(a, shape, n):
    a = np.frombuffer(a.get_obj(), dtype=np.uint8)
    a.shape = shape

    a[n][1] = (1+n)
    #a[...] = np.array([[0,1,0],[0,2,0]])
    a[n][2] = 2*(n+2)

    print('Array inside function:',a)
    a.shape = shape[0]*shape[1]

if __name__ == '__main__':
    main()