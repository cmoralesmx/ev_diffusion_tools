import multiprocessing
import numpy as np

#global data


def numbers(r):
    #global data
    data = np.load('data.npy', 'r+')
    for i in range(9):
        data[r][i] = r * 10 + i
    return


if __name__ == '__main__':
    #global data

    data = np.zeros([5,9])
    np.save('data.npy', data)

    jobs = []
    for i in range(5):
        p = multiprocessing.Process(target=numbers, args=(i,))
        jobs.append(p)
        p.start()

        p.join()
    data = np.load('data.npy', 'r+')
    
    print(data[0,:])
    print(data[1,:])
    print(data[2,:])
    print(data[3,:])
    print(data[4,:])
