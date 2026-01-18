import heapq
import numpy as np


def random_select(indexes, n_sample):
    return np.random.choice(range(len(indexes)), n_sample, replace=False)

def random_select_array(indexes, n_sample):
    select_indexes = np.random.choice(range(len(indexes)), n_sample, replace=False)
    return list(map(indexes.__getitem__, list(select_indexes)))

def mu_law_encoding(data, mu):
    return np.sign(data) * np.log(1 + mu * np.abs(data)) / np.log(mu + 1)

def mu_law_decoding(data, mu):
    return np.sign(data) * (np.exp(np.abs(data) * np.log(mu + 1)) - 1) / mu

def mergesort(list_of_lists, key=None):
    heap = []
    for i, itr in enumerate(iter(pl) for pl in list_of_lists):
        try:
            item = next(itr)
            toadd = (key(item), i, item, itr) if key else (item, i, itr)
            heap.append(toadd)
        except StopIteration:
            pass
    heapq.heapify(heap)

    if key:
        while heap:
            _, idx, item, itr = heap[0]
            yield item, itr
            try:
                item = next(itr)
                heapq.heapreplace(heap, (key(item), idx, item, itr))
            except StopIteration:
                heapq.heappop(heap)

    else:
        while heap:
            item, idx, itr = heap[0]
            yield item, itr
            try:
                heapq.heapreplace(heap, (next(itr), idx, itr))
            except StopIteration:
                heapq.heappop(heap)
