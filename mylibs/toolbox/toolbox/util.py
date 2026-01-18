import logging
import functools
import time
import uuid
from datetime import datetime
from itertools import repeat


def generate_uuid():
    return str(uuid.uuid4())

def generate_datetime():
    return datetime.fromtimestamp(time.time()).strftime(f'%Y%m%d_%H%M%S_%f')

def get_current_datetime():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]

def timeit(func):
    @functools.wraps(func)
    def wrapper(*args, **kargs):
        start = time.time()
        result = func(*args, **kargs)
        stop = time.time()
        print(stop - start)
        return result
    return wrapper

'''
example:
@notNoneAttrs('var1_name', 'var2_name', ...)
def func():
    pass
'''
def notNoneAttrs(*notNoneAttrs):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            for notNoneAttr in notNoneAttrs:
                if getattr(self, notNoneAttr) is None:
                    return None
            return func(self, *args, **kwargs)
        return wrapper
    return decorator

def create_simple_logger(logger_name='Test'):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    streamHandler = logging.StreamHandler()
    basic_format = '%(asctime)s [%(levelname)s] %(filename)s [line:%(lineno)d] %(message)s'
    datefmt = '%Y-%m-%d %H:%M:%S'
    streamHandler.setFormatter(logging.Formatter(fmt=basic_format, datefmt=datefmt))
    logger.addHandler(streamHandler)
    return logger

def inf_loop(data_loader):
    ''' wrapper function for endless data loader. '''
    for loader in repeat(data_loader):
        yield from loader

def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False