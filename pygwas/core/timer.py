#!/usr/bin/env python
# coding: utf-8
"""
    Used for timing functionality
"""


import time
import datetime 

class Timer(object):
    def __init__(self):
        self.start()

    def start(self):
        self._start = time.clock()

    def stop(self,is_readable=False):
        self._end = time.clock()
        self._secs = self._end - self._start
	if is_readable: 
            return self.readable()
        return self._secs

    def readable(self):
        return str(datetime.timedelta(seconds = self._secs))
