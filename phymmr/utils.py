# -*- coding: utf-8 -*-
# © 2022 GPLv3+ PhyMMR Team
import os
from threading import Thread
from queue import Queue

class ConcurrentLogger(Thread):
    def __init__(self, inq: Queue):
        super().__init__(daemon=True)
        self.inq = inq

    def run(self):
        while True:
            message, verbosity, reqv = self.inq.get()            
            printv(message, verbosity, reqv)
            self.inq.task_done()

    def __call__(self, msg: str, verbosity: int, reqv=1):
        self.inq.put((msg, verbosity, reqv))

def printv(msg, verbosity, reqverb=1) -> None:
    if verbosity >= reqverb:
        print(msg)


def gettempdir():
    if os.path.exists("/run/shm"):
        return "/run/shm"
    elif os.path.exists("/dev/shm"):
        return "/dev/shm"
    return None
