# -*- coding: utf-8 -*-
# © 2022 GPLv3+ PhyMMR Team
import os


def printv(msg, verbosity, reqverb=1) -> None:
    if verbosity >= reqverb:
        print(msg)


def gettempdir():
    if os.path.exists("/run/shm"):
        return "/run/shm"
    elif os.path.exists("/dev/shm"):
        return "/dev/shm"
    return None
