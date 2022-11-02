# -*- coding: utf-8 -*-
# © 2022 GPLv3+ bicobus <bicobus@keemail.me>
"""Simple helper class to measure the time a piece of code takes to execute.

It can either work with 2 time point, or a series of them.

Example:

    from timekeeper import TimeKeeper
    tike = TimeKeeper(mode=KeeperMode.DIRECT)
    ...
    print("Time taken: {} seconds".format(tike.differencial())

To take multiple point, in the use case of a generator.

Example:

    from timekeeper import TimeKeeper
    tike = TimeKeeper(mode=KeeperMode.SUM)
    ...
    for something in thatthing:
        ...
        tike.time1 = "now"
    print("Loop has taken {} seconds".format(tike.differencial())
"""
from enum import IntFlag, auto
from time import time


class KeeperBadTime(Exception):
    pass


class KeeperBadMode(Exception):
    pass


class KeeperMode(IntFlag):
    SUM = auto()
    DIRECT = auto()


class TimeKeeper:
    def __init__(self, mode: KeeperMode, auto_t1=True):
        if not isinstance(mode, IntFlag) or mode not in KeeperMode:
            raise KeeperBadMode()
        self.mode = mode
        self._t1 = None if mode is KeeperMode.DIRECT else []
        if auto_t1:
            self.time1 = "now"

    @property
    def time1(self):
        if self.mode is KeeperMode.SUM:
            return sum(self._t1[1::])
        return self._t1

    @time1.setter
    def time1(self, t1):
        issum = self.mode is KeeperMode.SUM
        if t1 == "now":
            if issum:
                try:
                    t = time()
                    self._t1.append(t - self._t1[0])
                    self._t1[0] = t
                except IndexError:
                    self._t1.append(time())
                return
            self._t1 = time()
        elif isinstance(t1, float):
            self._t1 = t1
        else:
            raise KeeperBadTime(f"Unable to understand value '{t1}'.")

    @time1.deleter
    def time1(self):
        del self._t1

    def differencial(self, t2=None):
        if self.mode is KeeperMode.SUM:
            return self.time1
        if not t2:
            return time() - self.time1
        if not isinstance(t2, float):
            raise KeeperBadTime(f"Second timestamp must be a float, got {t2} instead.")
        if not t2 >= self.time1:
            raise KeeperBadTime("Second timestamp is lower than the first.")
        return t2 - self.time1


if __name__ == "__main__":
    from time import sleep
    tike = TimeKeeper(KeeperMode.DIRECT)
    tiks = 0
    print("Ticking to 10.")
    print("Mode: direct")
    while tiks < 10:
        sleep(1)
        tiks += 1
        print(f"{tiks}...")
    print(f"Time elapsed {tike.differencial()} seconds.")
    print()
    print("Mode: sum")
    tike = TimeKeeper(KeeperMode.SUM)
    tiks = 0
    while tiks < 10:
        tike.time1 = "now"
        tiks += 1
        print(f"{tiks}...")
        sleep(1)
        tike.time1 = "now"

    print(f"Time elapsed {tike.differencial()} seconds.")
