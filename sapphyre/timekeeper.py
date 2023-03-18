# -*- coding: utf-8 -*-
# © 2022 GPLv3+ Sapphyre Team
"""Simple helper class to measure the time a piece of code takes to execute.

It can either work with 2 time point, or a series of them.

Example:

    from timekeeper import TimeKeeper
    tike = TimeKeeper(mode=KeeperMode.DIRECT)
    ...
    print("Time taken: {} seconds".format(tike.differential())

To take multiple point, in the use case of a generator.

Example:

    from timekeeper import TimeKeeper
    tike = TimeKeeper(mode=KeeperMode.SUM)
    ...
    for something in thatthing:
        ...
        tike.timer1("now")
    print("Loop has taken {} seconds".format(tike.differential())
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
        self._t2 = None
        if auto_t1:
            self.timer1("now")
            self._t2 = self._t1

    @property
    def time1(self):
        if self.mode is KeeperMode.SUM:
            return sum(self._t1[1::])
        return self._t1

    @time1.setter
    def time1(self, t1):
        self.timer1(t1)

    @time1.deleter
    def time1(self):
        del self._t1

    def timer1(self, t1):
        if t1 == "now":
            if self.mode is KeeperMode.SUM:
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

    def lap(self):
        t2 = self._t2
        self._t2 = time()
        return self._t2 - t2

    def differential(self, t2=None):
        if self.mode is KeeperMode.SUM:
            return self.time1
        if not t2:
            return time() - self.time1
        if not isinstance(t2, float):
            raise KeeperBadTime(f"Second timestamp must be a float, got {t2} instead.")
        if not t2 >= self.time1:
            raise KeeperBadTime("Second timestamp is lower than the first.")
        return t2 - self.time1

    def format(self, string):
        """Format and return a string containing {total}."""
        return string.format(total=self.differential())

    def format2(self, string):
        """format and return a string containing {total} and {lap}"""
        return string.format(total=self.differential(), lap=self.lap())


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
    print(f"Time elapsed {tike.differential()} seconds.")
    print()
    print("Mode: sum")
    tike = TimeKeeper(KeeperMode.SUM)
    tiks = 0
    while tiks < 10:
        tike.timer1("now")
        tiks += 1
        print(f"{tiks}...")
        sleep(1)
        tike.timer1("now")

    print(f"Time elapsed {tike.differential()} seconds.")
