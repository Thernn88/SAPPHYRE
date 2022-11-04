

def printv(msg, verbosity, reqverb=1) -> None:
    if verbosity >= reqverb :
        print(msg)
