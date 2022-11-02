from wrap_rocks import wrap_rocks


def create_pointer(key, path):
    global rockdb_pointers
    if not rockdb_pointers.get(path):
        rockdb_pointers[key] = wrap_rocks.RocksDB(path)


def get_rock(key):
    global rockdb_pointers
    return rockdb_pointers[key]


rockdb_pointers = {}
