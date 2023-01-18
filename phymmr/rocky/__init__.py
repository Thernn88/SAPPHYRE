from wrap_rocks import wrap_rocks


def create_pointer(key, path):
    global rockdb_pointers
    if not rockdb_pointers.get(key):
        rockdb_pointers[key] = wrap_rocks.RocksDB(path)


def close_pointer(key):
    global rockdb_pointers
    if rockdb_pointers.get(key):
        rockdb_pointers.pop(key)


def get_rock(key):
    global rockdb_pointers
    return rockdb_pointers[key]


rockdb_pointers = {}
