from wrap_rocks import RocksDB
from os import path
from .utils import printv

def recover_ref_db():
    print("TODO")
    return True

def recover_taxa(nt_db, out):
    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
        
    with open(out, "wb") as f:
        for i in recipe:
            f.write(nt_db.get_bytes(f"ntbatch:{i}"))
    print("Done!")
    

def run_process(args, input_path):
    if path.isfile(input_path):
        printv("ERROR: Input path must be a folder.", args.verbose, 0)
        return False
    
    try:
        if path.exists(path.join(input_path, "rocksdb")):
            db = RocksDB(path.join(input_path, "rocksdb", "sequences", "nt"))
        else:
            db = RocksDB(input_path)
    except:
        printv("ERROR: Input path must be a valid RocksDB folder.", args.verbose, 0)
        return False
    
    #Try to grab a ref db data
    if db.get("getall:taxoninset") is None:
        out = path.join(input_path, path.basename(input_path))
        recover_taxa(db, out)
    else:
        recover_ref_db()

    return True

def main(args):
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    results = []
    for input_path in args.INPUT:
        results.append(run_process(args, input_path))

    return all(results)