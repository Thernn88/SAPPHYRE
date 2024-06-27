### Clustering code

from functools import cached_property
from sapphyre_tools import get_overlap

class cluster_rec:
    def __init__(self, node, start, end, seq_data_coords, strand) -> None:
        self.node = node
        self.start = start
        self.end = end
        self.seq_data_coords = seq_data_coords
        self.strand = strand
    
    @cached_property
    def get_ids(self): 
        return node_to_ids(self.node)
    
    def get_first_id(self):
        return min(self.get_ids)

def determine_direction(start, end, current_start, current_end, current_direction, debug = False):
    this_direction = None
    if start == current_start and end == current_end:
        this_direction =  current_direction
    else:
        if start == current_start or current_direction == "reverse":
            if end >= current_end:
                this_direction = "forward"
            else:
                this_direction = "reverse"
        else:
            if start >= current_start:
                this_direction = "forward"
            else:
                this_direction = "reverse"

    if current_direction == "bi" or this_direction == "bi" or this_direction == current_direction:
        return this_direction
    return None

def quick_rec(header, frame, seq, start, end):
    data_cols = {i for i, let in enumerate(seq[start:end], start) if let != "-"}
    strand = None if frame is None else "-" if frame < 0 else "+"
    return cluster_rec(header, start, end, data_cols, strand)

def node_to_ids(node):
    if "NODE_" in node:
        node = node.replace("NODE_","")
    return list(map(lambda x: int(x.split("_")[0]), node.split("&&")))
   
def get_min_distance(a: set, b: set) -> int:
    """Get the minimum distance between two sets."""
    min_distance = float("inf")
    for a_obj in a:
        for b_obj in b:
            distance = abs(a_obj - b_obj)
            if distance < min_distance:
                min_distance = distance
    return min_distance
    
    
def within_distance(a: set, b: set, distance: int) -> bool:
    """Check if any object from set a is within distance of set b."""
    for a_obj in a:
        for b_obj in b:
            if abs(a_obj - b_obj) <= distance:
                return True
    return False


def finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage):

    if len(current_cluster) >= 2:
        cluster_data_cols = set()
        for rec in current_cluster:
            cluster_data_cols.update(rec.seq_data_coords)
            
        cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)
        clusters.append((min(current_indices), max(current_indices), cluster_coverage))
    elif len(current_cluster) == 1:
        cluster_rec = current_cluster[0]
        cluster_coverage = len(cluster_rec.seq_data_coords.intersection(ref_coords)) / len(ref_coords)
        
        if cluster_coverage > req_seq_coverage:
            clusters.append((min(current_indices), max(current_indices), cluster_coverage))
        else:
            kicks.add(cluster_rec.node)


def cluster_ids(ids, max_id_distance, max_gap, ref_coords, req_seq_coverage = 0.5):
    clusters = []
    kicks = set()
    ids.sort(key=lambda x: x.get_first_id())
    current_cluster = None
    
    for x, rec in enumerate(ids):
        
        if current_cluster is None:
            current_cluster = [rec]
            current_indices = set(rec.get_ids)
            current_direction = "bi"
        else:   
            passed = False
            cluster_distance = get_min_distance(current_indices, rec.get_ids)
            passed_id_distance = None
            passed_distance = None
            if rec.strand == current_cluster[-1].strand:
                if cluster_distance <= max_id_distance:
                    current_rec = current_cluster[-1]
                    cluster_overlap = get_overlap(rec.start, rec.end, current_rec.start, current_rec.end, -(max_gap)+1)
                    
                    if cluster_overlap:
                        cluster_pos_distance = min(
                            abs(rec.start - current_rec.start), 
                            abs(rec.start - current_rec.end), 
                            abs(rec.end - current_rec.start), 
                            abs(rec.end - current_rec.end)
                        )

                        this_direction = determine_direction(rec.start, rec.end, current_rec.start, current_rec.end, current_direction)
                        if this_direction:
                            passed_id_distance = cluster_distance
                            passed_distance = cluster_pos_distance
                            passed_direction = this_direction
                            passed = True
                
            if passed:
                # not last id
                if x != len(ids) - 1:
                    next_rec = ids[x + 1]
                    if within_distance(rec.get_ids, next_rec.get_ids, max_id_distance):
                        next_overlap = get_overlap(rec.start, rec.end, next_rec.start, next_rec.end, -(max_gap)+1)
                        if next_overlap:
                            next_distance = get_min_distance(rec.get_ids, next_rec.get_ids)
                                
                            next_amount = min(
                                abs(rec.start - next_rec.start), 
                                abs(rec.start - next_rec.end), 
                                abs(rec.end - next_rec.start), 
                                abs(rec.end - next_rec.end)
                            )
                            
                            next_direction = determine_direction(next_rec.start, next_rec.end, rec.start, rec.end, current_direction)


                            if next_direction and passed_direction and next_direction != passed_direction and next_distance < passed_id_distance and next_amount < passed_distance:
                                passed = False
                    
            if passed:
                current_cluster.append(rec)
                current_indices.update(rec.get_ids)
                if passed_direction != "bi":
                    current_direction = passed_direction
            else:
                if cluster_distance == 0:
                    new_cluster = [rec]
                    new_indices = set(rec.get_ids)
                    kick_occured = True
                    while kick_occured:
                        kick_occured = False
                        for rec_in in current_cluster:
                            if get_min_distance(new_indices, rec_in.get_ids) == 0:
                                current_cluster.remove(rec_in)
                                current_indices.difference_update(rec_in.get_ids)
                                new_cluster.append(rec_in)
                                new_indices.update(rec_in.get_ids)
                                kick_occured = True
                            
                    if current_cluster:
                        finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage)
                    
                    current_cluster = new_cluster
                    current_indices = new_indices
                    current_direction = determine_direction(rec.start, rec.end, current_rec.start, current_rec.end, "bi")
                else:
                    finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage)
                    current_cluster = [rec]
                    current_indices = set(rec.get_ids)
                    current_direction = "bi"
    
    if current_cluster:
        finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage)
                
    return clusters, kicks

###