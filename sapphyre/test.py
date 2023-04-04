from multiprocessing.pool import ThreadPool as Pool

# from multiprocessing import Pool


def process_line(l):
    header, ref_header_pair, frame = line.split("\t")[:3]
    if hash(header) != current_header:
        current_header = hash(header)
        ref_variation_filter = set()

    ref_gene, ref_taxon, _ = taxon_lookup[ref_header_pair]

    target_has_hit.add(ref_header_pair)

    ref_variation_key = header + ref_taxon
    if not ref_variation_key in ref_variation_filter:
        ref_variation_filter.add(ref_variation_key)
        rextaxon_count[ref_taxon] = rextaxon_count.get(ref_taxon, 0) + 1
    header_lines.setdefault(header + frame, []).append(
        line.strip() + f"\t{ref_gene}\t{ref_taxon}"
    )


def get_next_line():
    with open("test.csv", "r") as f:
        for line in f:
            yield line


f = get_next_line()

t = Pool(processes=8)

for i in f:
    t.map(process_line, (i,))
t.close()
t.join()
