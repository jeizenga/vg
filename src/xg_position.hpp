#ifndef VG_XG_POS_HPP_INCLUDED
#define VG_XG_POS_HPP_INCLUDED

#include "vg.pb.h"
#include "types.hpp"
#include "xg.hpp"
#include "lru_cache.h"
#include "utility.hpp"
#include "json2pb.h"
#include <gcsa/gcsa.h>
#include <iostream>

/** \file 
 * Functions for working with cached Positions and `pos_t`s.
 */

namespace vg {

using namespace std;

// xg/position traversal helpers with caching
// used by the Sampler and by the Mapper
string xg_node_sequence(id_t id, const xg::XG* xgidx);
/// Get the length of a Node from an xg::XG index, with cacheing of deserialized nodes.
size_t xg_node_length(id_t id, const xg::XG* xgidx);
/// Get the node start position in the sequence vector
int64_t xg_node_start(id_t id, const xg::XG* xgidx);
/// Get the character at a position in an xg::XG index, with cacheing of deserialized nodes.
char xg_pos_char(pos_t pos, const xg::XG* xgidx);
/// Get the characters at positions after the given position from an xg::XG index, with cacheing of deserialized nodes.
map<pos_t, char> xg_next_pos_chars(pos_t pos, const xg::XG* xgidx);
set<pos_t> xg_next_pos(pos_t pos, bool whole_node, const xg::XG* xgidx);
int64_t xg_distance(pos_t pos1, pos_t pos2, int64_t maximum, const xg::XG* xgidx);
set<pos_t> xg_positions_bp_from(pos_t pos, int64_t distance, bool rev, const xg::XG* xgidx);
//void xg_graph_context(VG& graph, const pos_t& pos, int length, xg::XG* xgidx);
Node xg_node(id_t id, const xg::XG* xgidx);
vector<Edge> xg_edges_on_start(id_t id, const xg::XG* xgidx);
vector<Edge> xg_edges_on_end(id_t id, const xg::XG* xgidx);
map<string, vector<pair<size_t, bool> > > xg_alignment_path_offsets(const Alignment& aln, bool just_min, bool nearby, const xg::XG* xgidx);
void xg_annotate_with_initial_path_positions(Alignment& aln, bool just_min, bool nearby, const xg::XG* xgidx);

}

#endif
