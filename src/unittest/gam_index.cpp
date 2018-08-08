///
///  \file genome_state.cpp
///
///  Unit tests for the GAMIndex which indexes seekable GAM files by node ID
///

#include <iostream>
#include "catch.hpp"
#include "../gam_index.hpp"
#include "../utility.hpp"


namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("GAMindex windowing works correctly", "[gam][gamindex]") {

    for (size_t i = 0; i < 10000000; i += 83373) {
        REQUIRE(GAMIndex::window_of_id(i) == i / 256);
    }

}

TEST_CASE("GAMindex binning works on a large number", "[gam][gamindex]") {

    for (id_t to_bin : {0ULL, 1ULL, 10ULL, 0xFFFFFFFFFFFFFFFFULL, 0xFACEDEADCAFEBEEFULL}) {

        auto bins = GAMIndex::bins_of_id(to_bin);

        // We need one bin per bit for all bits but the last one, plus one top-level 0 bin for everything.
        REQUIRE(bins.size() == CHAR_BIT * sizeof(id_t));

        // The bins should end with the least specific bin (0)
        REQUIRE(bins.back() == 0);
        
        // Bin levels go from 0 to bits-1.
        // The first bin should be 2^(bits - 1) - 1 + id >> 1
        REQUIRE(bins.front() == ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 1)) - 1 + ((GAMIndex::bin_t)to_bin >> 1));
        
        // The first bin is the most specific bin, which is the bin of the number and itself.
        REQUIRE(bins.front() == GAMIndex::common_bin(to_bin, to_bin));
    }
    
}

TEST_CASE("GAMindex binning works on adjacent numbers", "[gam][gamindex]") {
    
    // The common bin of an even and the next odd number should be the two numbers right shifted by one, pulss an offset.
    // For an odd and the next even, it should be them shifted by 2, pluss a smaller offset.
    
    // what offsets do we expect for bins based on their size?
    auto SIZE_2_OFFSET = ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 1)) - 1;
    auto SIZE_4_OFFSET = ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 2)) - 1;
    
    for (size_t i = -10; i < 10000; i++) {
        auto bin_found = GAMIndex::common_bin(i, i + 1);
        auto bin_found2 = GAMIndex::common_bin(i + 1, i);
        
        // Should work in any order
        REQUIRE(bin_found == bin_found2);
        
        if (i == -1) {
            // The common bin between -1 and 0 has to be bin 0, because that's where the discontinuity falls.
            // Not a problem because we don't use negative node IDs in real life.
            REQUIRE(bin_found == 0);
        } else {
            if (i % 2 == 0) {
                // Even number and next odd
                REQUIRE(bin_found == SIZE_2_OFFSET + ((GAMIndex::bin_t)i >> 1));
            } else {
                // Odd number and next even
                REQUIRE(bin_found == SIZE_4_OFFSET + ((GAMIndex::bin_t)i >> 2));
            }
        }
    }
    

}

TEST_CASE("GAMindex can look up inserted ranges", "[gam][gamindex]") {
    // Make an empty index
    GAMIndex index;

    // Add some ID-sorted groups
    index.add_group(1, 5, 0, 100);
    index.add_group(3, 7, 100, 200);
    index.add_group(6, 9, 200, 300);
    // Being sorted by lowest ID doesn't mean you are always sorted by highest ID
    index.add_group(7, 8, 300, 400);
    index.add_group(100, 110, 400, 500);
    index.add_group(1000, 1005, 500, 600);
    
    // Look for node 1
    auto found = index.find(1);
    
    // We should find the run from 0 to 100, or a set of runs encompassing that
    REQUIRE(found.size() > 0);
    REQUIRE(found.front().first <= 0);
    REQUIRE(found.back().second >= 100);
    for (size_t i = 0; i + 1 < found.size(); i++) {
        // Successive ranges shouldn't overlap, but may abut.
        REQUIRE(found[i].second <= found[i + 1].first);
    }

    // Look for node 7
    found = index.find(7);
    
    // It could occur as early as 100 or as late as before 400
    REQUIRE(found.size() > 0);
    REQUIRE(found.front().first <= 100);
    REQUIRE(found.back().second >= 400);
    for (size_t i = 0; i + 1 < found.size(); i++) {
        // Successive ranges shouldn't overlap, but may abut.
        REQUIRE(found[i].second <= found[i + 1].first);
    }

    
    // Look for node 500 which nothing can touch or be near
    found = index.find(500);
    
    REQUIRE(found.size() == 0);
    
    // Look for node 1000 which should benefit from the windowing
    found = index.find(1000);
    
    // We should find runs encompassing the run we added
    REQUIRE(found.size() > 0);
    REQUIRE(found.front().first <= 500);
    REQUIRE(found.back().second >= 600);
    
    // This should be the only thing in its window, so really we shouldn't find anything too early
    REQUIRE(found.front().first == 500);
    
    
}

TEST_CASE("GAMindex can be serialized and deserialized and still work", "[gam][gamindex]") {
    // Make an empty index
    GAMIndex build_index;

    // Add some ID-sorted groups
    build_index.add_group(1, 5, 0, 100);
    build_index.add_group(3, 7, 100, 200);
    build_index.add_group(6, 9, 200, 300);
    // Being sorted by lowest ID doesn't mean you are always sorted by highest ID
    build_index.add_group(7, 8, 300, 400);
    build_index.add_group(100, 110, 400, 500);
    build_index.add_group(1000, 1005, 500, 600);
    
    stringstream buffer;
    build_index.save(buffer);
    
    // Make another index and load it from the buffer
    GAMIndex index;
    index.load(buffer);
    
    // Look for node 1
    auto found = index.find(1);
    
    // We should find the run from 0 to 100, or a set of runs encompassing that
    REQUIRE(found.size() > 0);
    REQUIRE(found.front().first <= 0);
    REQUIRE(found.back().second >= 100);
    for (size_t i = 0; i + 1 < found.size(); i++) {
        // Successive ranges shouldn't overlap, but may abut.
        REQUIRE(found[i].second <= found[i + 1].first);
    }

    // Look for node 7
    found = index.find(7);
    
    // It could occur as early as 100 or as late as before 400
    REQUIRE(found.size() > 0);
    REQUIRE(found.front().first <= 100);
    REQUIRE(found.back().second >= 400);
    for (size_t i = 0; i + 1 < found.size(); i++) {
        // Successive ranges shouldn't overlap, but may abut.
        REQUIRE(found[i].second <= found[i + 1].first);
    }

    
    // Look for node 500 which nothing can touch or be near
    found = index.find(500);
    
    REQUIRE(found.size() == 0);
    
    // Look for node 1000 which should benefit from the windowing
    found = index.find(1000);
    
    // We should find runs encompassing the run we added
    REQUIRE(found.size() > 0);
    REQUIRE(found.front().first <= 500);
    REQUIRE(found.back().second >= 600);
    
    // This should be the only thing in its window, so really we shouldn't find anything too early
    REQUIRE(found.front().first == 500);
    
    
}

TEST_CASE("GAMIndex can work with ProtobufIterator cursors", "[gam][gamindex]") {
    // First we will fill this file with groups of alignments
    stringstream file;
    
    id_t next_id = 1;
    
    // Define a function to stamp out a group of Alignments
    auto make_group = [&](size_t count) {
    
        vector<Alignment> group;
    
        for (size_t i = 0; i < count; i++) {
            // Make a one-node alignment to each node, in order.
            group.emplace_back();
            Alignment& aln = group.back();
            auto* mapping = aln.mutable_path()->add_mapping();
            mapping->mutable_position()->set_node_id(next_id);
            next_id++;
            
            // Give the alignment some data to make it big ish
            aln.set_sequence(random_sequence(100));
        }
        
        stream::write_buffered(file, group, 0);
    };
    
    for (size_t group_number = 0; group_number < 100; group_number++) {
#ifdef debug
        cerr << "Make group " << group_number << endl;
#endif
        make_group(100);
    }
    
#ifdef debug
    cerr << "Data is " << file.str().size() << " bytes" << endl;
#endif

    // Make a cursor to read the file
    GAMIndex::cursor_t cursor(file);
    
    // Index the file
    GAMIndex index;
    index.index(cursor);
    
    // The index should be pretty small, even though we have a lot of groups
    stringstream index_data;
    index.save(index_data);
#ifdef debug
    cerr << "Index data is " << index_data.str().size() << " bytes" << endl;
#endif
    REQUIRE(index_data.str().size() < 10000);
    
    // Remember every range we look up
    vector<pair<id_t, id_t>> ranges;
    // And the number of alignments we find
    size_t total_found = 0;
    
    for (size_t start = 1; start < next_id; start += 345) {
        // Look up a series of ranges
        auto last = start + 9;
        
        // Remember each range we look up
        ranges.emplace_back(start, last);
        
        vector<id_t> seen;
        
        // Collect the visited nodes of all the alignments
        index.find(cursor, start, last, [&](const Alignment& found) {
            seen.push_back(found.path().mapping(0).position().node_id());
            total_found++;
        });
        
#ifdef debug
        cerr << "Found " << seen.size() << " alignments for range " << start << "-" << last << endl;
#endif
        
        // Make sure we found just the matching reads.
        if (last >= next_id) {
            REQUIRE(seen.size() == next_id - start);
            REQUIRE(seen.front() == start);
            REQUIRE(seen.back() == next_id - 1);
        } else {
            REQUIRE(seen.size() == 10);
            REQUIRE(seen.front() == start);
            REQUIRE(seen.back() == last);
        }
        
    }
    
    // Make sure we find the same alignment count when querying the ranges together, because they don't overlap.
    size_t recovered = 0;
    index.find(cursor, ranges, [&](const Alignment& found) {
        recovered++;
    });
    
    REQUIRE(recovered == total_found);
    
}


}
}
