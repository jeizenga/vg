/// \file aligner.cpp
///  
/// Unit tests for the basic methods of the Aligner class. See also:
/// pinned_alignment.cpp.
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include "../vg.pb.h"
#include "../gssw_aligner.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Aligner respects the full length bonus at both ends", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Aligner respects the full length bonus for a single base read", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("G");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("bonus is collected twice even though both ends are one match") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
}

TEST_CASE("Aligner works when end bonus is granted to a match at the start of a node", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    // Make sure aligner runs
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("bonus is collected twice") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Full-length bonus can hold down the left end", "[aligner][alignment][mapping]") {
    VG graph;
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTGCTGAAGT");
    
    string read = string("AATGCTGAAGT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("left end is detatched without bonus") {
        REQUIRE(aln1.path().mapping_size() == 1);
        REQUIRE(aln1.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln1.path().mapping(0).position().offset() == 2);
        REQUIRE(aln1.path().mapping(0).edit_size() == 2);
        REQUIRE(aln1.path().mapping(0).edit(0).from_length() == 0);
        REQUIRE(aln1.path().mapping(0).edit(0).sequence() == "AA");
    }
    
    SECTION("left end is attached with bonus") {
        REQUIRE(aln2.path().mapping_size() == 1);
        REQUIRE(aln2.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln2.path().mapping(0).position().offset() == 0);
        REQUIRE(aln2.path().mapping(0).edit_size() == 3);
        REQUIRE(aln2.path().mapping(0).edit(0).from_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(0).to_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(0).sequence() == "");
    }
}

TEST_CASE("Full-length bonus can hold down the right end", "[aligner][alignment][mapping]") {
    VG graph;
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTGCTGAAGT");
    
    string read = string("AGTGCTGAAAT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("right end is detatched without bonus") {
        REQUIRE(aln1.path().mapping_size() == 1);
        REQUIRE(aln1.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln1.path().mapping(0).position().offset() == 0);
        REQUIRE(aln1.path().mapping(0).edit_size() == 2);
        REQUIRE(aln1.path().mapping(0).edit(1).from_length() == 0);
        REQUIRE(aln1.path().mapping(0).edit(1).sequence() == "AT");
    }
    
    SECTION("right end is attached with bonus") {
        REQUIRE(aln2.path().mapping_size() == 1);
        REQUIRE(aln2.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln2.path().mapping(0).position().offset() == 0);
        REQUIRE(aln2.path().mapping(0).edit_size() == 3);
        REQUIRE(aln2.path().mapping(0).edit(2).from_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(2).to_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(2).sequence() == "");
    }
}

TEST_CASE("Full-length bonus can attach Ns", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("NNNNCTGANNN");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("bonused alignment ends in full-length match/mismatches") {
        REQUIRE(aln2.path().mapping_size() == 3);
        REQUIRE(mapping_from_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_to_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_from_length(aln2.path().mapping(2)) == 6);
        REQUIRE(mapping_to_length(aln2.path().mapping(2)) == 6);
    }
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Full-length bonus can attach to Ns", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("NNNG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGANNN");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("bonused alignment ends in full-length match/mismatches") {
        REQUIRE(aln2.path().mapping_size() == 3);
        REQUIRE(mapping_from_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_to_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_from_length(aln2.path().mapping(2)) == 6);
        REQUIRE(mapping_to_length(aln2.path().mapping(2)) == 6);
    }
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Full-length bonus can attach Ns to Ns", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("NNNG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGANNN");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("NNNGCTGANNN");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph, true, false);
    aligner_2.align(aln2, graph.graph, true, false);
    
    SECTION("bonused alignment ends in full-length match/mismatches") {
        REQUIRE(aln2.path().mapping_size() == 3);
        REQUIRE(mapping_from_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_to_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_from_length(aln2.path().mapping(2)) == 6);
        REQUIRE(mapping_to_length(aln2.path().mapping(2)) == 6);
    }
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Full-length bonus is applied to both ends by rescoring", "[aligner][alignment][scoring]") {
    
    string aln_str = R"({"sequence":"ACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGTGTGGTGGCATGCACCTGTAATCCCAGCTACTGGGCATGCTGAGGTAGCAGAATCGCTTGAACCCAGGAGGAACCGGTTGCAGTGAGCCGAGATTGTGCCACTCCACTCCAG","path":{"mapping":[{"position":{"node_id":2048512,"offset":21},"edit":[{"from_length":4,"to_length":4}],"rank":1},{"position":{"node_id":2048514},"edit":[{"from_length":1,"to_length":1}],"rank":2},{"position":{"node_id":2048515},"edit":[{"from_length":3,"to_length":3}],"rank":3},{"position":{"node_id":2048517},"edit":[{"from_length":1,"to_length":1}],"rank":4},{"position":{"node_id":2048518},"edit":[{"from_length":32,"to_length":32}],"rank":5},{"position":{"node_id":2048519},"edit":[{"from_length":32,"to_length":32}],"rank":6},{"position":{"node_id":2048520},"edit":[{"from_length":8,"to_length":8}],"rank":7},{"position":{"node_id":2048521},"edit":[{"from_length":1,"to_length":1}],"rank":8},{"position":{"node_id":2048523},"edit":[{"from_length":24,"to_length":24}],"rank":9},{"position":{"node_id":2048524},"edit":[{"from_length":1}],"rank":10},{"position":{"node_id":2048526},"edit":[{"from_length":2},{"from_length":3,"to_length":3},{"to_length":3,"sequence":"CCG"},{"from_length":27,"to_length":27}],"rank":11},{"position":{"node_id":2048527},"edit":[{"from_length":9,"to_length":9}],"rank":12}]},"fragment":[{"name":"21","length":413}]})";

    Alignment aln;
    json2pb(aln, aln_str.c_str(), aln_str.size());
    
    // Make an aligner with a full lenth bonus of 5
    Aligner aligner1(1, 4, 6, 1, 5);
    // And one with no bonus
    Aligner aligner2(1, 4, 6, 1, 0);

    REQUIRE(!softclip_start(aln));
    REQUIRE(!softclip_end(aln));

    // Normal score would be 129
    REQUIRE(aligner2.score_ungapped_alignment(aln) == 129);
    // And with a full length bonus at each end it's 139.
    REQUIRE(aligner1.score_ungapped_alignment(aln) == 139);
}

TEST_CASE("BaseAligner mapping quality estimation is robust", "[aligner][alignment][mapping][mapq]") {
    
    vector<double> scaled_scores;
    size_t max_idx;
    
    SECTION("exact mapping quality is robust") {
        
        // Empty vector disallowed
        
        SECTION("a 1-element positive vector has its element chosen") {
            scaled_scores = {10};
            BaseAligner::maximum_mapping_quality_exact(scaled_scores, &max_idx);
            REQUIRE(max_idx == 0);
        }
        
        SECTION("a 1-element zero vector has its element chosen") {
            scaled_scores = {0};
            BaseAligner::maximum_mapping_quality_exact(scaled_scores, &max_idx);
            REQUIRE(max_idx == 0);
        }
        
        SECTION("a 1-element negative vector has its element chosen") {
            scaled_scores = {-10};
            BaseAligner::maximum_mapping_quality_exact(scaled_scores, &max_idx);
            REQUIRE(max_idx == 0);
        }
        
        SECTION("a multi-element vector has a maximal element chosen") {
            scaled_scores = {1, 5, 2, 5, 4};
            BaseAligner::maximum_mapping_quality_exact(scaled_scores, &max_idx);
            REQUIRE(max_idx >= 1);
            REQUIRE(max_idx != 2);
            REQUIRE(max_idx <= 3);
        }
    }
    
    SECTION("inexact mapping quality is robust") {

        // Empty vector disallowed

        SECTION("a 1-element positive vector has its element chosen") {
            scaled_scores = {10};
            BaseAligner::maximum_mapping_quality_approx(scaled_scores, &max_idx);
            REQUIRE(max_idx == 0);
        }
        
        SECTION("a 1-element zero vector has its element chosen") {
            scaled_scores = {0};
            BaseAligner::maximum_mapping_quality_approx(scaled_scores, &max_idx);
            REQUIRE(max_idx == 0);
        }
        
        SECTION("a 1-element negative vector has its element chosen") {
            scaled_scores = {-10};
            BaseAligner::maximum_mapping_quality_approx(scaled_scores, &max_idx);
            REQUIRE(max_idx == 0);
        }
        
        SECTION("a multi-element vector has a maximal element chosen") {
            scaled_scores = {1, 5, 2, 5, 4};
            BaseAligner::maximum_mapping_quality_approx(scaled_scores, &max_idx);
            REQUIRE(max_idx >= 1);
            REQUIRE(max_idx != 2);
            REQUIRE(max_idx <= 3);
        }
    
    }
    
}
   
}
}
        
