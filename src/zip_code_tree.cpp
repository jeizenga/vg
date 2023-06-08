#define DEBUG_ZIP_CODE_TREE

#include "zip_code_tree.hpp"

#include "crash.hpp"

using namespace std;
namespace vg {

ZipCodeTree::ZipCodeTree(vector<Seed>& seeds, const SnarlDistanceIndex& distance_index) :
    seeds(seeds) {

    /*
    Constructor for the ZipCodeTree
    Takes a vector of seeds and constructs the tree

    Tree construction is done by first sorting the seeds along chains/snarls
    Then, adding each seed, snarl/chain boundary, and distance to zip_code_tree
    Finally (optionally), the tree is refined to take out unnecessary edges
    */

    //////////////////// Sort the seeds

    //A vector of indexes into seeds
    //To be sorted along each chain/snarl the snarl tree
    vector<size_t> seed_indices (seeds.size(), 0);
    for (size_t i = 0 ; i < seed_indices.size() ; i++) {
        seed_indices[i] = i;
    }

    //Sort the indices
    std::sort(seed_indices.begin(), seed_indices.end(), [&] (const size_t& a, const size_t& b) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Comparing seeds " << seeds[a].pos << " and " << seeds[b].pos << endl;
#endif
        //Comparator returning a < b
        size_t depth = 0;

        //Keep track of the orientation of each seed
        //Everything should be sorted according to the orientation in the top-level structure,
        //so if things are traversed backwards, reverse the orientation
        bool a_is_reversed = false;
        bool b_is_reversed = false;
        while (depth < seeds[a].zipcode_decoder->max_depth() &&
               depth < seeds[b].zipcode_decoder->max_depth() &&
               ZipCodeDecoder::is_equal(*seeds[a].zipcode_decoder, *seeds[b].zipcode_decoder, depth)) {
            if (seeds[a].zipcode_decoder->get_is_reversed_in_parent(depth)) {
                a_is_reversed = !a_is_reversed;
            }
            if (seeds[b].zipcode_decoder->get_is_reversed_in_parent(depth)) {
                b_is_reversed = !b_is_reversed;
            }
            depth++;
        }

        //Check the orientations one last time
        if (seeds[a].zipcode_decoder->get_is_reversed_in_parent(depth)) {
            a_is_reversed = !a_is_reversed;
        }
        if (seeds[b].zipcode_decoder->get_is_reversed_in_parent(depth)) {
            b_is_reversed = !b_is_reversed;
        }
#ifdef DEBUG_ZIP_CODE_TREE
        //cerr << "\t different at depth " << depth << endl;
#endif
        //Either depth is the last thing in a or b, or they are different at this depth


        if ( ZipCodeDecoder::is_equal(*seeds[a].zipcode_decoder, *seeds[b].zipcode_decoder, depth)) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tthey are on the same node" << endl;
#endif
            //If they are equal, then they must be on the same node

            size_t offset1 = is_rev(seeds[a].pos)
                           ? seeds[a].zipcode_decoder->get_length(depth) - offset(seeds[a].pos) - 1
                           : offset(seeds[a].pos);
            size_t offset2 = is_rev(seeds[b].pos)
                           ? seeds[b].zipcode_decoder->get_length(depth) - offset(seeds[b].pos) - 1
                           : offset(seeds[b].pos);
            if (!a_is_reversed) {
                //If they are in a snarl or they are facing forward on a chain, then order by
                //the offset in the node
                return offset1 < offset2;
            } else {
                //Otherwise, the node is facing backwards in the chain, so order backwards in node
                return offset2 < offset1;
            }
        }  else if (depth == 0) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            //cerr << "\tThey are on different connected components" << endl;
#endif
            //If they are on different connected components, sort by connected component
            return seeds[a].zipcode_decoder->get_distance_index_address(0) < seeds[b].zipcode_decoder->get_distance_index_address(0);
            
        }  else if (seeds[a].zipcode_decoder->get_code_type(depth-1) == CHAIN || seeds[a].zipcode_decoder->get_code_type(depth-1) == ROOT_CHAIN) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            //cerr << "\t they are children of a common chain" << endl;
#endif
            //If a and b are both children of a chain
            size_t offset_a = seeds[a].zipcode_decoder->get_offset_in_chain(depth);
            size_t offset_b = seeds[b].zipcode_decoder->get_offset_in_chain(depth);
            if ( offset_a == offset_b) {
                //If they have the same prefix sum, then the snarl comes first
                //They will never be on the same child at this depth
                return seeds[a].zipcode_decoder->get_code_type(depth) != NODE && seeds[b].zipcode_decoder->get_code_type(depth) == NODE;  
            } else {
                return offset_a < offset_b;
            }
        } else if (seeds[a].zipcode_decoder->get_code_type(depth-1) == REGULAR_SNARL) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            //cerr << "\t they are children of a common regular snarl" << endl;
#endif
            //If the parent is a regular snarl, then sort by order along the parent chain
            size_t offset1 = is_rev(seeds[a].pos) 
                           ? seeds[a].zipcode_decoder->get_length(depth) - offset(seeds[a].pos) - 1
                           : offset(seeds[a].pos); 
            size_t offset2 = is_rev(seeds[b].pos) 
                           ? seeds[b].zipcode_decoder->get_length(depth) - offset(seeds[b].pos) - 1
                           : offset(seeds[b].pos);
            if (a_is_reversed) {
                return offset1 < offset2;
            } else {
                return offset2 < offset1;
            }
        } else {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            //cerr << "\t they are children of a common irregular snarl" << endl;
#endif
            //Otherwise, they are children of an irregular snarl
            //Sort by the distance to the start of the irregular snarl
            size_t distance_to_start_a = seeds[a].zipcode_decoder->get_distance_to_snarl_start(depth);
            size_t distance_to_start_b = seeds[b].zipcode_decoder->get_distance_to_snarl_start(depth);
            if (distance_to_start_a == distance_to_start_b) {
                //If they are equi-distant to the start of the snarl, then put the one that is
                //farther from the end first

                return seeds[a].zipcode_decoder->get_distance_to_snarl_end(depth) >
                         seeds[b].zipcode_decoder->get_distance_to_snarl_end(depth);
            } else {
                return distance_to_start_a < distance_to_start_b;
            }
        } 
    });

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Sorted positions:" << endl;
    for (const size_t& i : seed_indices) {
        cerr << seeds[i].pos << endl;
    }
#endif

    //seed_indices is now sorted roughly along snarls and chains


    ///////////////////// Build the tree

    //For children of snarls, we need to remember the siblings and start bound that came before them
    //so we can record their distances
    //This holds the indices (into zip_code_tree) of each seed or start of a chain,
    // and each start and child chain start of a snarl
    //The children are stored at the depth of their parents. For example, for a root chain,
    //the vector at index 0 would have the chain start, seeds that are on the chain, and the start
    //of snarls on the chain. Similarly, for a top-level snarl at depth 1, the second vector would contain
    //the starts of chains at depth 2 
    //For the children of a chain, the value is the prefix sum in the chain (relative to the orientation 
    //of the top-level chain, not necessarily the chain itself)
    //For the children of a snarl, the value is the index of the seed
    struct child_info_t {
        tree_item_type_t type;  //the type of the item
        size_t value;  //A value associated with the item, could be offset in a chain, index of the seed
    };
    vector<vector<child_info_t>> sibling_indices_at_depth;

    /* The tree will hold all seeds and the bounds of snarls and chains
       For each chain, there must be a distance between each element of the chain (seeds and snarls)
       For each snarl, each element (chain or boundary) is preceded by the distances to everything
         before it in the snarl.
    */

    for (size_t i = 0 ; i < seed_indices.size() ; i++) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "At " << i << "st/nd/th seed: " << seeds[seed_indices[i]].pos << endl;
#endif

        //1. First, find the lowest common ancestor with the previous seed.
        //2. To finish the ancestors of the previous seed that are different from this one,
        //   walk up the snarl tree from the previous max depth and mark the end of the ancestor,
        //   adding distances for snarl ends 
        //3. To start anything for this seed, start from the first ancestor that is different
        //   and walk down the snarl tree, adding distances for each ancestor

        Seed& current_seed = seeds[seed_indices[i]];

        size_t current_max_depth = current_seed.zipcode_decoder->max_depth();
        //Make sure sibling_indices_at_depth has enough spaces for this zipcode
        while (sibling_indices_at_depth.size() < current_max_depth+1) {
            sibling_indices_at_depth.emplace_back();
        }

        //Get the previous seed (if this isn't the first one)
        Seed& previous_seed = i == 0 ? current_seed : seeds[seed_indices[i-1]];
        //And the previous max depth
        size_t previous_max_depth = i == 0 ? 0 : previous_seed.zipcode_decoder->max_depth();

        //Remember the orientation for the seeds at the current depth
        //We start the first traversal (2) from previous_max_depth
        //The second traversal (3) starts from first_different_ancestor_depth 
        //This one is for the first traversal, so it will be for previous_max_depth
        bool previous_is_reversed = false;
        //This is for the second traversal, find it when finding first_different_ancestor_depth
        bool current_is_reversed = false;


        //Find the depth at which the two seeds are on different snarl tree nodes
        size_t first_different_ancestor_depth = 0;
        bool same_node = false;
        size_t max_depth = std::min(current_max_depth, previous_max_depth);

        for (size_t depth = 0 ; depth <= max_depth ; depth++) {
            first_different_ancestor_depth = depth;
            current_is_reversed = current_seed.zipcode_decoder->get_is_reversed_in_parent(depth)
                                    ? !current_is_reversed : current_is_reversed;
            if (i != 0) {
                previous_is_reversed = previous_seed.zipcode_decoder->get_is_reversed_in_parent(depth)
                                        ? !previous_is_reversed : previous_is_reversed;
            }
            cerr << "At depth " << depth << " is reversed? " << current_is_reversed << endl;
            if (!ZipCodeDecoder::is_equal(*current_seed.zipcode_decoder, 
                        *previous_seed.zipcode_decoder, depth)) {
                break;
            } else if (depth == max_depth) {
                same_node = true;
            }
        }
        if (previous_max_depth > current_max_depth) {
            //We might need to update previous_is_reversed
            for (size_t depth = max_depth ; depth <= previous_max_depth ; depth++) {
                previous_is_reversed = previous_seed.zipcode_decoder->get_is_reversed_in_parent(depth)
                                    ? !previous_is_reversed : previous_is_reversed;
            }
        }
        if (i == 0) { 
            same_node = false;
        }
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\tthe depth of the first ancestor different than the previous seed is " << first_different_ancestor_depth << endl;
        cerr << "\tWalk up the snarl tree from depth " << previous_max_depth << " and close any snarl/chains" << endl;
#endif

        //Now, close anything that ended at the previous seed, starting from the leaf of the previous seed
        //If there was no previous seed, then the loop is never entered
        for (int depth = previous_max_depth ; !same_node && i!=0 && depth >= first_different_ancestor_depth && depth >= 0 ; depth--) {
            code_type_t previous_type = previous_seed.zipcode_decoder->get_code_type(depth);
            cerr << "At depth " << depth << " previous type was " << previous_type << endl;
            if (previous_type == CHAIN || previous_type == ROOT_CHAIN || previous_type == ROOT_NODE) {
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t\tclose a chain at depth " << depth << endl;
#endif
                //If this is the end of a chain, then add the distance from the last child to the end

                //If this is reversed, then the distance should be the distance to the start of 
                //the chain. Otherwise, the distance to the end
                //The value that got stored in sibling_indices_at_depth was the prefix sum
                //traversing the chain according to its orientation in the tree, so either way
                //the distance is the length of the chain - the prefix sum
                // TODO: When we get C++20, change this to emplace_back aggregate initialization
                if (previous_type == CHAIN) {
                    //Only add the distance for a non-root chain
                    zip_code_tree.push_back({EDGE, 
                        SnarlDistanceIndex::minus(previous_seed.zipcode_decoder->get_length(depth),
                                                  sibling_indices_at_depth[depth].back().value)});
                }

                zip_code_tree.push_back({CHAIN_END, std::numeric_limits<size_t>::max()});


            } else if (previous_type == REGULAR_SNARL || previous_type == IRREGULAR_SNARL) { 
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t\tclose a snarl at depth " << depth << endl;
#endif
                //If this is the end of the snarl, then we need to save the distances to 
                //all previous children of the snarl

                zip_code_tree.resize(zip_code_tree.size() + sibling_indices_at_depth[depth].size());

                for (size_t sibling_i = 0 ; sibling_i < sibling_indices_at_depth[depth].size() ; sibling_i++) {
                    const auto& sibling = sibling_indices_at_depth[depth][sibling_i];
                    if (sibling.type == SNARL_START) {
                        //First, the distance between ends of the snarl, which is the length
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            previous_seed.zipcode_decoder->get_length(depth)};
                            cerr << "Add distance to snarl start" << endl;
                    } else {
                        //For the rest of the children, find the distance from the child to
                        //the end
                        //If the child is reversed relative to the top-level chain, then get the distance to start
                                cerr << "Add distance to child with index " << sibling.value<< endl;
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            previous_is_reversed 
                                ? seeds[sibling.value].zipcode_decoder->get_distance_to_snarl_start(depth)
                                :seeds[sibling.value].zipcode_decoder->get_distance_to_snarl_end(depth)};

                    }
                }
                //Note the count of children and the end of the snarl
                zip_code_tree.push_back({NODE_COUNT, sibling_indices_at_depth[depth].size()-1});
                zip_code_tree.push_back({SNARL_END, std::numeric_limits<size_t>::max()});
            }
            //Update previous_is_reversed to the one before this
            previous_is_reversed = (depth > 0 && previous_seed.zipcode_decoder->get_is_reversed_in_parent(depth-1))
                                        ? !previous_is_reversed : previous_is_reversed;

            //Clear the list of children of the thing at this level
            sibling_indices_at_depth[depth].clear();
        }
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\tWalk down the snarl tree from depth " << first_different_ancestor_depth << " to " << current_max_depth  << " and open any snarl/chains" << endl;
#endif

        //Now go through everything that started a new snarl tree node going down the snarl tree
        //For each new snarl or seed in a chain, add the distance to the thing preceding it in the chain
        //For each new chain in a snarl, add the distance to everything preceding it in the snarl
        //If this is the same node as the previous, then first_different_ancestor_depth is the depth 
        //of the node
        for (size_t depth = first_different_ancestor_depth ; depth <= current_max_depth ; depth++) {
            code_type_t current_type = current_seed.zipcode_decoder->get_code_type(depth);
            cerr << "At depth " << depth << endl;

            if (current_type == NODE || current_type == REGULAR_SNARL || current_type == IRREGULAR_SNARL
                || current_type == ROOT_NODE) {
                //For these things, we need to remember the offset in the node/chain

                if (current_type == ROOT_NODE && sibling_indices_at_depth[depth].empty()) {
                    //If this is a root-level node and the first time we've seen it,
                    //then open the node
                    zip_code_tree.push_back({CHAIN_START, std::numeric_limits<size_t>::max()});
                    sibling_indices_at_depth[depth].push_back({CHAIN_START, 0});
                }

                ///////////////// Get the offset in the parent chain (or node)
                size_t current_offset;

                //If we're traversing this chain backwards, then the offset is the offset from the end
                bool current_parent_is_reversed = current_seed.zipcode_decoder->get_is_reversed_in_parent(depth) 
                    ? !current_is_reversed : current_is_reversed;

                //First, get the prefix sum in the chain
                if (current_type == ROOT_NODE) {
                    //Which is 0 if this is just a node
                    current_offset = 0;
                } else {
                    //And the distance to the start or end of the chain if it's a node/snarl in a chain
                    current_offset = current_parent_is_reversed 
                            ? SnarlDistanceIndex::minus(current_seed.zipcode_decoder->get_length(depth-1) ,
                                                        SnarlDistanceIndex::sum(
                                                            current_seed.zipcode_decoder->get_offset_in_chain(depth),
                                                            current_seed.zipcode_decoder->get_length(depth))) 
                            : current_seed.zipcode_decoder->get_offset_in_chain(depth);
                }

                if (depth == current_max_depth) {
                    //If this is a node, then add the offset of the position in the node
                    current_offset = SnarlDistanceIndex::sum(current_offset, 
                        current_is_reversed != is_rev(current_seed.pos)
                            ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                            : offset(current_seed.pos)+1);
                }

                /////////////////////// Get the offset of the previous thing in the parent chain/node
                size_t previous_offset = depth == 0 ? sibling_indices_at_depth[depth][0].value 
                                                    : sibling_indices_at_depth[depth-1][0].value;

#ifdef DEBUG_ZIP_CODE_TREE
                if (depth > 0) {
                    assert(sibling_indices_at_depth[depth-1].size() == 1);
                }
                cerr << current_offset << " " << previous_offset << endl;
                assert(current_offset >= previous_offset);
#endif

                ///////////////////// Record the distance from the previous thing in the chain/node
                if (!(depth == 0 && sibling_indices_at_depth[depth][0].type == CHAIN_START) &&
                    !(depth == 1 && current_seed.zipcode_decoder->get_code_type(depth-1) == ROOT_CHAIN &&
                       sibling_indices_at_depth[depth-1][0].type == CHAIN_START)) {
                    //for everything except the first thing in a root node, or root chain
                    zip_code_tree.push_back({EDGE, current_offset-previous_offset});
                }

                /////////////////////////////Record this thing in the chain
                if (current_type == NODE || current_type == ROOT_NODE) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t\tContinue node/chain with seed " << seeds[seed_indices[i]].pos << " at depth " << depth << endl;
#endif
                    //If this was a node, just remember the seed
                    zip_code_tree.push_back({SEED, seed_indices[i]});
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t\tOpen new snarl at depth " << depth << endl;
#endif
                    //If this was a snarl, record the start of the snarl
                    zip_code_tree.push_back({SNARL_START, std::numeric_limits<size_t>::max()});

                    //Remember the start of the snarl
                    sibling_indices_at_depth[depth].push_back({SNARL_START, std::numeric_limits<size_t>::max()});

                    //For finding the distance to the next thing in the chain, the offset
                    //stored should be the offset of the end bound of the snarl, so add the 
                    //length of the snarl
                    current_offset = SnarlDistanceIndex::sum(current_offset,
                        current_seed.zipcode_decoder->get_length(depth));

                }

                //Remember this thing for the next sibling in the chain
                if (depth == 0) {
                    sibling_indices_at_depth[depth].pop_back();
                    sibling_indices_at_depth[depth].push_back({SEED, current_offset}); 
                } else {
                    sibling_indices_at_depth[depth-1].pop_back();
                    //THis may or may not be a seed but it doesn't matter, as long as its a child of a chain
                    sibling_indices_at_depth[depth-1].push_back({SEED, current_offset}); 
                }
            } else {
                //Otherwise, this is a chain or root chain
                //If it is a chain, then it is the child of a snarl, so we need to find distances
                //to everything preceding it in the snarl
                assert(current_type == CHAIN || current_type == ROOT_CHAIN);
                if (sibling_indices_at_depth[depth].size() == 0) {
                    //If this is the start of a new chain
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t\tOpen new chain at depth " << depth << endl;
#endif

                    //For each sibling in the snarl, record the distance from the sibling to this
                    if (current_type == CHAIN) {
                        //If this is the start of a non-root chain, then it is the child of a snarl and 
                        //we need to find the distances to the previous things in the snarl

                        //The distances will be added in reverse order that they were found in
                        zip_code_tree.resize(zip_code_tree.size() + sibling_indices_at_depth[depth-1].size());
                        for ( size_t sibling_i = 0 ; sibling_i < sibling_indices_at_depth[depth-1].size() ; sibling_i++) {
                            const auto& sibling = sibling_indices_at_depth[depth-1][sibling_i];
                            if (sibling.type == SNARL_START) {
                                cerr << "Add distance to sibling start" << endl;
                                zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = 
                                 {EDGE, 
                                    current_is_reversed
                                        ? current_seed.zipcode_decoder->get_distance_to_snarl_end(depth)
                                        : current_seed.zipcode_decoder->get_distance_to_snarl_start(depth)};
                            } else {
                                //Otherwise, the previous thing was another child of the snarl
                                //and we need to record the distance between these two
                                //TODO: This can be improved for simple snarls
                                size_t distance;
                                if (current_type == CHAIN && 
                                    current_seed.zipcode_decoder->get_code_type(depth-1) == REGULAR_SNARL) {
                                    //If this is the child of a regular snarl, then the distance between
                                    //any two chains is inf
                                    distance = std::numeric_limits<size_t>::max();
                                } else {
                                    net_handle_t snarl_handle = current_seed.zipcode_decoder->get_net_handle(depth-1, &distance_index);
                                    size_t rank2 = current_seed.zipcode_decoder->get_rank_in_snarl(depth);
                                    size_t rank1 = seeds[sibling.value].zipcode_decoder->get_rank_in_snarl(depth);
                                    //TODO: idk about this distance- I think the orientations need to change
                                    distance = distance_index.distance_in_snarl(snarl_handle, rank1, false, rank2, false);
                                }
                                zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE, distance};
                            }

                        }
                    }

                    //Now record the start of this chain
                    zip_code_tree.push_back({CHAIN_START, std::numeric_limits<size_t>::max()});

                    //Remember the start of the chain, with the prefix sum value
                    sibling_indices_at_depth[depth].push_back({CHAIN_START, 0});

                    //And, if it is the child of a snarl, then remember the chain as a child of the snarl
                    if (depth != 0) {
                        sibling_indices_at_depth[depth-1].push_back({CHAIN_START,
                                                                     seed_indices[i]});
                    }
                }

                if (current_type == CHAIN && depth == current_max_depth) {
                    //If this is a trivial chain, then also add the seed and the distance to the 
                    //thing before it
                    size_t current_offset = current_is_reversed
                            ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                            : offset(current_seed.pos)+1;

                    zip_code_tree.push_back({EDGE, current_offset - sibling_indices_at_depth[depth].back().value}); 
                    zip_code_tree.push_back({SEED, seed_indices[i]}); 

                    //And update sibling_indices_at_depth to remember this child
                    sibling_indices_at_depth[depth].pop_back();
                    sibling_indices_at_depth[depth].push_back({SEED, current_offset});
                    
                }
            }
            
            //Finished with this depth, so update current_is_reversed to be for the next ancestor
            current_is_reversed = depth < current_max_depth && current_seed.zipcode_decoder->get_is_reversed_in_parent(depth+1)
                                    ? !current_is_reversed : current_is_reversed;
        }


    }
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Close any snarls or chains that remained open" << endl;
#endif

    // Now close anything that remained open
    const Seed& last_seed = seeds[seed_indices.back()];
    size_t last_max_depth = last_seed.zipcode_decoder->max_depth();
    print_self();

    //Find out if this seed is reversed at the leaf of the snarl tree (the node)
    bool last_is_reversed = false;
    for (size_t depth = 0 ; depth <= last_max_depth ; depth++) {
        if (last_seed.zipcode_decoder->get_is_reversed_in_parent(depth)) {
            last_is_reversed = !last_is_reversed;
        }
    }
    for (int depth = last_max_depth ; depth >= 0 ; depth--) {
        cerr << "At depth " << depth << endl;
        print_self();
        if (sibling_indices_at_depth[depth].size() > 0) {
            code_type_t last_type = last_seed.zipcode_decoder->get_code_type(depth);
            if (last_type == CHAIN || last_type == ROOT_CHAIN || last_type == ROOT_NODE) {
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t\tclose a chain at depth " << depth << endl;
#endif
                //If this is the end of a chain, then add the distance from the last child to the end

                //If this is reversed, then the distance should be the distance to the start of 
                //the chain. Otherwise, the distance to the end
                //The value that got stored in sibling_indices_at_depth was the prefix sum
                //traversing the chain according to its orientation in the tree, so either way
                //the distance is the length of the chain - the prefix sum
                // TODO: When we get C++20, change this to emplace_back aggregate initialization
                if (last_type == CHAIN) {
                    zip_code_tree.push_back({EDGE, 
                        SnarlDistanceIndex::minus(last_seed.zipcode_decoder->get_length(depth),
                                                  sibling_indices_at_depth[depth].back().value)});
                }

                zip_code_tree.push_back({CHAIN_END, std::numeric_limits<size_t>::max()});

            } else if (last_type == REGULAR_SNARL || last_type == IRREGULAR_SNARL) { 
#ifdef DEBUG_ZIP_CODE_TREE
               cerr << "\t\tclose a snarl at depth " << depth << endl;
#endif
                //If this is the end of the snarl, then we need to save the distances to 
                //all previous children of the snarl

                zip_code_tree.resize(zip_code_tree.size() + sibling_indices_at_depth[depth].size());

                for (size_t sibling_i = 0 ; sibling_i < sibling_indices_at_depth[depth].size() ; sibling_i++) {
                    const auto& sibling = sibling_indices_at_depth[depth][sibling_i];
                    if (sibling.type == SNARL_START) {
                        cerr << "Add length " << endl;
                        //First, the distance between ends of the snarl, which is the length
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            last_seed.zipcode_decoder->get_length(depth)};
                    } else {
                        //For the rest of the children, find the distance from the child to
                        //the end
                        //If the child is reversed relative to the top-level chain, then get the distance to start
                        cerr << "Get distance to snarl start" << endl;
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            last_is_reversed 
                                ? seeds[sibling.value].zipcode_decoder->get_distance_to_snarl_start(depth)
                                : seeds[sibling.value].zipcode_decoder->get_distance_to_snarl_end(depth)};
                    }
                }
                //Note the count of children and the end of the snarl
                zip_code_tree.push_back({NODE_COUNT, sibling_indices_at_depth[depth].size()-1});
                zip_code_tree.push_back({SNARL_END, std::numeric_limits<size_t>::max()});
            }
        }
        //Update last_is_reversed to the one before this
        last_is_reversed = (depth > 0 && last_seed.zipcode_decoder->get_is_reversed_in_parent(depth-1))
                                    ? !last_is_reversed : last_is_reversed;
    }
}

void ZipCodeTree::print_self() const {
    for (const tree_item_t item : zip_code_tree) {
        if (item.type == SEED) {
            cerr << seeds[item.value].pos;
        } else if (item.type == SNARL_START) {
            cerr << "(";
        } else if (item.type == SNARL_END) {
            cerr << ")";
        } else if (item.type == CHAIN_START) {
            cerr << "[";
        } else if (item.type == CHAIN_END) {
            cerr << "]";
        } else if (item.type == EDGE) {
            cerr << " " << item.value << " ";
        } else if (item.type == NODE_COUNT) {
            cerr << " " << item.value;
        } else {
            throw std::runtime_error("[zip tree]: Trying to print a zip tree item of the wrong type");
        }
    }
    cerr << endl;
}


ZipCodeTree::iterator::iterator(vector<tree_item_t>::const_iterator it, vector<tree_item_t>::const_iterator end) : it(it), end(end) {
    // Nothing to do!
}

auto ZipCodeTree::iterator::operator++() -> iterator& {
    ++it;
    while (it != end && it->type != SEED) {
        // Advance to the next seed, or the end.
        ++it;
    }
    return *this;
}

auto ZipCodeTree::iterator::operator==(const iterator& other) const -> bool {
    // Ends don't matter for comparison.
    return it == other.it;
}
    
auto ZipCodeTree::iterator::operator*() const -> size_t {
    return it->value;
}

auto ZipCodeTree::iterator::remaining_tree() const -> size_t {
    return end - it;
}

auto ZipCodeTree::begin() const -> iterator {
    return iterator(zip_code_tree.begin(), zip_code_tree.end());
}

auto ZipCodeTree::end() const -> iterator {
    return iterator(zip_code_tree.end(), zip_code_tree.end());
}

ZipCodeTree::reverse_iterator::reverse_iterator(vector<tree_item_t>::const_reverse_iterator it, vector<tree_item_t>::const_reverse_iterator rend, size_t distance_limit) : it(it), rend(rend), distance_limit(distance_limit), stack(), current_state(S_START) {
    while (it != rend && !tick()) {
        // Skip ahead to the first seed we actually want to yield, or to the end of the data.
        ++it;
    }
    // As the end of the constructor, the iterator points to a seed that has been ticked and yielded, or is rend.
}

auto ZipCodeTree::reverse_iterator::operator++() -> reverse_iterator& {
    // Invariant: the iterator points to a seed that has been ticked and yielded, or to rend.
    if (it != rend) {
        ++it;
    }
    while (it != rend && !tick()) {
        // Skip ahead to the next seed we actually want to yield, or to the end of the data.
        ++it;
    }
    return *this;
}

auto ZipCodeTree::reverse_iterator::operator==(const reverse_iterator& other) const -> bool {
    // Ends and other state don't matter for comparison.
    return it == other.it;
}

auto ZipCodeTree::reverse_iterator::operator*() const -> std::pair<size_t, size_t> {
    // We are always at a seed, so show that seed
    crash_unless(it != rend);
    crash_unless(it->type == SEED);
    crash_unless(!stack.empty());
    // We know the running distance to this seed will be at the top of the stack.
    return {it->value, stack.top()};
}

auto ZipCodeTree::reverse_iterator::push(size_t value) -> void {
    stack.push(value);
}

auto ZipCodeTree::reverse_iterator::pop() -> size_t {
    size_t value = stack.top();
    stack.pop();
    return value;
}

auto ZipCodeTree::reverse_iterator::top() -> size_t& {
    return stack.top();
}

auto ZipCodeTree::reverse_iterator::dup() -> void {
    push(stack.top());
}

auto ZipCodeTree::reverse_iterator::depth() const -> size_t {
    return stack.size();
}

auto ZipCodeTree::reverse_iterator::reverse(size_t depth) -> void {
    // We reverse by moving from a stack to a queue and back.
    // TODO: would using a backing vector and STL algorithms be better?
    std::queue<size_t> queue;
    for (size_t i = 0; i < depth; i++) {
        queue.push(stack.top());
        stack.pop();
    }
    for (size_t i = 0; i < depth; i++) {
        stack.push(queue.front());
        queue.pop();
    }
}

auto ZipCodeTree::reverse_iterator::state(State new_state) -> void {
    current_state = new_state;
}

auto ZipCodeTree::reverse_iterator::halt() -> void {
    it = rend;
}

auto ZipCodeTree::reverse_iterator::tick() -> bool {
    switch (current_state) {
    case S_START:
        // Stack is empty and we must be at a seed to start at.
        switch (it->type) {
        case SEED:
            push(0);
            state(S_SCAN_CHAIN);
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_SCAN_CHAIN:
        // Stack has at the top the running distance along the chain, and under
        // that running distances to use at the other chains in the snarl, and
        // under that running distances to use for the other chains in the
        // snarl's parent snarl, etc.
        switch (it->type) {
        case SEED:
            // Emit seed here with distance at top of stack.
            return true;
            break;
        case SNARL_END:
            // Running distance along chain is on stack, and will need to be added to all the stored distances.
            push(0); // Depth of stack that needs reversing after we read all the distances into it
            state(S_STACK_SNARL); // Stack up pre-made scratch distances for all the distances right left of here
            break;
        case CHAIN_START:
            if (depth() == 1) {
                // We never entered the parent snarl of this chain, so stack up
                // the distances left of here as options added to the
                // distance along this chain.
                //
                // Running distance along chain is on stack, and will need to
                // be added to all the stored distances.
                push(0); // Depth of stack that needs reversing after we read all the distances into it
                state(S_STACK_SNARL);
            } else {
                // We did enter the parent snarl already.
                // Discard the running distance along this chain, which no longer matters.
                pop();
                // Running distance for next chain, or running distance to cross the snarl, will be under it.
                state(S_SCAN_SNARL);
            }
            break;
        case EDGE:
            // Distance between things in a chain.
            // Add value into running distance.
            top() += it->value;
            if (top() > distance_limit) {
                // Skip over the rest of this chain
                if (depth() == 1) {
                    // We never entered the parent snarl of this chain.
                    // So if the distance along the chain is too much, there are not going to be any results with a smaller distance.
                    halt();
                } else {
                    // We need to try the next thing in the parent snarl, so skip the rest of the chain.
                    // We're skipping in 0 nested snarls right now.
                    push(0);
                    state(S_SKIP_CHAIN);
                }
            }
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_STACK_SNARL:
        // Stack has at the top the number of edges we have stacked up, and
        // under that the running distance along the parent chain, and under
        // that the stacked running distances for items in the snarl.
        switch (it->type) {
        case EDGE:
            // Swap top 2 elements to bring parent running distance to the top
            reverse(2);
            // Duplicate it
            dup();
            // Add in the edge value to make a running distance for the thing this edge is for
            top() += it->value;
            // Flip top 3 elements, so now edge count is on top, over parent running distance, over edge running distance.
            reverse(3);
            // Add 1 to the edge count
            top()++;
            break;
        case CHAIN_END:
            // Bring parent running distance above edge count
            reverse(2);
            // Throw it out
            pop();
            // Re-order all the edge running distances so we can pop them in the order we encounter the edge targets.
            reverse(pop());
            if (top() > distance_limit) {
                // Running distance is already too high so skip over the chain
                push(0);
                state(S_SKIP_CHAIN);
            } else {
                // Do the chain
                state(S_SCAN_CHAIN);
            }
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_SCAN_SNARL:
        // Stack has at the top running distances to use for each chain still
        // to be visited in the snarl, and under those the same for the snarl
        // above that, etc.
        switch (it->type) {
        case SNARL_START:
            // Stack holds running distance along parent chain plus edge
            // distance to cross the snarl, or running distance out of chain we
            // started in plus distance to exit the snarl.
            //
            // This is the right running distance to use for the parent chain now.
            // So go back to scanning the parent chain.
            state(S_SCAN_CHAIN);
            break;
        case CHAIN_END:
            // We've encountered a chain to look at, and the running distance
            // into the chain is already on the stack.
            if (top() > distance_limit) {
                // Running distance is already too high so skip over the chain
                push(0);
                state(S_SKIP_CHAIN);
            } else {
                // Do the chain
                state(S_SCAN_CHAIN);
            }
            break;
        case EDGE:
            // We've found edge data in the snarl, but we already know the
            // running distances to everythign we will encounter, so we ignore
            // it.
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_SKIP_CHAIN:
        /// Stack has the nesting level of child snarls we are reading over
        /// until we get back to the level we want to skip past the chain
        /// start.
        /// Under that is the running distance along the chain being skipped.
        /// And under that it has the running distance for ther next thing in
        /// the snarl, which had better exist or we shouldn't be trying to skip
        /// the chain, we should have halted.
        switch (it->type) {
        case SEED:
            // We don't emit seeds until the chain is over
            return false;
            break;
        case SNARL_START:
            // We might now be able to match chain starts again
            top() -= 1;
            break;
        case SNARL_END:
            // We can't match chain starts until we leave the snarl
            top() += 1;
            break;
        case CHAIN_START:
            if (top() == 0) {
                // This is the start of the chain we were wanting to skip.
                pop();
                // We definitely should have entered the parent snarl of the chain, or we would have halted instead of trying to skip the rest of the chain.
                crash_unless(depth() > 1);
                // Discard the running distance along this chain, which no longer matters.
                pop();
                // Running distance for next chain, or running distance to cross the snarl, will be under it.
                state(S_SCAN_SNARL);
            }
            // Otherwise this is the start of a chain inside a child snarl we are skipping over and we ignore it.
            break;
        case CHAIN_END:
            // Ignore chain ends
            break;
        case EDGE:
            // Ignore edge values
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    default:
        throw std::domain_error("Unimplemented state " + std::to_string(current_state)); 
    }
    // Unless we yield something, we don't yield anything.
    return false;
}

auto ZipCodeTree::look_back(const iterator& from, size_t distance_limit) const -> reverse_iterator {
    return reverse_iterator(zip_code_tree.rbegin() + from.remaining_tree(), zip_code_tree.rend(), distance_limit);
}
auto ZipCodeTree::rend() const -> reverse_iterator {
    return reverse_iterator(zip_code_tree.rend(), zip_code_tree.rend(), 0);
}


}