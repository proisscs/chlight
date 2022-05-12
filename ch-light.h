#ifndef FILE_CHL_SEEN
#define FILE_CHL_SEEN

#include <iostream>
#include <vector>
#include <random>
#include "c.h"
#include "Timer.h"
#include "CHParser.h"
#include <fstream>
#include <limits>

using namespace std;

struct PQElementCH
{
    Rank rank;
    Distance distance;
    Level min_rank;
    bool forward;

    bool operator()(const PQElementCH &lhs, const PQElementCH &rhs)
    {
        return lhs.distance > rhs.distance;
    }
};

class CHLight
{
public:
    CHLight()
    {
    }

    void constructCHLight(string file, int kept_shortcuts_in_percentage = 0)
    {
        vector<EdgeType> edge_list;
        vector<NodeID> rankToNodeID;
        vector<Rank> nodeIDToRank;
        long shortcut_counter = 0;
        Timer timer;
        {
            cout << "Loading CH..." << endl;
            CHParser ch_fmi;
            ch_fmi.readFromFMIFile(file);
            cout << "Finished." << endl;
            timer.start();
            _numNodes = ch_fmi.nofNodes();
            nodeIDToRank.resize(_numNodes);
            rankToNodeID = ch_fmi.rankToNodeID;
            vector<pair<long, Rank>> osm_rank_pairs(_numNodes);
            for (Rank i = 0; i < _numNodes; i++)
            {
                NodeID nodeID = rankToNodeID.at(i);
                osm_rank_pairs[i] = make_pair(ch_fmi.nodeList[nodeID].osmID, i);
            }
            _osm_id_to_rank = getOSMIDToRankMapping(osm_rank_pairs);
            _min_ranks.resize(_numNodes);
            for (Rank i = 0; i < _numNodes; i++)
            {
                Rank rank = _osm_id_to_rank[i];
                NodeID nodeID = rankToNodeID.at(rank);
                nodeIDToRank.at(nodeID) = i;
                _osm_id_to_rank[i] = i;
                int level = ch_fmi.nodeIDToLevel[nodeID];
                _min_ranks[i] = level <= numeric_limits<Level>::max() ? level : numeric_limits<Level>::max();
            }

            edge_list = ch_fmi.getEdges();
            for (EdgeType &edge : edge_list)
            {
                edge.source = nodeIDToRank[edge.source];
                edge.target = nodeIDToRank[edge.target];
                if (edge.child_1 != c::NO_ENTRY)
                {
                    shortcut_counter++;
                }
            }
        }
        long num_kept_shortcuts = shortcut_counter / 100 * kept_shortcuts_in_percentage;
        vector<bool> keep_shortcut(edge_list.size(), false);
        vector<EdgeID> index_to_edgeid(edge_list.size());
        { // Calc max ranks
            cout << "Start calculating max ranks..." << endl;
            sort(edge_list.begin(), edge_list.end(),
                 [&](const EdgeType &A, const EdgeType &B) -> bool
                 {
                     Level rank1 = _min_ranks[A.source] < _min_ranks[A.target] ? _min_ranks[A.source] : _min_ranks[A.target];
                     Level rank2 = _min_ranks[B.source] < _min_ranks[B.target] ? _min_ranks[B.source] : _min_ranks[B.target];
                     return rank1 > rank2;
                 });
            for (EdgeID edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                index_to_edgeid[edge_list[edgeid].index] = edgeid;
            }
            vector<bool> is_shortcut_processed(edge_list.size(), false);
            shortcut_counter = 0;
            for (EdgeID edgeid = edge_list.size() - 1; edgeid >= 0 && shortcut_counter < num_kept_shortcuts; edgeid--)
            {
                if (edge_list[edgeid].child_1 != c::NO_ENTRY)
                {
                    shortcut_counter++;
                    keep_shortcut[edge_list[edgeid].index] = true;
                }
            }
            _max_ranks = _min_ranks;
            for (EdgeID edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                const EdgeType &edge = edge_list[edgeid];
                if (!is_shortcut_processed[edge.index] && !keep_shortcut[edge.index])
                {
                    is_shortcut_processed[edge.index] = true;
                    if (edge.child_1 != c::NO_ENTRY)
                    {
                        Level cur_level = _min_ranks[edge.source] < _min_ranks[edge.target] ? _min_ranks[edge.source] : _min_ranks[edge.target];
                        vector<EdgeID> edge_queue;
                        long edge_queue_index = 0;
                        if (!is_shortcut_processed.at(edge.child_1))
                        {
                            edge_queue.push_back(edge.child_1);
                            is_shortcut_processed[edge.child_1] = true;
                        }
                        if (!is_shortcut_processed.at(edge.child_2))
                        {
                            edge_queue.push_back(edge.child_2);
                            is_shortcut_processed[edge.child_2] = true;
                        }
                        while (edge_queue_index < edge_queue.size())
                        {
                            const EdgeType &new_edge = edge_list[index_to_edgeid[edge_queue[edge_queue_index]]];
                            if (new_edge.child_1 != c::NO_ENTRY && !keep_shortcut[new_edge.index])
                            {
                                if (!is_shortcut_processed[new_edge.child_1])
                                {
                                    edge_queue.push_back(new_edge.child_1);
                                    is_shortcut_processed[new_edge.child_1] = true;
                                }
                                if (!is_shortcut_processed[new_edge.child_2])
                                {
                                    edge_queue.push_back(new_edge.child_2);
                                    is_shortcut_processed[new_edge.child_2] = true;
                                }
                            }
                            else
                            {
                                if (_max_ranks[new_edge.source] < cur_level)
                                {
                                    _max_ranks[new_edge.source] = cur_level;
                                }
                                if (_max_ranks[new_edge.target] < cur_level)
                                {
                                    _max_ranks[new_edge.target] = cur_level;
                                }
                            }
                            edge_queue_index++;
                        }
                    }
                }
            }
        }

        { // Save edges
            _fwd_edges.clear();
            _bwd_edges.clear();
            _fwd_edges.resize(_numNodes);
            _bwd_edges.resize(_numNodes);
            if (num_kept_shortcuts > 0)
            {
                _fwd_shortcuts.resize(_numNodes);
                _bwd_shortcuts.resize(_numNodes);
            }
            _numEdges = 0;
            _numShortcuts = 0;
            int count_unnecessary_edges = 0;

            sort(edge_list.begin(), edge_list.end(), compareEdgeTypes);
            for (EdgeID edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                index_to_edgeid[edge_list[edgeid].index] = edgeid;
            }
            cout << "Start initializing edges..." << endl;
            long current_rank = 0;
            for (long edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                const EdgeType &edge = edge_list[edgeid];
                if ((edge.child_1 == c::NO_ENTRY || keep_shortcut[edge.index]) && _min_ranks[edge.source] <= _max_ranks[edge.target] && edge.source != edge.target)
                {
                    if (edgeid == 0 || edge_list[edgeid - 1].source != edge.source || edge_list[edgeid - 1].target != edge.target)
                    {
                        if (edge.child_1 == c::NO_ENTRY)
                        {
                            _fwd_edges[edge.source].push_back(make_pair(edge.target, edge.weight));
                            _numEdges++;
                        }
                        else
                        {
                            _fwd_shortcuts[edge.source].push_back(Shortcut{
                                edge.target,
                                edge.weight,
                                edge_list.at(index_to_edgeid[edge.child_1]).target,
                            });
                            _numShortcuts++;
                        }
                    }
                    else
                    {
                        count_unnecessary_edges++;
                    }
                }
            }
            sort(edge_list.begin(), edge_list.end(), compareEdgeTypesBackwards);
            for (EdgeID edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                index_to_edgeid[edge_list[edgeid].index] = edgeid;
            }
            current_rank = 0;
            for (long edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                const EdgeType &edge = edge_list[edgeid];
                if ((edge.child_1 == c::NO_ENTRY || keep_shortcut[edge.index]) && _min_ranks[edge.target] <= _max_ranks[edge.source] && edge.source != edge.target)
                {
                    if (edgeid == 0 || edge_list[edgeid - 1].source != edge.source || edge_list[edgeid - 1].target != edge.target)
                    {
                        if (edge.child_1 == c::NO_ENTRY)
                        {
                            _bwd_edges[edge.target].push_back(make_pair(edge.source, edge.weight));
                            _numEdges++;
                        }
                        else
                        {
                            _bwd_shortcuts[edge.target].push_back(Shortcut{
                                edge.source,
                                edge.weight,
                                edge_list.at(index_to_edgeid[edge.child_1]).target,
                            });
                            _numShortcuts++;
                        }
                    }
                    else
                    {
                        count_unnecessary_edges++;
                    }
                }
            }

            cout << "Construction finished.\nNumber of edges: " << _numEdges << "\nShortcuts: " << _numShortcuts << "\nUnnecessary edges: " << count_unnecessary_edges << endl;
        }
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        _added_by_bwd.resize(_numNodes);
        _added_by_fwd.resize(_numNodes);
        _is_initialized = true;
        timer.stop();
        cout << "Construction took " << timer.secs() << " seconds." << endl;
        // DEBUG
        cout << "num edges: " << _numEdges << endl;
        // DEBUG END
    }

    bool isInitialized()
    {
        return _is_initialized;
    }

    Distance getDistance(Rank source, Rank target)
    {
        if (source == target)
        {
            return 0;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        fill(_bwdDistances.begin(), _bwdDistances.end(), c::NO_ENTRY);
        Distance best_dist = c::NO_ENTRY;
        priority_queue<PQElementCH, vector<PQElementCH>, PQElementCH> pQ;
        pQ.push(PQElementCH{
            source,
            0,
            _min_ranks[source],
            true,
        });
        pQ.push(PQElementCH{
            target,
            0,
            _min_ranks[target],
            false,
        });
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        while (!pQ.empty())
        {
            PQElementCH next_element = pQ.top();
            Rank rank = next_element.rank;
            Distance distance = next_element.distance;
            if (best_dist != c::NO_ENTRY && distance >= best_dist)
            {
                break;
            }
            Level min_rank = next_element.min_rank;
            pQ.pop();
            if (next_element.forward)
            {
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _bwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                        }
                    }
                    for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _fwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                        {
                            _fwdDistances[edge.first] = new_distance;
                            if (_max_ranks[edge.first] >= min_rank)
                            {
                                pQ.push(PQElementCH{
                                    edge.first,
                                    new_distance,
                                    min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                    true,
                                });
                            }
                        }
                    }
                }
            }
            else
            {
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _fwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                        }
                    }
                    for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _bwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_bwdDistances[edge.first] == c::NO_ENTRY || new_distance < _bwdDistances[edge.first])
                        {
                            _bwdDistances[edge.first] = new_distance;
                            if (_max_ranks[edge.first] >= min_rank)
                            {
                                pQ.push(PQElementCH{
                                    edge.first,
                                    new_distance,
                                    min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                    false,
                                });
                            }
                        }
                    }
                }
            }
        }
        return best_dist;
    }

    long getNumberVisitedNodes(Rank source, Rank target)
    {
        long counter_visited_nodes = 0;
        bool has_shortcuts = _fwd_shortcuts.size() == _numNodes;
        if (source == target)
        {
            return 0;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        fill(_bwdDistances.begin(), _bwdDistances.end(), c::NO_ENTRY);
        Distance best_dist = c::NO_ENTRY;
        priority_queue<PQElementCH, vector<PQElementCH>, PQElementCH> pQ;
        pQ.push(PQElementCH{
            source,
            0,
            _min_ranks[source],
            true,
        });
        pQ.push(PQElementCH{
            target,
            0,
            _min_ranks[target],
            false,
        });
        counter_visited_nodes += 2;
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        while (!pQ.empty())
        {
            PQElementCH next_element = pQ.top();
            Rank rank = next_element.rank;
            Distance distance = next_element.distance;
            if (best_dist != c::NO_ENTRY && distance >= best_dist)
            {
                break;
            }
            Level min_rank = next_element.min_rank;
            pQ.pop();
            if (next_element.forward)
            {
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _bwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                        }
                    }
                    for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _fwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                        {
                            if (_max_ranks[edge.first] >= min_rank)
                            {
                                pQ.push(PQElementCH{
                                    edge.first,
                                    new_distance,
                                    min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                    true,
                                });
                                if (_fwdDistances[edge.first] == c::NO_ENTRY)
                                {
                                    counter_visited_nodes++;
                                }
                            }
                            _fwdDistances[edge.first] = new_distance;
                        }
                    }
                    if (has_shortcuts)
                    {
                        for (int edgeid = 0; edgeid < _fwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _fwd_shortcuts[rank][edgeid];
                            Distance new_distance = edge.distance + distance;
                            if (_fwdDistances[edge.node] == c::NO_ENTRY || new_distance < _fwdDistances[edge.node])
                            {
                                if (_max_ranks[edge.node] >= min_rank)
                                {
                                    pQ.push(PQElementCH{
                                        edge.node,
                                        new_distance,
                                        min_rank < _min_ranks[edge.node] ? _min_ranks[edge.node] : min_rank,
                                        true,
                                    });
                                    if (_fwdDistances[edge.node] == c::NO_ENTRY)
                                    {
                                        counter_visited_nodes++;
                                    }
                                }
                                _fwdDistances[edge.node] = new_distance;
                            }
                        }
                    }
                }
            }
            else
            {
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _fwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                        }
                    }
                    for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _bwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_bwdDistances[edge.first] == c::NO_ENTRY || new_distance < _bwdDistances[edge.first])
                        {
                            if (_max_ranks[edge.first] >= min_rank)
                            {
                                pQ.push(PQElementCH{
                                    edge.first,
                                    new_distance,
                                    min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                    false,
                                });
                                if (_bwdDistances[edge.first] == c::NO_ENTRY)
                                {
                                    counter_visited_nodes++;
                                }
                            }
                            _bwdDistances[edge.first] = new_distance;
                        }
                    }
                    if (has_shortcuts)
                    {
                        for (int edgeid = 0; edgeid < _bwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _bwd_shortcuts[rank][edgeid];
                            Distance new_distance = edge.distance + distance;
                            if (_bwdDistances[edge.node] == c::NO_ENTRY || new_distance < _bwdDistances[edge.node])
                            {
                                if (_max_ranks[edge.node] >= min_rank)
                                {
                                    pQ.push(PQElementCH{
                                        edge.node,
                                        new_distance,
                                        min_rank < _min_ranks[edge.node] ? _min_ranks[edge.node] : min_rank,
                                        false,
                                    });
                                    if (_bwdDistances[edge.node] == c::NO_ENTRY)
                                    {
                                        counter_visited_nodes++;
                                    }
                                }
                                _bwdDistances[edge.node] = new_distance;
                            }
                        }
                    }
                }
            }
        }
        return counter_visited_nodes;
    }

    vector<Edge> getShortestPath(Rank source, Rank target)
    {
        vector<Edge> shortest_path;
        if (source == target)
        {
            return shortest_path;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        fill(_bwdDistances.begin(), _bwdDistances.end(), c::NO_ENTRY);
        Distance best_dist = c::NO_ENTRY;
        Rank mid_point = c::NO_ENTRY;
        priority_queue<PQElementCH, vector<PQElementCH>, PQElementCH> pQ;
        pQ.push(PQElementCH{
            source,
            0,
            _min_ranks[source],
            true,
        });
        pQ.push(PQElementCH{
            target,
            0,
            _min_ranks[target],
            false,
        });
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        while (!pQ.empty())
        {
            _counter++;
            PQElementCH next_element = pQ.top();
            Rank rank = next_element.rank;
            Distance distance = next_element.distance;
            if (best_dist != c::NO_ENTRY && distance > best_dist)
            {
                break;
            }
            Level min_rank = next_element.min_rank;
            pQ.pop();
            if (next_element.forward)
            {
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _bwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist || (new_dist == best_dist && _min_ranks[rank] > _min_ranks[mid_point]))
                        {
                            best_dist = new_dist;
                            mid_point = rank;
                        }
                    }
                    for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _fwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                        {
                            _fwdDistances[edge.first] = new_distance;
                            if (_max_ranks[edge.first] >= min_rank)
                            {
                                _added_by_fwd[edge.first] = rank;
                                pQ.push(PQElementCH{
                                    edge.first,
                                    new_distance,
                                    min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                    true,
                                });
                            }
                        }
                    }
                }
            }
            else
            {
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _fwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist || (new_dist == best_dist && _min_ranks[rank] > _min_ranks[mid_point]))
                        {
                            best_dist = new_dist;
                            mid_point = rank;
                        }
                    }
                    for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _bwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_bwdDistances[edge.first] == c::NO_ENTRY || new_distance < _bwdDistances[edge.first])
                        {
                            _bwdDistances[edge.first] = new_distance;
                            if (_max_ranks[edge.first] >= min_rank)
                            {
                                _added_by_bwd[edge.first] = rank;
                                pQ.push(PQElementCH{
                                    edge.first,
                                    new_distance,
                                    min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                    false,
                                });
                            }
                        }
                    }
                }
            }
        }
        { // Construct shortest path
            if (mid_point != c::NO_ENTRY)
            {
                vector<Edge> fwd_shortest_path;
                vector<Edge> bwd_shortest_path;
                fwd_shortest_path.reserve(_numNodes / 1000);
                bwd_shortest_path.reserve(_numNodes / 1000);
                Rank cur_rank = mid_point;
                while (cur_rank != source)
                {
                    Rank prev_rank = _added_by_fwd[cur_rank];
                    for (const Edge &edge : _fwd_edges[prev_rank])
                    {
                        if (edge.first == cur_rank)
                        {
                            fwd_shortest_path.push_back(edge);
                            break;
                        }
                    }
                    cur_rank = prev_rank;
                }
                cur_rank = mid_point;
                while (cur_rank != target)
                {
                    Rank prev_rank = _added_by_bwd[cur_rank];
                    for (const Edge &edge : _bwd_edges[prev_rank])
                    {
                        if (edge.first == cur_rank)
                        {
                            bwd_shortest_path.push_back(make_pair(prev_rank, edge.second));
                            break;
                        }
                    }
                    cur_rank = prev_rank;
                }
                shortest_path.reserve(fwd_shortest_path.size() + bwd_shortest_path.size());
                for (int edgeid = fwd_shortest_path.size() - 1; edgeid >= 0; edgeid--)
                {
                    shortest_path.push_back(fwd_shortest_path[edgeid]);
                }
                shortest_path.insert(shortest_path.end(), bwd_shortest_path.begin(), bwd_shortest_path.end());
            }
        }

        return shortest_path;
    }

    vector<Edge> getShortestPathWithSoD(Rank source, Rank target)
    {
        vector<Edge> shortest_path;
        if (source == target)
        {
            return shortest_path;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        fill(_bwdDistances.begin(), _bwdDistances.end(), c::NO_ENTRY);
        Distance best_dist = c::NO_ENTRY;
        Rank mid_point = c::NO_ENTRY;
        priority_queue<PQElementCH, vector<PQElementCH>, PQElementCH> pQ;
        pQ.push(PQElementCH{
            source,
            0,
            _min_ranks[source],
            true,
        });
        pQ.push(PQElementCH{
            target,
            0,
            _min_ranks[target],
            false,
        });
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        while (!pQ.empty())
        {
            _counter++;
            PQElementCH next_element = pQ.top();
            Rank rank = next_element.rank;
            Distance distance = next_element.distance;
            if (best_dist != c::NO_ENTRY && distance > best_dist)
            {
                break;
            }
            Level min_rank = next_element.min_rank;
            pQ.pop();
            if (next_element.forward)
            {
                if (_fwdDistances[rank] == distance)
                {
                    bool should_be_stalled = false;
                    for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _bwd_edges[rank][edgeid];
                        if (_fwdDistances[edge.first] != c::NO_ENTRY)
                        {
                            Distance new_distance = _fwdDistances[edge.first] + edge.second;
                            if (new_distance < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }
                    if (!should_be_stalled)
                    {
                        if (_bwdDistances[rank] != c::NO_ENTRY)
                        {
                            Distance new_dist = _bwdDistances[rank] + distance;
                            if (best_dist == c::NO_ENTRY || new_dist < best_dist || (new_dist == best_dist && _min_ranks[rank] > _min_ranks[mid_point]))
                            {
                                best_dist = new_dist;
                                mid_point = rank;
                            }
                        }
                        for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _fwd_edges[rank][edgeid];
                            Distance new_distance = edge.second + distance;
                            if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                            {
                                _fwdDistances[edge.first] = new_distance;
                                if (_max_ranks[edge.first] >= min_rank)
                                {
                                    _added_by_fwd[edge.first] = rank;
                                    pQ.push(PQElementCH{
                                        edge.first,
                                        new_distance,
                                        min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                        true,
                                    });
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                if (_bwdDistances[rank] == distance)
                {
                    bool should_be_stalled = false;
                    for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _fwd_edges[rank][edgeid];
                        if (_bwdDistances[edge.first] != c::NO_ENTRY)
                        {
                            Distance new_distance = _bwdDistances[edge.first] + edge.second;
                            if (new_distance < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }
                    if (!should_be_stalled)
                    {
                        if (_fwdDistances[rank] != c::NO_ENTRY)
                        {
                            Distance new_dist = _fwdDistances[rank] + distance;
                            if (best_dist == c::NO_ENTRY || new_dist < best_dist || (new_dist == best_dist && _min_ranks[rank] > _min_ranks[mid_point]))
                            {
                                best_dist = new_dist;
                                mid_point = rank;
                            }
                        }
                        for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _bwd_edges[rank][edgeid];
                            Distance new_distance = edge.second + distance;
                            if (_bwdDistances[edge.first] == c::NO_ENTRY || new_distance < _bwdDistances[edge.first])
                            {
                                _bwdDistances[edge.first] = new_distance;
                                if (_max_ranks[edge.first] >= min_rank)
                                {
                                    _added_by_bwd[edge.first] = rank;
                                    pQ.push(PQElementCH{
                                        edge.first,
                                        new_distance,
                                        min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                        false,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
        { // Construct shortest path
            if (mid_point != c::NO_ENTRY)
            {
                vector<Edge> fwd_shortest_path;
                vector<Edge> bwd_shortest_path;
                fwd_shortest_path.reserve(_numNodes / 1000);
                bwd_shortest_path.reserve(_numNodes / 1000);
                Rank cur_rank = mid_point;
                while (cur_rank != source)
                {
                    Rank prev_rank = _added_by_fwd[cur_rank];
                    for (const Edge &edge : _fwd_edges[prev_rank])
                    {
                        if (edge.first == cur_rank)
                        {
                            fwd_shortest_path.push_back(edge);
                            break;
                        }
                    }
                    cur_rank = prev_rank;
                }
                cur_rank = mid_point;
                while (cur_rank != target)
                {
                    Rank prev_rank = _added_by_bwd[cur_rank];
                    for (const Edge &edge : _bwd_edges[prev_rank])
                    {
                        if (edge.first == cur_rank)
                        {
                            bwd_shortest_path.push_back(make_pair(prev_rank, edge.second));
                            break;
                        }
                    }
                    cur_rank = prev_rank;
                }
                shortest_path.reserve(fwd_shortest_path.size() + bwd_shortest_path.size());
                for (int edgeid = fwd_shortest_path.size() - 1; edgeid >= 0; edgeid--)
                {
                    shortest_path.push_back(fwd_shortest_path[edgeid]);
                }
                shortest_path.insert(shortest_path.end(), bwd_shortest_path.begin(), bwd_shortest_path.end());
            }
        }

        return shortest_path;
    }

    vector<Edge> getShortestPathWithShortcuts(Rank source, Rank target)
    {
        vector<Edge> shortest_path;
        if (source == target)
        {
            return shortest_path;
        }
        vector<Rank> fwd_visited_nodes;
        // fwd_visited_nodes.reserve(_numNodes / 100);
        priority_queue<PQElementCH, vector<PQElementCH>, PQElementCH> pQ;
        pQ.push(PQElementCH{
            source,
            0,
            _min_ranks[source],
            true,
        });
        pQ.push(PQElementCH{
            target,
            0,
            _min_ranks[target],
            false,
        });
        fwd_visited_nodes.push_back(source);
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        vector<Rank> bwd_visited_nodes;
        // bwd_visited_nodes.reserve(_numNodes / 100);
        bwd_visited_nodes.push_back(target);
        Distance best_dist = c::NO_ENTRY;
        Rank midpoint = c::NO_ENTRY;
        bool forward = true;
        while (!pQ.empty())
        {
            _counter++;
            PQElementCH next_element = pQ.top();
            pQ.pop();
            forward = next_element.forward;
            Rank rank = next_element.rank;
            Level min_rank = next_element.min_rank;
            Distance distance = next_element.distance;
            if (best_dist != c::NO_ENTRY && distance > best_dist)
            {
                break;
            }
            if (forward)
            {
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = distance + _bwdDistances[rank];
                        if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                        {
                            best_dist = new_dist;
                            midpoint = rank;
                        }
                    }
                    bool should_be_stalled = false;
                    /*for (const Shortcut &edge : _bwd_shortcuts[rank])
                    {
                        if (_fwdDistances[edge.node] != c::NO_ENTRY && _fwdDistances[edge.node] + edge.distance < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (const Edge &edge : _bwd_edges[rank])
                        {
                            if (_fwdDistances[edge.first] != c::NO_ENTRY && _fwdDistances[edge.first] + edge.second < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }*/
                    if (!should_be_stalled || true)
                    {
                        for (int edgeid = 0; edgeid < _fwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _fwd_shortcuts[rank][edgeid];
                            Distance new_dist = distance + edge.distance;
                            if (_fwdDistances[edge.node] == c::NO_ENTRY || _fwdDistances[edge.node] > new_dist)
                            {
                                /*if (_fwdDistances[edge.node] == c::NO_ENTRY)
                                {
                                    fwd_visited_nodes.push_back(edge.node);
                                }*/
                                _fwdDistances[edge.node] = new_dist;
                                pQ.push(PQElementCH{
                                    edge.node,
                                    new_dist,
                                    min_rank < _min_ranks[edge.node] ? _min_ranks[edge.node] : min_rank,
                                    true,
                                });
                                _added_by_fwd[edge.node] = rank;
                            }
                        }
                        for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _fwd_edges[rank][edgeid];
                            Distance new_dist = distance + edge.second;
                            if (_fwdDistances[edge.first] == c::NO_ENTRY || _fwdDistances[edge.first] > new_dist)
                            {
                                /*if (_fwdDistances[edge.first] == c::NO_ENTRY)
                                {
                                    fwd_visited_nodes.push_back(edge.first);
                                }*/
                                _fwdDistances[edge.first] = new_dist;
                                pQ.push(PQElementCH{
                                    edge.first,
                                    new_dist,
                                    min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                    true,
                                });
                                _added_by_fwd[edge.first] = rank;
                            }
                        }
                    }
                }
            }
            else
            {
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = distance + _fwdDistances[rank];
                        if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                        {
                            best_dist = new_dist;
                            midpoint = rank;
                        }
                    }
                    bool should_be_stalled = false;
                    /*for (const Shortcut &edge : _fwd_shortcuts[rank])
                    {
                        if (_bwdDistances[edge.node] != c::NO_ENTRY && _bwdDistances[edge.node] + edge.distance < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (const Edge &edge : _fwd_edges[rank])
                        {
                            if (_bwdDistances[edge.first] != c::NO_ENTRY && _bwdDistances[edge.first] + edge.second < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }*/
                    if (!should_be_stalled || true)
                    {
                        for (int edgeid = 0; edgeid < _bwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _bwd_shortcuts[rank][edgeid];
                            if (min_rank <= _max_ranks[edge.node])
                            {
                                Distance new_dist = distance + edge.distance;
                                if (_bwdDistances[edge.node] == c::NO_ENTRY || _bwdDistances[edge.node] > new_dist)
                                {
                                    /*if (_bwdDistances[edge.node] == c::NO_ENTRY)
                                    {
                                        bwd_visited_nodes.push_back(edge.node);
                                    }*/
                                    _bwdDistances[edge.node] = new_dist;
                                    pQ.push(PQElementCH{
                                        edge.node,
                                        new_dist,
                                        min_rank < _min_ranks[edge.node] ? _min_ranks[edge.node] : min_rank,
                                        false,
                                    });
                                    _added_by_bwd[edge.node] = rank;
                                }
                            }
                        }
                        for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _bwd_edges[rank][edgeid];
                            if (min_rank <= _max_ranks[edge.first])
                            {
                                Distance new_dist = distance + edge.second;
                                if (_bwdDistances[edge.first] == c::NO_ENTRY || _bwdDistances[edge.first] > new_dist)
                                {
                                    /*if (_bwdDistances[edge.first] == c::NO_ENTRY)
                                    {
                                        bwd_visited_nodes.push_back(edge.first);
                                    }*/
                                    _bwdDistances[edge.first] = new_dist;
                                    pQ.push(PQElementCH{
                                        edge.first,
                                        new_dist,
                                        min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                        false,
                                    });
                                    _added_by_bwd[edge.first] = rank;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (fwd_visited_nodes.size() > 0)
        {
            fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        }
        else
        {
            for (auto rank : fwd_visited_nodes)
            {
                _fwdDistances[rank] = c::NO_ENTRY;
            }
        }
        if (bwd_visited_nodes.size() > 0)
        {
            fill(_bwdDistances.begin(), _bwdDistances.end(), c::NO_ENTRY);
        }
        else
        {
            for (auto rank : bwd_visited_nodes)
            {
                _bwdDistances[rank] = c::NO_ENTRY;
            }
        }
        if (midpoint != c::NO_ENTRY)
        {
            struct TmpEdge
            {
                Rank source;
                Rank target;
            };
            shortest_path.reserve(10000);
            vector<TmpEdge> fwd_ch_path;
            fwd_ch_path.reserve(10000);
            Rank cur_node = midpoint;
            while (cur_node != source)
            {
                fwd_ch_path.push_back(TmpEdge{
                    _added_by_fwd[cur_node],
                    cur_node,
                });
                cur_node = _added_by_fwd[cur_node];
            }
            vector<TmpEdge> bwd_ch_path;
            bwd_ch_path.reserve(10000);
            cur_node = midpoint;
            while (cur_node != target)
            {
                bwd_ch_path.push_back(TmpEdge{
                    cur_node,
                    _added_by_bwd[cur_node],
                });
                cur_node = _added_by_bwd[cur_node];
            }
            vector<TmpEdge> queue;
            queue.reserve(10000);
            for (int i = bwd_ch_path.size() - 1; i >= 0; i--)
            {
                queue.push_back(bwd_ch_path[i]);
            }
            for (int i = 0; i < fwd_ch_path.size(); i++)
            {
                queue.push_back(fwd_ch_path[i]);
            }
            int queue_index = queue.size();
            while (queue_index > 0)
            {
                queue_index--;
                const TmpEdge cur_edge = queue[queue_index];
                bool is_forward = cur_edge.source < cur_edge.target;
                const vector<Shortcut> &shortcuts = is_forward ? _fwd_shortcuts[cur_edge.source] : _bwd_shortcuts[cur_edge.target];
                bool edge_found = false;
                Rank other_side_node = is_forward ? cur_edge.target : cur_edge.source;
                for (const Shortcut &edge : shortcuts)
                {
                    if (edge.node == other_side_node)
                    {
                        edge_found = true;
                        queue[queue_index] = TmpEdge{
                            edge.midnode,
                            cur_edge.target,
                        };
                        queue_index++;
                        if (queue_index == queue.size())
                        {
                            queue.push_back(TmpEdge{
                                cur_edge.source,
                                edge.midnode,
                            });
                        }
                        else
                        {
                            queue[queue_index] = TmpEdge{
                                cur_edge.source,
                                edge.midnode,
                            };
                        }
                        queue_index++;
                        break;
                    }
                }
                if (!edge_found)
                {
                    const vector<Edge> &edges = is_forward ? _fwd_edges[cur_edge.source] : _bwd_edges[cur_edge.target];
                    for (const Edge &edge : edges)
                    {
                        if (edge.first == other_side_node)
                        {
                            edge_found = true;
                            shortest_path.push_back(make_pair(cur_edge.target, edge.second));
                            break;
                        }
                    }
                }
            }
        }
        return shortest_path;
    }

    vector<pair<Rank, Rank>> getUpwardSearch(Rank source)
    {
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        priority_queue<PQElementCH, vector<PQElementCH>, PQElementCH> pQ;
        pQ.push(PQElementCH{
            source,
            0,
            _min_ranks[source],
            true,
        });
        _fwdDistances[source] = 0;
        _added_by_fwd[source] = source;
        while (!pQ.empty())
        {
            _counter++;
            PQElementCH next_element = pQ.top();
            Rank rank = next_element.rank;
            Distance distance = next_element.distance;
            Level min_rank = next_element.min_rank;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                {
                    const Edge &edge = _fwd_edges[rank][edgeid];
                    Distance new_distance = edge.second + distance;
                    if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                    {
                        _fwdDistances[edge.first] = new_distance;
                        if (_max_ranks[edge.first] >= min_rank)
                        {
                            _added_by_fwd[edge.first] = rank;
                            pQ.push(PQElementCH{
                                edge.first,
                                new_distance,
                                min_rank < _min_ranks[edge.first] ? _min_ranks[edge.first] : min_rank,
                                true,
                            });
                        }
                    }
                }
            }
        }
        vector<pair<Rank, Rank>> upward_search;
        for (Rank rank = 0; rank < _numNodes; rank++)
        {
            if (_fwdDistances[rank] != c::NO_ENTRY)
            {
                upward_search.push_back(make_pair(rank, _added_by_fwd[rank]));
            }
        }
        return upward_search;
    }

    long nofNodes()
    {
        return _numNodes;
    }

    long nofEdges()
    {
        return _numEdges;
    }

    long getSizeOfCH()
    {
        if (_numShortcuts == 0)
        {
            return _numEdges * sizeof(Edge) + 2 * sizeof(Rank) * _numNodes + 2 * sizeof(Level) * _numNodes;
        }
        return _numEdges * sizeof(Edge) + _numShortcuts * sizeof(Shortcut) + 4 * sizeof(Rank) * _numNodes + 2 * sizeof(Level) * _numNodes;
    }

    void resetCounter()
    {
        _counter = 0;
    }

    long getCounter()
    {
        return _counter;
    }

    Rank osmIDToRank(long osmID)
    {
        return _osm_id_to_rank[osmID];
    }

private:
    long _numNodes;
    long _numEdges;
    long _numShortcuts;
    long _counter = 0;
    bool _is_initialized = false;
    vector<vector<Edge>> _fwd_edges;
    vector<vector<Edge>> _bwd_edges;
    vector<vector<Shortcut>> _fwd_shortcuts;
    vector<vector<Shortcut>> _bwd_shortcuts;
    vector<Level> _max_ranks;
    vector<Level> _min_ranks;

    vector<Distance> _fwdDistances;
    vector<Distance> _bwdDistances;
    vector<Rank> _added_by_fwd;
    vector<Rank> _added_by_bwd;
    vector<Rank> _osm_id_to_rank;
};

#endif